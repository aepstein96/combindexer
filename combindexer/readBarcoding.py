from collections import defaultdict, deque
from contextlib import ExitStack
from multiprocessing import Pool
from functools import partial
import gzip
import glob
import sys
import os
import fire
import csv
from string import Formatter


# For a given barcode, returns list of barcodes a given Levenshtein distance away. This is essential for efficient barcode correction.
def getCombinedVariants(
    seq, 
    distance,
    post_context="",
    correct_snps=True,
    correct_indels=True,
    bases={'A', 'C', 'G', 'T', 'N'}
):
    """
    Generates all unique, fixed-length variants of a sequence within a given
    combined correction distance.

    This version has been corrected to remove the flawed pre_context logic.
    All deletions are now correctly compensated by the post_context.
    
    Returns a map of {variant: (post_context, distance)}.
    """
    # The state now tracks variant, post_context, and distance
    all_variants_map = {seq: (post_context, 0)}
    last_generation = [(seq, post_context)]

    for d in range(distance):
        current_generation = []
        
        for var, post in last_generation:
            
            # 1. Generate SNPs (these don't change post_context)
            if correct_snps:
                for i in range(len(var)):
                    original_base = var[i]
                    for new_base in bases:
                        if new_base != original_base:
                            new_variant = var[:i] + new_base + var[i+1:]
                            if new_variant not in all_variants_map:
                                all_variants_map[new_variant] = (post, d + 1)
                                current_generation.append((new_variant, post))

            # 2. Generate Indels
            if correct_indels:
                # Deletions within the slice (compensated by post_context)
                for i in range(len(var)):
                    deleted_in_slice = var[:i] + var[i+1:]
                    if post:
                        new_variant = deleted_in_slice + post[0]
                        if new_variant not in all_variants_map:
                            # New post_context is one char shorter
                            new_post = post[1:]
                            all_variants_map[new_variant] = (new_post, d + 1)
                            current_generation.append((new_variant, new_post))
                    else: # Fallback: No post_context, pad with all bases
                        for b in bases:
                            new_variant = deleted_in_slice + b
                            if new_variant not in all_variants_map:
                                all_variants_map[new_variant] = ("", d + 1)
                                current_generation.append((new_variant, ""))
                
                # Insertions anywhere (compensated by truncation)
                for i in range(len(var) + 1):
                    for b in bases:
                        inserted = var[:i] + b + var[i:]
                        new_variant = inserted[:len(seq)]
                        if new_variant not in all_variants_map:
                            # New post_context is the truncated char + old post
                            new_post = inserted[-1] + post
                            all_variants_map[new_variant] = (new_post, d + 1)
                            current_generation.append((new_variant, new_post))
        
        last_generation = current_generation
        
    return all_variants_map


# For a list of correct barcodes, creates a dictionary mapping uncorrected barcodes to correct barcodes
def getCorrDict(correct_bars, default=None, ambiguous=None, **kwargs):
    out_dict = defaultdict(lambda: default) # returns default if a non-matching barcode is inputted
    ambiguous_variants = set()
    for correct_bar in correct_bars:
        for bar in getCombinedVariants(correct_bar, **kwargs):
            if bar not in ambiguous_variants:
                if bar in out_dict:
                    if out_dict[bar] != correct_bar:
                        ambiguous_variants.add(bar)
                        if ambiguous:
                            out_dict[bar] = ambiguous
                        else:
                            del out_dict[bar]
                else:
                    out_dict[bar] = correct_bar

    return out_dict

def findCorrectBarcode(var_seq, candidates):
    if len(candidates) == 1:
        return candidates[0][0]

    min_dist = min(c[1] for c in candidates)
    best_matches = [c[0] for c in candidates if c[1] == min_dist]

    if len(best_matches) == 1:
        return best_matches[0]
    
    return False

# Returns reverse complement of sequence
# Not as flexible as Bio.Seq, but probably faster and easier to use
def revComp(seq, comp_dict={'A':'T','G':'C','C':'G','T':'A','N':'N'}):
    return ''.join([comp_dict[base] for base in seq[::-1]])

# Handles a single input read file, including the barcode positions and UMIs contained within it
class InputRead:
    def __init__(self, name, sample, commands, file_pattern, trim_partner=None, trim_partner_len=17, min_length=20, debug_mode=False):
        self.name = name
        self.sample = sample
        self.commands = commands
        self.file_pattern = file_pattern
        self.trim_partner = trim_partner # should be in order of position
        self.trim_partner_len = trim_partner_len
        self.min_length = min_length
        
        if debug_mode:
            print("Input pattern: %s. Sample: %s" % (self.file_pattern, sample))
            
        self.file_paths = glob.glob(self.file_pattern.replace("{SAMPLE}", sample))
        if self.file_paths:
            self.file_path = self.file_paths[0]
        else:
            raise ValueError("%s not found" % self.file_pattern.replace("{SAMPLE}", sample))
    
    def __enter__(self):
        if self.file_path.endswith('.gz'):
            self.input = gzip.open(self.file_path, 'rt') #rt makes it output text and not binary
        else:
            self.input = open(self.file_path, 'rt')
        self.head = 'placeholder' # required to start the while loop
        
    
    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.input.close()
    
    # Moves on to the next read sequence in the input file
    def nextRead(self):
        # Read in data from file
        self.head = self.input.readline().rstrip() # remove newline marker, it will be put back later
        self.seq = self.input.readline().rstrip()
        self.sep = self.input.readline().rstrip()
        self.qual = self.input.readline().rstrip()
    
    # For a given read sequence, returns the barcodes and UMIs contained within it
    def getBarcodes(self, barcode_dict, debug_mode=False):
            
        start = 0
        cur_bars = {}
        UMI = ''
        UMI_qual = ''
        
        if not self.commands or not self.commands[0]: # this is what happens if that read has no barcodes in it
            self.bar_seq = ''
            return UMI, UMI_qual, cur_bars
        
        command_queue = deque(self.commands)
        trim_pos = 0
        while command_queue:
            command_name = command_queue.popleft()
            
            if command_name.startswith('Goto:'):
                start = int(command_name[5:])
                
            elif command_name.startswith('Gap:'):
                start += int(command_name[4:])
                
            elif command_name.startswith('UMI:'):
                UMI_len = int(command_name[4:])
                cur_UMI = self.seq[start:start+UMI_len]
                cur_UMI_qual = self.qual[start:start+UMI_len]
                if 'N' in cur_UMI:
                    return False # this is kind of inelegant, maybe there is a better way to do this
                
                UMI += cur_UMI
                UMI_qual += cur_UMI_qual
                start += UMI_len
                
            elif command_name == 'Trim':
                trim_pos = start    
    
            elif command_name.startswith('Bar:'): # a barcode
                bar_name = command_name[4:]
                bar_pos = barcode_dict[bar_name]
                potential_bar_seq = self.seq[start:start+bar_pos.length]
                bar = bar_pos.getBar(potential_bar_seq)
                
                if not bar:
                    if debug_mode:
                        print("No match for %s. Start: %d. File: %s. Read sequence: %s. Test sequence: %s" % (bar_name, start, self.file_path, self.seq, potential_bar_seq))
                    return False
                
                start += bar.length
                cur_bars[bar_name] = bar
                if bar.next_commands:
                    command_queue = deque(bar.next_commands)
            else:
                raise ValueError("Command not found: %s. Valid commands: Goto:<pos>, Gap:<len>, UMI:<len>, Trim, Bar:<name>" % command_name)
                
        self.bar_seq = self.seq[:start] # used for trimming partner
        if debug_mode:
            print("All barcodes found! Read: %s, Barcodes: %s\n" % (self.name, str([bar.seq for bar in cur_bars.values()])))
        self.seq = self.seq[trim_pos:]
        self.qual = self.qual[trim_pos:]
            
        return UMI, UMI_qual, cur_bars
        
    def trimRevComp(self, partner_seq):
        if len(partner_seq) < self.trim_partner_len: # Don't trim based on too short a sequence
            return False
        
        search_seq = revComp(partner_seq)[:self.trim_partner_len]
        search_seq_pos = self.seq.find(search_seq)
        if search_seq_pos != -1:
            self.seq = self.seq[:search_seq_pos]
            self.qual = self.qual[:search_seq_pos]
        
        return True


# General output file class, handles writing to either FASTQ or CSV files
class OutputFile:
    def __init__(self, name, sample, file_pattern, header_field_str, sequence_field_str, default_write=True, default_qual='I'):
        self.name = name
        self.sample = sample
        self.out_path = file_pattern.replace("{SAMPLE}", sample)
        self.out = None
        self.default_write = default_write
        self.default_qual = default_qual
        
        # Parse header and sequence field strings (assembly will be done during writing)
        self.header_field_str = header_field_str
        self.sequence_field_str = sequence_field_str
        
        # Determine output type based on file extension
        if self.out_path.endswith('.gz'):
            self.open_func = gzip.open
        else:
            self.open_func = open
        
        if any(self.out_path.lower().endswith(ext) for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']):
            self.output_type = 'fastq'
        else:
            self.output_type = 'csv'
            
    def __enter__(self):
        self.out = self.open_func(self.out_path, 'wt')
        if self.output_type == 'csv':
            self.out.write(self.header_field_str.replace('{', '').replace('}', '') + ',' + self.sequence_field_str.replace('{', '').replace('}', '') + '\n')
        return self
    
    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.out.close()
        
    def write(self, current_barcodes, current_reads, current_UMI = '', current_UMI_qual = None):
        if not self.out:
            raise ValueError("Can't write to output file because it's not open")
        
        # Assemble and parse header fields
        parsed = Formatter().parse(self.header_field_str)
        header_string = ''
        for literal, barcode_pos, field, _ in parsed: # slight modification to standard f-string format for my own purposes
            if literal:
                header_string += literal
            if barcode_pos == 'UMI':
                header_string += current_UMI
                
            elif field == 'PlateWell':
                plate = current_barcodes[barcode_pos].plate
                well = current_barcodes[barcode_pos].well
                if plate and well:
                    header_string += (plate + '-' + well)
                else:
                    header_string += (plate + well)
                
            elif field == 'Condition':
                header_string += current_barcodes[barcode_pos].condition
                
            elif field == 'ReadType':
                header_string += current_barcodes[barcode_pos].read_type
                
            elif field == 'Seq':
                header_string += current_barcodes[barcode_pos].seq
                
            else:
                raise ValueError("Invalid field: %s" % field)
            
            
        # Assemble and parse sequence fields
        parsed = Formatter().parse(self.sequence_field_str)
        sequence_string = ''
        qual_string = ''
        for literal, output_type, output_name, _ in parsed: # slight modification to standard f-string format for my own purposes
            if literal:
                sequence_string += literal
                qual_string += self.default_qual * len(literal)
                
            if output_type == 'UMI':
                sequence_string += current_UMI
                if current_UMI_qual is None:
                    current_UMI_qual = self.default_qual * len(current_UMI)
                qual_string += current_UMI_qual
            else:
                # Extract length if needed
                if '(' in output_name:
                    output_name, output_len = output_name.split('(')
                    output_len = int(output_len.rstrip(')'))
                else:
                    output_len = None # use entire equence
                    
                if output_type == 'Read':
                    sequence_string += current_reads[output_name].seq[:output_len]
                    qual_string += current_reads[output_name].qual[:output_len]
                elif output_type == 'Barcode':
                    sequence_string += current_barcodes[output_name].seq[:output_len]
                    
                    if output_len is None:
                        output_len = len(current_barcodes[output_name].seq)
                    qual_string += self.default_qual * output_len
                
        
        if self.output_type == 'fastq':
            self.out.write('@' + header_string + '\n')
            self.out.write(sequence_string + '\n')
            self.out.write('+\n')
            self.out.write(qual_string + '\n')
            
        elif self.output_type == 'csv':
            self.out.write(header_string + ',' + sequence_string + '\n')

# This is a single barcode, e.g. "CATGATA", corresponding to well A1, condition HeLa, read type oligo-dT, etc.
class Barcode:
    def __init__(self, seq='', plate='', well='', condition='', read_type='', next_commands=[], output_files=''):
        self.seq = seq # full sequence, including part that is not used
        self.matching_seq = seq # sequence that is used for matching; will later be truncated to the length used for matching
        self.plate = plate
        self.well = well
        self.condition = condition
        self.read_type = read_type
        self.next_commands = next_commands
        self.length = len(seq)
        
        # Parse output files as semicolon-separated set
        if output_files:
            self.output_targets = {target.strip() for target in output_files.split(';') if target.strip()}
        else:
            self.output_targets = set()
            

# This is a barcode position, e.g. "RT". It contains a list of barcodes, each of which is a Barcode object.
class BarcodePos:
    def __init__(self, bar_list, length_to_use, plate_list, corr_dist, orientation, correct_snps, correct_indels, debug_mode=False):
        
        # Get length of shortest barcode, this is what will be used for disambiguation
        # Also see if there is a default for "No match"
        self.no_match = None
        self.ambiguous = None
        min_length = length_to_use
        true_bars = [] # doesn't include no match
        for bar in bar_list:
            
            # Define no match and ambiguous barcodes
            if bar.seq == 'No match':
                self.no_match = bar
                continue
            if bar.seq == 'Ambiguous':
                self.ambiguous = bar
                continue

            # Define the shortest barcode length
            if bar.length < min_length:
                min_length = bar.length
                
            # Add the barcode to the list of true barcodes
            true_bars.append(bar)
        
        self.length = min_length
        self.plate_list = plate_list
        self.orientation = orientation

        # Make barcode correction dictionary
        potential_corrections = defaultdict(list)
        for bar in true_bars:
            if plate_list:
                if bar.plate not in plate_list:
                    continue
                    
            if orientation == 'reverse':
                seq = revComp(bar.seq)
            else:
                seq = bar.seq
            
            post_context = seq[min_length:]
            bar_seq = seq[:min_length]
            bar.matching_seq = bar_seq

            
            variants_map = getCombinedVariants(bar_seq, corr_dist, post_context, correct_snps, correct_indels)
            for var_seq, (_, dist) in variants_map.items():
                potential_corrections[var_seq].append((bar, dist))

        self.corr_dict = defaultdict(lambda: self.no_match)
        for var_seq, candidates in potential_corrections.items():
            correct_bar = findCorrectBarcode(var_seq, candidates)
            if correct_bar:
                self.corr_dict[var_seq] = correct_bar
            elif self.ambiguous:
                self.corr_dict[var_seq] = self.ambiguous
                        
    
    # Returns barcode information
    def getBar(self, bar_seq):
        return self.corr_dict[bar_seq]


# Make dictionary of barcode positions. Each barcode position is a BarcodePos object.
def makeBarDict(barcode_folder, debug_mode=False):
    bar_dict = {}
    for bar_file in os.listdir(barcode_folder):
        if bar_file.endswith('.csv') and bar_file.startswith('Barcode-'):
            name = bar_file.split('-')[1].split('.')[0]
            with open(os.path.join(barcode_folder, bar_file), 'rt') as f:
                length_to_use = f.readline().split(',')[1].strip()
                if length_to_use == 'all':
                    length_to_use = float('inf')
                else:
                    length_to_use = int(length_to_use)
                    
                plates_to_use = f.readline().split(',')[1].strip()
                if plates_to_use == 'all':
                    plate_list = None
                else:
                    plate_list = plates_to_use.split(';')
                    
                corr_dist = f.readline().split(',')[1].strip()
                if corr_dist: # deals with 0 becoming blank in certain cases
                    corr_dist = int(corr_dist)
                else:
                    corr_dist = 0
                
                correct_snps = f.readline().split(',')[1].strip()
                if correct_snps.lower() == 'true':
                    correct_snps = True
                else:
                    correct_snps = False 
                    
                correct_indels = f.readline().split(',')[1].strip()
                if correct_indels.lower() == 'true':
                    correct_indels = True
                else:
                    correct_indels = False
                    
                orientation = f.readline().split(',')[1].strip()
                f.readline()
                f.readline()
                
                bar_list = []
                for line in f:
                    parts = line.strip().split(',')
                    if len(parts) < 3:
                        continue
                    
                    seq = parts[0]
                    plate = parts[1]
                    well = parts[2]
                    condition = parts[3] if len(parts) > 3 else ''
                    read_type = parts[4] if len(parts) > 4 else ''
                    next_command_string = parts[5] if len(parts) > 5 else ''
                    output_files = parts[6] if len(parts) > 6 else ''
                    
                    if next_command_string:
                        next_commands = next_command_string.split(';')
                    else:
                        next_commands = None
                    bar_list.append(Barcode(seq, plate, well, condition, read_type, next_commands, output_files))                                                                                                        
            
            if debug_mode:
                print("%d barcodes found for %s" % (len(bar_list), name))
            bar_dict[name] = BarcodePos(bar_list, length_to_use, plate_list, corr_dist, orientation, correct_snps, correct_indels, debug_mode)
    return bar_dict

# Make dictionary of input reads. Each input read is an InputRead object.
def makeInputReadDict(sample, reads_file, input_folder, debug_mode=False):
    # Get input read information
    input_read_dict = {}
    with open(reads_file, 'rt') as f:
        f.readline() # skips header

        # Iterate through input reads
        for line in f:
            parts = line.strip().split(',')
            name = parts[0]
            file_pattern = os.path.join(input_folder, parts[1])
            command_str = parts[2]
            trim_partner_str = parts[3] if len(parts) > 3 else ""
            min_length = parts[4] if len(parts) > 4 else ""
                        
            commands = command_str.strip().split(';') if command_str else []
            if trim_partner_str:
                trim_partner, trim_partner_len = trim_partner_str.strip().split(':')
                trim_partner_len = int(trim_partner_len)
            else:
                trim_partner = None
                trim_partner_len = None
            
            if min_length:
                min_length = int(min_length)
            else:
                min_length = None  # Use default from InputRead constructor
                
            input_read = InputRead(name, sample, commands, file_pattern, trim_partner, trim_partner_len, min_length, debug_mode)
            input_read_dict[name] = input_read

    return input_read_dict

# Make dictionary of output files. Each output file is an OutputFile object.
def makeOutputFiles(sample, outputs_file, output_folder):
    output_files = {}
    default_output_files = []
    with open(outputs_file, 'rt') as f:
        f.readline() # skip header
        reader = csv.reader(f)
        for line in reader:
            name, file_pattern, default_write_setting, header_field_str, sequence_field_str = line[:5]
            if default_write_setting.lower() == 'yes':
                default_output_files.append(name)

            # Prepend output_folder to file_pattern
            full_file_pattern = os.path.join(output_folder, file_pattern)
            output_file = OutputFile(name, sample, full_file_pattern, header_field_str, sequence_field_str)
            output_files[name] = output_file
            
    return output_files, default_output_files

# Process a single sample, including all input reads and output files
def barcodeReadsSample(sample, barcode_folder, input_folder, output_folder, debug_mode=False):
    # Make input read dictionary
    reads_file = os.path.join(barcode_folder, "Reads.csv")
    input_read_dict = makeInputReadDict(sample, reads_file, input_folder, debug_mode)
            
    # Make output file dictionary
    outputs_file = os.path.join(barcode_folder, "Outputs.csv")
    output_file_dict, default_output_files = makeOutputFiles(sample, outputs_file, output_folder)
    
    # Make barcode position dictionary
    bar_pos_dict = makeBarDict(barcode_folder, debug_mode)                                                          
 
    # Get lists of input reads and output files
    input_reads = list(input_read_dict.values())
    output_files = list(output_file_dict.values())
    
    # Open all input reads and output files, and run through input reads
    with ExitStack() as stack:
        # Add input reads and output files to context stack
        for input_read in input_reads:
            stack.enter_context(input_read)
        
        for output_file in output_files:
            stack.enter_context(output_file)
            
        # Move to next input read in each file
        i = 0
        while input_reads[0].head and (not debug_mode or i < 100):
            i += 1
            for input_read in input_reads:
                input_read.nextRead()
            
            # Get barcodes and UMIs
            cur_bars = {bar_name:Barcode() for bar_name in bar_pos_dict.keys()}
            cur_UMI = ""
            UMI_qual = ""
            input_read_bad=False
            
            for input_read in input_reads:
                input_read_output = input_read.getBarcodes(bar_pos_dict, debug_mode)
                if input_read_output == False:
                    input_read_bad = True
                    break
                else:
                    input_read_UMI, input_read_UMI_qual, input_read_bars = input_read_output
                
                for name, bar in input_read_bars.items():
                    cur_bars[name] = bar
                cur_UMI += input_read_UMI
                UMI_qual += input_read_UMI_qual
                
            
            # If any of the input reads didn't work out: continue to next line
            if input_read_bad:
                continue
                
            # Trim input reads with reverse complement
            for input_read in input_reads:
                if input_read.trim_partner:
                    partner_seq = input_read_dict[input_read.trim_partner].bar_seq
                    input_read.trimRevComp(partner_seq)
                    if len(input_read.seq) < input_read.min_length:
                        input_read_bad = True
                        break
            
            # If any input reads are too short after trimming, continue to next line
            if input_read_bad:
                continue
            
            # Determine output targets from barcodes using set union
            all_output_targets = set()
            for bar in cur_bars.values():
                all_output_targets |= bar.output_targets  # Set union
            
            all_output_targets = list(all_output_targets)
            
            if debug_mode:
                print("All output targets: %s" % all_output_targets)
            
            for output_target in all_output_targets:
                if output_target not in output_file_dict.keys():
                    raise ValueError("Output target %s not found in output file dictionary" % output_target)
            
            
            if not all_output_targets:
                all_output_targets = default_output_files
            
            # Write to output files
            if debug_mode:
                print("Writing to output files: %s" % all_output_targets)
            for output_target in all_output_targets:
                output_file_dict[output_target].write(cur_bars, input_read_dict, cur_UMI, UMI_qual)

# Distribute the processing of multiple samples across multiple cores
def barcodeReadsMulti(input_folder, output_folder, barcode_folder, cores, debug_mode=False):
    sample_list = []
    with open(os.path.join(barcode_folder,"Samples.csv"), 'rt') as f:
        f.readline()
        for line in f:
            sample_list.append(line.strip().split(',')[0])
            
    p = Pool(processes = int(cores))
    func = partial(barcodeReadsSample, barcode_folder=barcode_folder, input_folder=input_folder, output_folder=output_folder, debug_mode=debug_mode)
    p.map(func, sample_list)
    p.close()
    p.join()


if __name__ == '__main__':
    fire.Fire(barcodeReadsMulti)