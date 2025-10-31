from collections import defaultdict, deque
from contextlib import ExitStack
from multiprocessing import Pool
from functools import partial
from collections import Counter
import gzip
import glob
import sys
import os
import fire

# For a given barcode, returns list of barcodes a given Levenshtein distance away


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

# Returns reverse complement of sequence
# Not as flexible as Bio.Seq, but probably faster and easier to use
def revComp(seq, comp_dict={'A':'T','G':'C','C':'G','T':'A','N':'N'}):
    return ''.join([comp_dict[base] for base in seq[::-1]])


class Read:
    def __init__(self, name, sample, commands, in_pattern, out_patterns=None, trim_partner=None, trim_partner_len=17, min_length=20, debug_mode=False):
        self.name = name
        self.sample = sample
        self.commands = commands
        self.in_pattern = in_pattern
        self.out_patterns = out_patterns or []  # List of output file patterns
        self.trim_partner = trim_partner # should be in order of position
        self.trim_partner_len = trim_partner_len
        self.sample = None
        self.min_length = min_length
        self.out_files = []  # List of open file handles
        
        if debug_mode:
            print("Input pattern: %s. Sample: %s" % (self.in_pattern, sample))
            
        self.in_paths = glob.glob(self.in_pattern.replace("{SAMPLE}", sample))
        if self.in_paths:
            self.in_path = self.in_paths[0]
        else:
            raise ValueError("%s not found" % self.in_pattern.replace("{SAMPLE}", sample))
            
        # Assign multiple output file paths
        self.out_paths = []
        for pattern in self.out_patterns:
            if pattern:  # Skip empty patterns
                self.out_paths.append(pattern.replace("{SAMPLE}", sample))
            
    
    def __enter__(self):
        if self.in_path.endswith('.gz'):
            self.input = gzip.open(self.in_path, 'rt') #rt makes it output text and not binary
        else:
            self.input = open(self.in_path, 'rt')
        self.head = 'placeholder' # required to start the while loop
        
        # Open all output files
        self.out_files = []
        for out_path in self.out_paths:
            if out_path.endswith('.gz'):
                self.out_files.append(gzip.open(out_path, 'wt'))
            else:
                self.out_files.append(open(out_path, 'wt'))
    
    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.input.close()
        
        # Close all output files
        for out_file in self.out_files:
            out_file.close()
    
    def nextRead(self):
        # Read in data from file
        self.head = self.input.readline().rstrip() # remove newline marker, it will be put back later
        self.seq = self.input.readline().rstrip()
        self.sep = self.input.readline().rstrip()
        self.qual = self.input.readline().rstrip()
        
    def getBarcodes(self, barcode_dict, debug_mode=False):
            
        start = 0
        cur_bars = {}
        UMI = ''
        
        if not self.commands or not self.commands[0]: # this is what happens if that read has no barcodes in it
            self.bar_seq = ''
            return UMI, cur_bars
        
        command_queue = deque(self.commands)
        #print("getBarcodes triggered! Sample: %s. Read name: %s" % (self.sample, self.name))
        trim_pos = 0
        while command_queue:
            command_name = command_queue.popleft()
            
            if command_name.startswith('Goto:'):
                start = int(command_name[5:])
                
            elif command_name.startswith('Gap:'):
                start += int(command_name[4:])
                #print("Start after gap: %d" % start)
                
            elif command_name.startswith('UMI:'):
                UMI_len = int(command_name[4:])
                cur_UMI = self.seq[start:start+UMI_len]
                if 'N' in cur_UMI:
                    return False # this is kind of inelegant, maybe there is a better way to do this
                
                UMI += cur_UMI
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
                        print("No match for %s. Start: %d. File: %s. Read sequence: %s. Test sequence: %s" % (bar_name, start, self.in_path, self.seq, potential_bar_seq))
                    return False
                
                #print("Bar name: %s. Bar sequence: %s. Read bar sequence: %s. Length: %d" % (bar_name, bar.seq, self.seq[start:start+bar_pos.length], bar.length))
                
                start += bar.length
                cur_bars[bar_name] = bar
                if bar.next_commands:
                    #print("Found next commands! Remaining commands: %s" % str(bar.next_commands))
                    command_queue = deque(bar.next_commands)
            else:
                raise ValueError("Command not found: %s. Valid commands: Goto:<pos>, Gap:<len>, UMI:<len>, Trim, Bar:<name>" % command_name)
                
        self.bar_seq = self.seq[:start] # used for trimming partner
        if debug_mode:
            print("All barcodes found! Read: %s, Barcodes: %s\n" % (self.name, str([bar.seq for bar in cur_bars.values()])))
        self.seq = self.seq[trim_pos:]
        self.qual = self.qual[trim_pos:]
            
        return UMI, cur_bars
        
    def trimRevComp(self, partner_seq):
        if len(partner_seq) < self.trim_partner_len: # Don't trim based on too short a sequence
            return False
        
        search_seq = revComp(partner_seq)[:self.trim_partner_len]
        search_seq_pos = self.seq.find(search_seq)
        if search_seq_pos != -1:
            self.seq = self.seq[:search_seq_pos]
            self.qual = self.qual[:search_seq_pos]
        
        return True

    
    def write(self, file_numbers, header_data, seq=None, sep=None, qual=None):
        """Write FASTQ records to specified output file numbers"""
        if not self.out_files:
            return
            
        # Use defaults if not provided
        if seq is None:
            seq = self.seq
        if sep is None:
            sep = self.sep
        if qual is None:
            qual = self.qual
        
        # Extract identifier from original header
        identifier = self.head.split(' ')[0][1:]
        fastq_header = '@' + header_data + ',' + identifier
        
        # Write to all specified file numbers
        for file_number in file_numbers:
            # Convert to 0-indexed
            file_idx = file_number - 1
            
            if file_idx >= len(self.out_files):
                raise ValueError(f"Output file {file_number} not defined for read {self.name}")
            
            out_file = self.out_files[file_idx]
            
            # Write FASTQ record
            out_file.write(fastq_header + '\n')
            out_file.write(seq + '\n')
            out_file.write(sep + '\n')
            out_file.write(qual + '\n')
    
class Barcode:
    def __init__(self, seq='', plate='', well='', condition='', read_type='', next_commands=[], output_files=''):
        self.seq = seq # full sequence, including part that is not used
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
        
class BarcodePos:
    def __init__(self, bar_list, length_to_use, plate_list, corr_dist, orientation, correct_snps, correct_indels, debug_mode=False):
        
        # Get length of shortest barcode, this is what will be used for disambiguation
        # Also see if there is a default for "No match"
        self.no_match = None
        self.ambiguous = None
        min_length = length_to_use
        true_bars = [] # doesn't include no match
        for bar in bar_list:
            if bar.seq == 'No match':
                self.no_match = bar
                continue
            if bar.seq == 'Ambiguous':
                self.ambiguous = bar
                continue
                
            if bar.length < min_length:
                min_length = bar.length
                
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

                
            variants_map = getCombinedVariants(bar_seq, corr_dist, post_context, correct_snps, correct_indels)
            for var_seq, (_, dist) in variants_map.items():
                potential_corrections[var_seq].append((bar, dist))

        self.corr_dict = defaultdict(lambda: self.no_match)
        for var_seq, candidates in potential_corrections.items():
            correct_bar = findCorrectBarcode(var_seq, candidates)
            if correct_bar:
                if debug_mode:
                    print("Var seq: %s. Correct bar: %s" % (var_seq, correct_bar.seq))
                self.corr_dict[var_seq] = correct_bar
            else:
                if debug_mode:
                    print("Var seq: %s. No correct bar found" % var_seq)
                    
                if self.ambiguous:
                    self.corr_dict[var_seq] = self.ambiguous
                        
    
    # Returns barcode information
    def getBar(self, bar_seq):
        return self.corr_dict[bar_seq]


def findCorrectBarcode(var_seq, candidates):
    if len(candidates) == 1:
        return candidates[0][0]

    min_dist = min(c[1] for c in candidates)
    best_matches = [c[0] for c in candidates if c[1] == min_dist]

    if len(best_matches) == 1:
        return best_matches[0]
    
    return False


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

def makeReadDict(sample, barcode_folder, input_folder, output_folder, debug_mode=False):
    # Get read information
    read_dict = {}
    reads_to_align = []
    with open(os.path.join(barcode_folder, "Reads.csv"), 'rt') as f:
        f.readline() # skips header

        # Iterate through reads
        for line in f:
            parts = line.strip().split(',')
            name = parts[0]
            in_file_pattern = parts[1]
            command_str = parts[2]
            trim_partner_str = parts[3] if len(parts) > 3 else ""
            min_length = parts[4] if len(parts) > 4 else ""
            
            # All remaining parts are output file patterns (moved to the end)
            out_file_patterns = parts[5:] if len(parts) > 5 else []
            
            in_pattern = os.path.join(input_folder, in_file_pattern)
            
            # Process output patterns and add to reads_to_align
            out_patterns = []
            for pattern in out_file_patterns:
                if pattern:  # Skip empty patterns
                    out_pattern = os.path.join(output_folder, pattern)
                    out_patterns.append(out_pattern)
                    reads_to_align.append(pattern)
            
            commands = command_str.strip().split(';') if command_str else []
            if trim_partner_str:
                trim_partner, trim_partner_len = trim_partner_str.strip().split(':')
                trim_partner_len = int(trim_partner_len)
            else:
                trim_partner = None
                trim_partner_len = None
            
            if min_length:
                min_length = int(min_length)
                
            read = Read(name, sample, commands, in_pattern, out_patterns, trim_partner, trim_partner_len, min_length, debug_mode)
            read_dict[name] = read
    
    # Write out read pattern file
    with open(os.path.join(barcode_folder, "reads_to_align.txt"), 'wt') as f:
        for read in reads_to_align:
            f.write(read + "\n")

    return read_dict
            

def barcodeReadsSample(sample, barcode_folder, input_folder, output_folder, default_write_setting="1;csv", debug_mode=False):
    
    # Make dictionaries (temporary until I can fix the multiprocessing issue)
    read_dict = makeReadDict(sample, barcode_folder, input_folder, output_folder, debug_mode)
    for key, read in read_dict.items():
        if debug_mode:
            print("Read name: %s. Commands: %s" % (key, read.commands))
    

    bar_pos_dict = makeBarDict(barcode_folder, debug_mode)                                                          
            
    # Read in barcodes in the correct order
    # Use provided order if it exists, otherwise use alphabetical order

    ordered_bar_list = []
    if os.path.exists(os.path.join(barcode_folder, 'sheet_names.txt')):
        with open(os.path.join(barcode_folder, 'sheet_names.txt'), 'rt') as f:
            sheets = f.readline().strip().split(',')
            ordered_bar_list = [sheet.split('-')[1] for sheet in sheets if sheet.startswith('Barcode-')]
            if debug_mode:
                print("Ordered bar list:")
                print(ordered_bar_list)
    else:
        files = [f for f in os.listdir(barcode_folder) if f.startswith('Barcode-') and f.endswith('.csv')]
        ordered_bar_list = sorted([f.split('-')[1].split('.')[0] for f in files])
        if debug_mode:
            print("Ordered bar list:")
            print(ordered_bar_list)
 
    # Run through files
    reads = list(read_dict.values())
    
    # Set up CSV file
    csv_file_path = os.path.join(output_folder, f"{sample}_output.csv")
    csv_file = open(csv_file_path, 'wt')
    
    with ExitStack() as stack:
        # Add output files to context stack
        stack.enter_context(csv_file)
        for read in reads:
            stack.enter_context(read)
        
        # Write CSV header
        read_names = [read.name for read in reads]
        csv_header = 'wells_conditions,UMI,read_types,' + ','.join(read_names)
        csv_file.write(csv_header + '\n')
            
        # Move to next read in each file
        i = 0
        while reads[0].head and (not debug_mode or i < 100):
            i += 1
            for read in reads:
                read.nextRead()
            
            # Get barcodes and UMIs
            cur_bars = {bar_name:Barcode() for bar_name in ordered_bar_list}
            UMI = ""
            read_bad=False
            
            for read in reads:
                read_output = read.getBarcodes(bar_pos_dict, debug_mode)
                if read_output == False:
                    read_bad = True
                    break
                else:
                    read_UMI, read_bars = read_output
                
                for name, bar in read_bars.items():
                    cur_bars[name] = bar
                    #print("Barcode found! Name: %s. Plate: %s Well: %s. Condition: %s. Read type: %s" % (name, bar.plate, bar.well, bar.condition, bar.read_type))
                UMI += read_UMI
            
            # If any of the reads didn't work out: continue to next line
            if read_bad:
                continue
                
            # Trim reads with reverse complement
            for read in reads:
                if read.trim_partner:
                    partner_seq = read_dict[read.trim_partner].bar_seq
                    read.trimRevComp(partner_seq)
                    if len(read.seq) < read.min_length:
                        read_bad = True
                        break
                        
            if read_bad:
                continue
            
            # Determine output targets from barcodes using set union
            all_output_targets = set()
            for bar in cur_bars.values():
                all_output_targets |= bar.output_targets  # Set union
            
            # If no barcode specifies outputs, use default
            if not all_output_targets:
                default_targets = [target.strip() for target in default_write_setting.split(';') if target.strip()]
                all_output_targets = set(default_targets)
            
            # Separate FASTQ file numbers from CSV using set operations
            write_to_csv = 'csv' in {target.lower() for target in all_output_targets}
            
            # Get numeric targets only
            numeric_targets = all_output_targets - {'csv', 'CSV'}  # Remove csv variants
            fastq_file_numbers = []
            
            for target in numeric_targets:
                try:
                    fastq_file_numbers.append(int(target))
                except ValueError:
                    raise ValueError(f"Invalid output target '{target}'. Must be a number or 'csv'.")
            
            # Construct header from barcodes and UMIs
            cur_bars_ordered = list(cur_bars.values()) # Dictionaries are ordered in Python 3.7+
            plate_wells = ';'.join(["%s-%s" % (bar.plate, bar.well) for bar in cur_bars_ordered])
            conditions = ';'.join([bar.condition for bar in cur_bars_ordered])
            wells_conditions = '&'.join((plate_wells, conditions))
            read_types = ';'.join([bar.read_type for bar in cur_bars_ordered])
            
            # Construct header data
            header_data = ','.join((wells_conditions,UMI,read_types))
            
            # Write to FASTQ files
            if fastq_file_numbers:
                for read in reads:
                    if read.out_files:  # Only write if there are output files
                        read.write(fastq_file_numbers, header_data=header_data)
            
            # Write to CSV file
            if write_to_csv:
                # Add read sequences to the CSV row
                read_sequences = [read.seq for read in reads]
                csv_row = header_data + ',' + ','.join(read_sequences)
                csv_file.write(csv_row + '\n')
        
    
def barcodeReadsMulti(input_folder, output_folder, barcode_folder, cores, default_write_setting="1;csv", debug_mode=False):
    #read_dict = makeReadDict(barcode_folder, input_folder, output_folder)
    #barcode_dict = makeBarDict(barcode_folder)
    
    sample_list = []
    with open(os.path.join(barcode_folder,"Samples.csv"), 'rt') as f:
        f.readline()
        for line in f:
            sample_list.append(line.strip().split(',')[0])
            
    p = Pool(processes = int(cores))
    func = partial(barcodeReadsSample, barcode_folder=barcode_folder, input_folder=input_folder, output_folder=output_folder, default_write_setting=default_write_setting, debug_mode=debug_mode)
    result = p.map(func, sample_list)
    p.close()
    p.join()
                

if __name__ == '__main__':
    fire.Fire(barcodeReadsMulti)
    
    