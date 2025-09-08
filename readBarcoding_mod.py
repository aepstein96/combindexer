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
import itertools
import csv

def getCombinedVariants(
    mid_seq, 
    distance,
    correct_snps,
    correct_indels,
    post_context="", 
    bases={'A', 'C', 'G', 'T', 'N'}
):
    """
    Generates all unique, fixed-length variants of a sequence within a given
    combined correction distance.

    This version has been corrected to remove the flawed pre_context logic.
    All deletions are now correctly compensated by the post_context.
    """
    # The state now only needs to track the variant and its post_context
    all_variants_map = {mid_seq: post_context}
    last_generation = [(mid_seq, post_context)]

    for _ in range(distance):
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
                                all_variants_map[new_variant] = post
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
                            all_variants_map[new_variant] = new_post
                            current_generation.append((new_variant, new_post))
                    else: # Fallback: No post_context, pad with all bases
                        for b in bases:
                            new_variant = deleted_in_slice + b
                            if new_variant not in all_variants_map:
                                all_variants_map[new_variant] = ""
                                current_generation.append((new_variant, ""))
                
                # Insertions anywhere (compensated by truncation)
                for i in range(len(var) + 1):
                    for b in bases:
                        inserted = var[:i] + b + var[i:]
                        new_variant = inserted[:len(mid_seq)]
                        if new_variant not in all_variants_map:
                            # New post_context is the truncated char + old post
                            new_post = inserted[-1] + post
                            all_variants_map[new_variant] = new_post
                            current_generation.append((new_variant, new_post))
        
        last_generation = current_generation
        
    return set(all_variants_map.keys())


# For a list of correct barcodes, creates a dictionary mapping uncorrected barcodes to correct barcodes
def getCorrDict(bar_dict, default=None, correction_distance=0, correct_snps=False, correct_indels=False, 
                check_overlaps=False, ambiguous_barcode=None, 
                orientation='forward', bases={'A', 'C', 'G', 'T', 'N'}):
    
    corr_dict = defaultdict(lambda: default)
    ambiguous_variants = set()

    for canonical_slice, bar_obj in bar_dict.items():
        
        full_seq = bar_obj.seq
        if orientation == 'reverse':
            full_seq = revComp(full_seq)
        
        # The length is implicitly defined by the canonical slice
        length = len(canonical_slice)
        post_context = full_seq[length:]

        variants_to_check = getCombinedVariants(
            canonical_slice,
            correction_distance,
            correct_snps,
            correct_indels,
            post_context,
            bases
        )

        for var in variants_to_check:
            if var in ambiguous_variants:
                continue

            # If check_overlaps is on and we find a collision
            if check_overlaps and var in corr_dict and corr_dict[var] != canonical_slice:
                ambiguous_variants.add(var)
                # If an 'Ambiguous' barcode is defined, map to it, otherwise remove the key
                if ambiguous_barcode:
                    corr_dict[var] = 'Ambiguous' # Use a placeholder
                else:
                    if var in corr_dict:
                        del corr_dict[var]
            elif var not in corr_dict: # Don't overwrite an existing mapping
                corr_dict[var] = canonical_slice

    # Now, build the final dictionary, resolving placeholders to actual objects
    final_dict = defaultdict(lambda: default)
    for variant, result_key in corr_dict.items():
        if result_key == 'Ambiguous':
            final_dict[variant] = ambiguous_barcode
        else:
            final_dict[variant] = bar_dict[result_key]
            
    return final_dict

# Returns reverse complement of sequence
# Not as flexible as Bio.Seq, but probably faster and easier to use
def revComp(seq, comp_dict={'A':'T','G':'C','C':'G','T':'A','N':'N'}):
    return ''.join([comp_dict[base] for base in seq[::-1]])


class Read:
    def __init__(self, name, commands, in_pattern, out_pattern=None, trim_partner=None, trim_partner_len=17, min_length=20):
        self.name = name
        self.commands = commands
        self.in_pattern = in_pattern
        self.out_pattern = out_pattern
        self.trim_partner = trim_partner # should be in order of position
        self.trim_partner_len = trim_partner_len
        self.sample = None
        self.min_length = min_length
        
        # Paths are assigned through assignSample
        self.in_path = None
        self.out_path = None
       
    
    def assignSample(self, sample, output_full): # assigns the sample name, updating input and output paths
        if output_full:
            print("Input pattern: %s. Sample: %s" % (self.in_pattern, sample))
            
        self.sample = sample
        self.in_paths = glob.glob(self.in_pattern.replace("{SAMPLE}", sample))
        if self.in_paths:
            self.in_path = self.in_paths[0]
        else:
            raise ValueError("%s not found" % self.in_pattern.replace("{SAMPLE}", sample))
            
        if self.out_pattern:
            self.out_path = self.out_pattern.replace("{SAMPLE}", sample)
            
    
    def __enter__(self):
        if self.in_path.endswith('.gz'):
            self.input = gzip.open(self.in_path, 'rt') #rt makes it output text and not binary
        else:
            self.input = open(self.in_path, 'rt')
        self.head = 'placeholder' # required to start the while loop
        
        if self.out_path:
            if self.out_path.endswith('.gz'):
                self.out = gzip.open(self.out_path, 'wt')
            else:
                self.out = open(self.out_path, 'wt')
    
    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.input.close()
        
        if self.out_path:
            self.out.close()
    
    def nextRead(self):
        # Read in data from file
        self.head = self.input.readline().rstrip() # remove newline marker, it will be put back later
        self.seq = self.input.readline().rstrip()
        self.sep = self.input.readline().rstrip()
        self.qual = self.input.readline().rstrip()
        
    def getBarcodes(self, barcode_dict, output_full=False):
            
        start = 0
        cur_bars = {}
        UMI = ''
        
        if not self.commands[0]: # this is what happens if that read has no barcodes in it
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
                
                # Since all variants are now the same length, we can do a direct lookup.
                if start + bar_pos.length > len(self.seq):
                    if output_full:
                        print(f"Not enough sequence left for barcode {bar_name}")
                    return False

                read_bar_seq = self.seq[start:start+bar_pos.length]
                bar = bar_pos.getBar(read_bar_seq)
                
                if not bar:
                    if output_full:
                        print("No match for %s. Start: %d. File: %s. Read sequence: %s" % (bar_name, start, self.in_path, self.seq))
                    return False
                
                # print("Bar name: %s. Bar sequence: %s. Read bar sequence: %s. Length: %d" % (bar_name, bar.seq, read_bar_seq, bar.length))
                
                # Advance start position by the length of the barcode segment that was read.
                start += bar_pos.length

                cur_bars[bar_name] = bar
                if bar.next_commands:
                    #print("Found next commands! Remaining commands: %s" % str(bar.next_commands))
                    command_queue = deque(bar.next_commands)
            else:
                raise ValueError("Command not found: %s. Valid commands: Goto:<pos>, Gap:<len>, UMI:<len>, Trim, Bar:<name>" % command_name)
                
        self.bar_seq = self.seq[:start] # used for trimming partner
        if output_full:
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

    
    def write(self, head=None, seq=None, sep=None, qual=None):
        if self.out_path:
            self.out.write((head or self.head) + "\n")
            self.out.write((seq or self.seq) + "\n")
            self.out.write((sep or self.sep) + "\n")
            self.out.write((qual or self.qual) + "\n")
    
class Barcode:
    def __init__(self, seq='', plate='', well='', condition='', read_type='', gene='', next_commands=[]):
        self.seq = seq # full sequence, including part that is not used
        self.plate = plate
        self.well = well
        self.condition = condition
        self.read_type = read_type
        self.gene = gene
        self.next_commands = next_commands
        self.length = len(seq)
        
class BarcodePos:
    def __init__(self, bar_list, length, plate_list, correction_distance, 
                 correct_snps, correct_indels, check_overlaps, orientation):
        
        # Store parameters
        self.length = length
        self.plate_list = plate_list
        self.correction_distance = correction_distance
        self.correct_snps = correct_snps
        self.correct_indels = correct_indels
        self.check_overlaps = check_overlaps
        self.orientation = orientation
        self.no_match = None
        self.ambiguous = None
        
        # Find the "No match" and "Ambiguous" barcodes
        true_bars = []
        for bar in bar_list:
            if bar.seq == 'No match':
                self.no_match = bar
            elif bar.seq == 'Ambiguous':
                self.ambiguous = bar
            else:
                true_bars.append(bar)
                if len(bar.seq) < self.length:
                    self.length = len(bar.seq)
                
        if not true_bars:
             self.bar_dict = {}
             self.corr_dict = defaultdict(lambda: self.no_match)
             return

        # If length is 'all', find the minimum length from the barcode list
        if self.length == 'all':
            self.length = min(len(bar.seq) for bar in true_bars) if true_bars else 0
        else:
            self.length = int(self.length)
        
        # Make barcode dictionary mapping canonical slice to barcode object
        self.bar_dict = {}
        for bar in true_bars:
            if self.plate_list and self.plate_list[0] != 'all' and bar.plate not in self.plate_list:
                continue
            
            full_seq = bar.seq
            if self.orientation == 'reverse':
                full_seq = revComp(full_seq)

            if len(full_seq) < self.length:
                continue
            
            key_seq = full_seq[:self.length]
            self.bar_dict[key_seq] = bar
        
        # Make correction dictionary
        if self.correction_distance > 0:
            self.corr_dict = getCorrDict(
                self.bar_dict,
                default=self.no_match,
                correction_distance=self.correction_distance,
                correct_snps=self.correct_snps,
                correct_indels=self.correct_indels,
                check_overlaps=self.check_overlaps,
                ambiguous_barcode=self.ambiguous,
                orientation=self.orientation
            )
        else:
            self.corr_dict = defaultdict(lambda: self.no_match, self.bar_dict)


    # Returns barcode information
    def getBar(self, bar_seq):
        return self.corr_dict[bar_seq]

    
def makeBarDict(barcode_folder, output_full=False):
    
    barcode_pos_dict = {} # a dictionary of BarcodePos objects, keyed by the barcode name
    
    barcode_files = glob.glob(os.path.join(barcode_folder, "Barcode-*.csv"))
    
    for bar_file in barcode_files:
        name = os.path.basename(bar_file).replace("Barcode-", "").replace(".csv", "")
        
        with open(bar_file, 'r') as file:
            lines = file.readlines()
            
            # Parse settings from the first 7 lines
            settings_lines = [line.strip().split(',') for line in lines[:7]]
            length_to_use = settings_lines[0][1]
            if length_to_use == 'all':
                length_to_use = 'all'
            else:
                length_to_use = int(length_to_use)
            
            plates_to_use = settings_lines[1][1].split(';')
            correction_distance = int(settings_lines[2][1])
            orientation = settings_lines[3][1]
            correction_type = settings_lines[4][1].lower()
            correct_snps = 'snp' in correction_type
            correct_indels = 'indel' in correction_type
            check_overlaps = settings_lines[6][1].lower() == 'true'

            # The 8th line is the header for the barcode data
            header = [h.strip() for h in lines[7].split(',')]
            
            bar_list = []
            
            # Use DictReader on the rest of the lines
            reader = csv.DictReader(lines[8:], fieldnames=header)
            for row in reader:
                if not row.get('Barcode'):
                    continue
                
                bar_list.append(
                    Barcode(
                        seq=row.get('Barcode', ''),
                        plate=row.get('Plate', ''),
                        well=row.get('Well', ''),
                        condition=row.get('Condition', ''),
                        read_type=row.get('Read Type', ''),
                        gene=row.get('Gene', ''),
                        next_commands=row.get('Next Commands', '').split(';') if row.get('Next Commands') else []
                    )
                )

            # Create the BarcodePos object with all the data
            barcode_pos_dict[name] = BarcodePos(
                bar_list=bar_list,
                length=length_to_use,
                plate_list=plates_to_use,
                correction_distance=correction_distance,
                correct_snps=correct_snps,
                correct_indels=correct_indels,
                check_overlaps=check_overlaps,
                orientation=orientation
            )

    return barcode_pos_dict

def makeReadDict(barcode_folder, input_folder, output_folder):
    # Get read information
    read_dict = {}
    reads_to_align = []
    with open(os.path.join(barcode_folder, "Reads.csv"), 'rt') as f:
        f.readline() # skips header

        # Iterate through reads
        for line in f:
            name, in_file_pattern, out_file_pattern, command_str, trim_partner_str, min_length = line.strip().split(',')
            in_pattern = os.path.join(input_folder, in_file_pattern)
            if out_file_pattern:
                out_pattern = os.path.join(output_folder, out_file_pattern)
                reads_to_align.append(out_file_pattern)
            else:
                out_pattern = None
            
            commands = command_str.strip().split(';')
            if trim_partner_str:
                trim_partner, trim_partner_len = trim_partner_str.strip().split(':')
                trim_partner_len = int(trim_partner_len)
            else:
                trim_partner = None
                trim_partner_len = None
            
            if min_length:
                min_length = int(min_length)
                
            read = Read(name, commands, in_pattern, out_pattern, trim_partner, trim_partner_len, min_length)
            read_dict[name] = read
    
    # Write out read pattern file
    with open(os.path.join(barcode_folder, "reads_to_align.txt"), 'wt') as f:
        for read in reads_to_align:
            f.write(read + "\n")

    return read_dict
            
def barcodeReadsSample(sample, barcode_folder, input_folder, output_folder, gene_output_folder, output_full=False):
    
    # Make dictionaries (temporary until I can fix the multiprocessing issue)
    read_dict = makeReadDict(barcode_folder, input_folder, output_folder)
    for key, read in read_dict.items():
        if output_full:
            print("Read name: %s. Commands: %s" % (key, read.commands))
    
        
        
    barcode_dict = makeBarDict(barcode_folder, output_full)
    for key, bar in barcode_dict.items():
        if output_full:
            print("Barcode name: %s. Length: %s. Plates: %s. Orientation: %s" % (key, bar.length, bar.plate_list, bar.orientation))
            
        if key == 'SSSprimer' or key == 'Gene':
            if output_full:
                print("Corr dict: %s" % bar.corr_dict)
    # Read in barcodes in the correct order
    ordered_bar_list = []
    with open(os.path.join(barcode_folder, 'sheet_names.txt'), 'rt') as f:
        sheets = f.readline().strip().split(',')
        ordered_bar_list = [sheet.split('-')[1] for sheet in sheets if sheet.startswith('Barcode-')]
        if output_full:
            print("Ordered bar list:")
            print(ordered_bar_list)
 
    # Run through files
    reads = list(read_dict.values())
    gene_assigned_file = os.path.join(gene_output_folder, "%s.csv" % sample)
    #genes_assigned = Counter()
    with open(gene_assigned_file, 'wt') as f_gene:
        with ExitStack() as stack:
            for read in reads:
                read.assignSample(sample, output_full)
                stack.enter_context(read)
            
            # Move to next read in each file
            while reads[0].head:
            #for i in range(1000):
                for read in reads:
                    read.nextRead()
                
                # Get barcodes and UMIs
                cur_bars = {bar_name:Barcode() for bar_name in ordered_bar_list}
                UMI = ""
                read_bad=False
                
                for read in reads:
                    read_output = read.getBarcodes(barcode_dict, output_full)
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
                
                # Construct header from barcodes and UMIs
                cur_bars_ordered = list(cur_bars.values()) # Dictionaries are ordered in Python 3.7+
                plate_wells = ';'.join(["%s-%s" % (bar.plate, bar.well) for bar in cur_bars_ordered])
                conditions = ';'.join([bar.condition for bar in cur_bars_ordered])
                wells_conditions = '&'.join((plate_wells, conditions))
                read_types = ';'.join([bar.read_type for bar in cur_bars_ordered])
                genes = ';'.join([bar.gene for bar in cur_bars_ordered if bar.gene])
                
                
                # Write out aligned reads to aligned read file
                if genes:
                    #print("Genes found! Genes: %s" % genes)
                    f_gene.write(','.join((wells_conditions,UMI,read_types,genes)) + '\n')
                    #genes_assigned[(wells_conditions,UMI,read_types,genes)] += 1
                    #print(genes_assigned)
                    #print("Counter sum: %s" % sum(genes_assigned.values()))
                    
                
                # Write out reads to files
                else:
                    #print("No genes found! Genes: %s" % genes)
                    new_header = '@' + ','.join((wells_conditions,UMI,read_types))
                    for read in reads:
                        identifier = read.head.split(' ')[0][1:]
                        read.write(head=','.join((new_header, identifier)))
                
    # Write out gene assigned file
    #with open(gene_assigned_file, 'wt') as f:
    #    print("Genes assigned length: %s" % len(genes_assigned))
    #    for unique_read, num_dups in genes_assigned.items():
    #        #print(unique_read)
    #        read_str = ','.join(unique_read)
    #        f.write(','.join((read_str, str(num_dups))) + '\n')
        
    
def barcodeReadsMulti(input_folder, output_folder, barcode_folder, gene_output_folder, cores, output_full=False):
    #read_dict = makeReadDict(barcode_folder, input_folder, output_folder)
    #barcode_dict = makeBarDict(barcode_folder)
    
    sample_list = []
    with open(os.path.join(barcode_folder,"Samples.csv"), 'rt') as f:
        f.readline()
        for line in f:
            sample_list.append(line.strip().split(',')[0])
            
    p = Pool(processes = int(cores))
    func = partial(barcodeReadsSample, barcode_folder=barcode_folder, input_folder=input_folder, output_folder=output_folder, gene_output_folder=gene_output_folder, output_full=output_full)
    result = p.map(func, sample_list)
    p.close()
    p.join()
                

if __name__ == '__main__':
    #input_folder = sys.argv[1]
    #output_folder = sys.argv[2]
    #barcode_folder = sys.argv[3]
    #gene_output_folder = sys.argv[4]
    #cores = sys.argv[5]
    
    fire.Fire(barcodeReadsMulti)
    
    