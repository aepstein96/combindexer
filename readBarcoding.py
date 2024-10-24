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


def getLevRecursive(bar_list, distance=1, bases={'A','C','G','T','N'}, recursed=False):
    
    if not isinstance(bar_list, list):
        bar_list = [bar_list]
    
    if distance < 1:
        return bar_list
    
    if not recursed: # shouldn't be any duplicates from calling it previously
        if len(set(bar_list)) != len(bar_list):
            raise ValueError("Input list contains duplicates.")
    
    bar_length = len(bar_list[0])
    out_list = bar_list.copy()
    for bar in bar_list:
        if len(bar) != bar_length:
            return ValueError("Input barcodes are different lengths.")
        
        for pos in range(bar_length):
            for base in bases - {bar[pos]}: #don't replace base with itself
                out_list.append(''.join((bar[:pos], base, bar[pos+1:])))
                
    return getLevRecursive(out_list, distance-1, recursed=True)

# For a list of correct barcodes, creates a dictionary mapping uncorrected barcodes to correct barcodes
def getLevDict(correct_bars, default=None, **kwargs):
    out_dict = defaultdict(lambda: default) # returns default if a non-matching barcode is inputted
    for correct_bar in correct_bars:
        for bar in getLevRecursive(correct_bar, **kwargs):
            out_dict[bar] = correct_bar
        
    return out_dict

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
                bar = bar_pos.getBar(self.seq[start:start+bar_pos.length])
                
                if not bar:
                    if output_full:
                        print("No match for %s. Start: %d. File: %s. Read sequence: %s" % (bar_name, start, self.in_path, self.seq))
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
    def __init__(self, bar_list, length_to_use, plate_list, corr_dist, orientation):
        
        # Get length of shortest barcode, this is what will be used for disambiguation
        # Also see if there is a default for "No match"
        self.no_match = None
        min_length = length_to_use
        true_bars = [] # doesn't include no match
        for bar in bar_list:
            if bar.seq == 'No match':
                self.no_match = bar
                continue
                
            if bar.length < min_length:
                min_length = bar.length
                
            true_bars.append(bar)
        
        self.length = min_length
        self.plate_list = plate_list
        self.orientation = orientation
                
        # Make barcode dictionary
        self.bar_dict = defaultdict(lambda: self.no_match)
        for bar in true_bars:
            if plate_list:
                if bar.plate not in plate_list:
                    continue
                    
            if orientation == 'reverse':
                self.bar_dict[revComp(bar.seq)[:min_length]] = bar
            else:
                self.bar_dict[bar.seq[:min_length]] = bar
        
        # Make correction dictionary
        
        self.corr_dict = getLevDict(self.bar_dict.keys(), distance=corr_dist)
    
    # Returns barcode information
    def getBar(self, bar_seq):
        return self.bar_dict[self.corr_dict[bar_seq]]

    
def makeBarDict(barcode_folder, output_full=False):
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
                
                if output_full:
                    print("Corr dist: %d" % corr_dist)
                orientation = f.readline().split(',')[1].strip()
                f.readline()
                f.readline()
                
                bar_list = []
                for line in f:
                    seq, plate, well, condition, read_type, gene, next_command_string = line.strip().split(',')[:7]
                    
                    if not seq: # blank line at end of file
                        break
                    
                    if next_command_string:
                        next_commands = next_command_string.split(';')
                    else:
                        next_commands = None
                    bar_list.append(Barcode(seq, plate, well, condition, read_type, gene, next_commands))
                    
            bar_dict[name] = BarcodePos(bar_list, length_to_use, plate_list, corr_dist, orientation) 
    return bar_dict

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
                cur_bars_ordered = [cur_bars[bar_name] for bar_name in ordered_bar_list] # probably not necessary now that dictionaries are ordered in Python 3.7+
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
    
    