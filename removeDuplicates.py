
'''
This script is designed for SAM files that are sorted in the following order:
1) Read properties (which should be in the order of Condition,Cell Name,Read Type,UMI,Read Identifier). It is important that the UMI and read identifier come last.
2) Chromosome
3) Starting position

This should result in read pairs being put together. 

Then the duplicate removal script works as follows:
-Check if any of cell name, read type have changed. If so, we are now on a new batch; clear UMI list.
-If this is the first read from this batch, assume not a duplicate, skip ahead. Before that, write out all reads from the previous batch, along with their duplicate numbers. Add the duplicate numbers to the running counter of duplicate numbers, which will be written to the output file.
-If this is paired end, check that the next read matches (same read name including identifier). If not, start with that next read and continue (in case that next read is the first of a correct pair).
-Check if the starting position (and, if paired end, ending position) are within X bases of any stored reads
-If so, check if the UMI is <= 1 edit distance away. You can optimize this computationally at some point.
-If not a match, add the UMI, start/ending position (as a tuple) to the UMI list. Append the read name (and the mate name, if needed) to the list of reads to be written out.
-If a match, add a duplicate to the appropriate read.
'''
import sys
import os
from multiprocessing import Pool
from functools import partial
from collections import Counter
sys.setrecursionlimit(1000000)


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


def rm_dup_samfile(sam_name, input_folder, output_folder, paired=False, min_dist=30, UMI_correct_dist=1):
    input_file = os.path.join(input_folder, "%s.sam" % sam_name)
    output_file = os.path.join(output_folder, "%s.sam" % sam_name)
    output_csv = os.path.join(output_folder, "%s_duplicate_histogram.csv" % sam_name)
    
    with open(samfile, 'rt') as f_in, open(output_file, 'w') as f_out:
        
        line = f_in.readline()
        while line.startswith('@'):
            f_out.write(line)
            line = f.readline()
        
        cur_cell = ''
        cur_read_type = ''
        dup_counter = Counter()
        while line:

            name, flag, chrom, start = line.split('\t')[:4]
            cell, condition, read_type, UMI, identifier = name.split(',')
            
            if cell != cur_cell or read_type != cur_read_type: # new batch
                cur_UMIs = {}
                cur_cell = cell
                cur_read_type = read_type
                
                

            
            if paired: # THIS NEEDS TO BE FINISHED
                line2 = f.readline()
                name2, _, chrom2, start2 = line2.split('\t')[:4]
                if name != name2:
                    print("WARNING: improper read pairing! Read will be discarded.")
                    line = line2
                    continue
            
            if cell != prev_cell or read_type != prev_read_type:
                f_out.write(duplicate
                UMI_list = []
                prev_cell = cell
                prev_read_type = read_type
                duplicate_num = 0
                
                
            
            cell, UMI, condition, read_type = (((line.split('\t'))[0]).split(','))
            chrom_num = (line.split('\t'))[2]
            start_site = int((line.split('\t'))[3])

            if (abs(start_site-pre_site) <= min_dist) and (chrom_num == pre_chrom):    
                dup = False
                for each_barcode in pre_barcode:
                    if each_barcode == barcode_UMI:
                        dup = True
                        break
                if dup == False:
                    pre_dup_num = cur_dup_num
                    cur_dup_num = 1
                    f2.write(line)
                    pre_barcode.add(barcode_UMI)
                    f3.write('%d' %(pre_dup_num))
                    f3.write('\n')
                else:
                    cur_dup_num += 1               
            else:
                pre_dup_num = cur_dup_num
                cur_dup_num = 1
                f2.write(line)
                pre_chrom = chrom_num
                pre_site = start_site
                pre_barcode = set()
                pre_barcode.add(barcode_UMI)
                if (pre_dup_num != 0):
                    f3.write("%d" % (pre_dup_num))
                    f3.write('\n')

    
    '''
    #plot the histogram for the read duplication number
    dups = (pd.read_csv(output_file+'.csv', header=None))[0]
    fig = plt.figure()
    plt.hist(dups, bins=100)
    plt.xlabel("Duplication number")
    plt.ylabel("Read number")
    fig.savefig(output_file + '.png')
    '''

if __name__ == "__main__":
    input_folder=sys.argv[1]
    sampleID=sys.argv[2]
    output_folder=sys.argv[3]
    core=sys.argv[4]
    
    
    sample_list = []
    with open(sampleID, 'rt') as f:
        for line in f:
            sample_list.append(line.strip())
    
    p = Pool(processes = int(core))
    func = partial(rm_dup_samfile, input_folder=input_folder, output_folder=output_folder)
    result = p.map(func, sample_list)
    p.close()
    p.join()