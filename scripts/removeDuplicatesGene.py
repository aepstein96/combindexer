import sys
import os
from collections import Counter, defaultdict
from multiprocessing import Pool
from functools import partial
#from umi_tools import UMIClusterer

# NOTE: Input files should be sorted in the following order: 1) cell, 2) gene, 3) read type, 4) UMI
def removeDuplicates(sample, in_folder, out_folder):
    umi_counts = os.path.join(out_folder, "%s_umi_counts.csv" % sample)
    with open(os.path.join(in_folder, "%s.csv" % sample), 'rt') as f_in, open(umi_counts, 'wt') as f_umi:
        
        line = f_in.readline().strip().split(',')
        if len(line) < 3:
            return
        
        # Write initial header
        header = ['Cell','UMI','Read Type','Count'] + line[3:]
        umi_counts.write(','.join(header) + '\n')
        
        cur_cell, cur_umi, cur_read_type = line[:3]
        cur_umi_count = 1
        while len(line) >= 3:
            cell, umi, read_type = line[:3]
            #print("UMI: %s. Length: %d" % (umi, len(umi)))
            
            if umi != cur_umi or or cell != cur_cell or read_type != cur_read_type:
                list_to_join = [cur_cell, cur_umi, cur_read_type, str(cur_umi_count)] + line[3:]
                umi_counts.write(','.join(list_to_join) + '\n')
                cur_umi_count = 1
                
            else:
                cur_umi_count += 1
                
            line = f_in.readline().strip().split(',')
                

def removeDuplicatesMulti(in_folder, out_folder, cores):
    files = [os.path.splitext(f)[0] for f in os.listdir(in_folder)]
    p = Pool(processes = int(cores))
    func = partial(removeDuplicates, in_folder=in_folder, out_folder=out_folder)
    p.map(func, files)
    p.close()
    p.join()
                

if __name__ == '__main__':
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    cores = sys.argv[3]
    
    removeDuplicatesMulti(input_folder, output_folder, cores)
    
