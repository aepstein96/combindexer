import os
import glob
import gzip
import fire

def getReadNumbers(sample, input_folder, barcoded_folder, umis_folder):
    read_numbers = {'Raw': 0, 'Barcoded': 0, 'Assigned UMI': 0, 'Gene Identified': 0, 'Gene Identified (Unique >=1)': 0, 'Gene Identified (Unique >=2)': 0, 'Gene Identified (Unique >=3)': 0, 'Gene Identified (Unique >=4)': 0, 'Gene Identified (Unique >=5)': 0}
    
    # Counting raw reads
    print("Counting raw reads...")
    os.chdir(input_folder)
    fastq_files = glob.glob(f"{sample}_S*_I1_001.fastq.gz")
    if len(fastq_files) == 0:
        print(f"No fastq files found for {sample}")
        return read_numbers
    
    if len(fastq_files) > 1:
        raise ValueError(f"Multiple fastq files found for {sample}")
        
    fastq_file = fastq_files[0]
    with gzip.open(fastq_file, 'rt') as f:
        line = f.readline()
        num_lines = 0
        while line:
            num_lines += 1
            line = f.readline()
            
        read_numbers['Raw'] = num_lines // 4
    
    # Counting barcoded reads
    print("Counting barcoded reads...")
    barcoded_file = os.path.join(barcoded_folder, f"{sample}_output.csv")
    if not os.path.exists(barcoded_file):
        raise ValueError(f"Barcoded file {barcoded_file} not found for {sample}")
    
    with open(barcoded_file, 'r') as f:
        line = f.readline()
        num_lines = 0
        while line:
            num_lines += 1
            line = f.readline()
            
        read_numbers['Barcoded'] = num_lines - 1 # Subtracting 1 for header
    
    # Counting both total reads assigned UMIs and distinct UMIs
    print("Counting assigned UMIs...")
    umis_file = os.path.join(umis_folder, f"{sample}_output_umi_counts.csv")
    if not os.path.exists(umis_file):
        raise ValueError(f"UMIs file not found for {sample}")
    
    with open(umis_file, 'r') as f:
        line = f.readline().split(',') # Header
        line = f.readline().split(',')
        while len(line) >= 4:
            count = int(line[3])
            read_numbers['Assigned UMI'] += count
            
            # Checking if gene was identified
            RT_primer, _,  SSS_primer, genotype, RT_primer = line[2].strip(';').split(';')
            SSS_gene = SSS_primer.split('_')[0]
            genotype_gene = genotype.split('_')[0]
            if SSS_gene == genotype_gene:
                read_numbers['Gene Identified'] += count
                for i in range(1, 6):
                    if count >= i:
                        read_numbers[f'Gene Identified (Unique >={i})'] += 1
                        
            line = f.readline().split(',')

    return read_numbers

def main(sample, input_folder, barcoded_folder, umis_folder, output_file):does star
    print(f"Counting reads for {sample}...")
    read_numbers = getReadNumbers(sample, input_folder, barcoded_folder, umis_folder)
    with open(output_file, 'a') as f:
        f.write(f"{sample}")
        for key in read_numbers:
            f.write(f",{read_numbers[key]}")
        f.write("\n")

if __name__ == "__main__":
    fire.Fire(main)