import argparse
import subprocess
import os

def run_command(command):
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        exit(1)

def extract_single_cov_aligned_interval(samtools_depth_file):
    # requiring the depth file to be run with "-aa" option of samtools depth, (samtools depth -aa ...)
    # so that absolutely all positions are recorded in the file, no matter alignment was present or not
    list_coordinate = []
    # [[chr, start, end]]
    buff_scaln_chr = ''
    buff_scaln_start = -1
    buff_scaln_end = -1
    fr = open(samtools_depth_file, 'r')
    for line in fr:
        ls = line.strip().split("\t")
        chr = ls[0]
        pos = int(ls[1])
        cov = int(ls[2])
        if cov == 1:
            # start the interval, or extend
            if buff_scaln_start == -1:
                buff_scaln_chr = chr
                buff_scaln_start = pos
                buff_scaln_end = pos
            else:
                buff_scaln_end = pos
        else:
            # keep it quiet, or report the aligned coordinate that just ended
            if buff_scaln_start != -1:
                scaln_coord = [buff_scaln_chr, buff_scaln_start, buff_scaln_end]
                list_coordinate.append(scaln_coord)
                buff_scaln_chr = ''
                buff_scaln_start = -1
                buff_scaln_end = -1
    if buff_scaln_start != -1:
        scaln_coord = [buff_scaln_chr, buff_scaln_start, buff_scaln_end]
        list_coordinate.append(scaln_coord)
    fr.close()
    return list_coordinate

def extract_snps_simple(bcftools_snp_only_vcf_file):
    list_snp = []
    # [[chr, 1-based index, allele]]
    fr = open(bcftools_snp_only_vcf_file, 'r')
    for line in fr:
        if line.strip().startswith('#'):
            continue
        ls = line.strip().split("\t")
        chr = ls[0]
        pos = int(ls[1])
        alt_expr = ls[4]
        allele = alt_expr.split(',')[0]
        new_snp = [chr, pos, allele]
        list_snp.append(new_snp)
    fr.close()
    return list_snp

def write_align_interval(list_scaln_interval, outfile):
    fw = open(outfile, 'w')
    num_coord = len(list_scaln_interval)
    for aln_index in range(num_coord):
        aln_coord = list_scaln_interval[aln_index]
        fw.write(aln_coord[0] + "\t" + str(aln_coord[1]) + "\t" + str(aln_coord[2]) + "\n")
    fw.close()

def write_snp_var(list_snp_var, outfile):
    fw = open(outfile, 'w')
    num_snp = len(list_snp_var)
    for snp_index in range(num_snp):
        snp_def = list_snp_var[snp_index]
        fw.write(snp_def[0] + "\t" + str(snp_def[1]) + "\t" + snp_def[2] + "\n")
    fw.close()


def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="Script to handle genome data processing paths.")

    # Add the arguments
    parser.add_argument("--minimap", dest = "minimap_executable_path", required=False, default="minimap2", help="Path to the Minimap2 executable")
    parser.add_argument("--samtools", dest = "samtools_executable_path", required=False, default="samtools", help="Path to the Samtools executable")
    parser.add_argument("--bcftools", dest = "bcftools_executable_path", required=False, default="bcftools", help="Path to the Bcftools executable")
    parser.add_argument("--ref_mmi", dest = "ref_genome_index", required=True, help="Path to the reference genome minimap index file, built with -x asm20")
    parser.add_argument("--ref_fa", dest = "ref_genome_fasta", required=True, help="Path to the reference genome FASTA file")
    parser.add_argument("--input_fa", dest = "input_genome_fasta", required=True, help="Path to the input genome FASTA file")
    parser.add_argument("--out_prefix", dest = "output_prefix", required=True, help="Path to the prefix of output files")
    parser.add_argument("--threads", dest = "num_threads_str", required=False, default="4", help="number of threads to use")

    # Parse the arguments
    args = parser.parse_args()

    # Store the arguments in variables
    minimap_executable_path = args.minimap_executable_path
    samtools_executable_path = args.samtools_executable_path
    bcftools_executable_path = args.bcftools_executable_path
    ref_genome_index = args.ref_genome_index
    ref_genome_fasta = args.ref_genome_fasta
    input_genome_fasta = args.input_genome_fasta
    output_prefix = args.output_prefix
    num_threads_str = args.num_threads_str

    # Create output directory if output prefix contains new directory
    output_dir = os.path.dirname(output_prefix)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    # Define temporary files
    temp_sam = output_prefix + ".temp.align.sam"
    temp_bam = output_prefix + ".temp.align.bam"
    temp_bam_sort = output_prefix + ".temp.align.sorted.bam"
    temp_bam_sort_bai = output_prefix + ".temp.align.sorted.bam.bai"
    temp_depth = output_prefix + ".temp.align.depth"
    temp_bcf = output_prefix + ".temp.align.bcf"
    temp_snp_vcf = output_prefix + ".temp.align.snp.vcf"
    # Define final output files
    final_scaln_coord = output_prefix + ".on_ref_coord.single_cov_aln"
    final_snp = output_prefix + ".on_ref_coord.snp"

    # Commands to execute
    commands = [
        f"{minimap_executable_path} -x asm20 -t {num_threads_str} -o {temp_sam} -a --secondary=no {ref_genome_index} {input_genome_fasta}",
        f"{samtools_executable_path} view -b -o {temp_bam} --threads {num_threads_str} {temp_sam}",
        f"{samtools_executable_path} sort -o {temp_bam_sort} --threads {num_threads_str} {temp_bam}",
        f"{samtools_executable_path} index --threads {num_threads_str} {temp_bam_sort}",
        f"{samtools_executable_path} depth -aa -o {temp_depth} {temp_bam_sort}",
        f"{bcftools_executable_path} mpileup -f {ref_genome_fasta} -O u -o {temp_bcf} --threads {num_threads_str} --skip-indels {temp_bam_sort}",
        f"{bcftools_executable_path} view -o {temp_snp_vcf} -O v --threads {num_threads_str} --types snps {temp_bcf}"
    ]

    # Execute the commands
    for command in commands:
        print(f"Executing: {command}")
        run_command(command)

    print("All commands executed successfully.")

    # Extract the depth and SNP information in a minimal format
    # [[chr, start, end]]
    list_scaln_interval = extract_single_cov_aligned_interval(temp_depth)
    # [[chr, 1-based index, allele]]
    list_snp_var = extract_snps_simple(temp_snp_vcf)

    # Write final output
    write_align_interval(list_scaln_interval, final_scaln_coord)
    write_snp_var(list_snp_var, final_snp)

    os.remove(temp_bam)
    os.remove(temp_bam_sort)
    os.remove(temp_bam_sort_bai)
    os.remove(temp_sam)
    os.remove(temp_bcf)
    os.remove(temp_snp_vcf)
    os.remove(temp_depth)

if __name__ == "__main__":
    main()