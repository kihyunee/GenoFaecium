import os
import argparse


def count_diffs(str1, str2):
  """Counts the number of character differences between two equal-length strings.
  Args:
      str1: The first string.
      str2: The second string.
  Returns:
      The number of character differences between the strings.
  """
  if len(str1) != len(str2):
    raise ValueError("Strings must have equal length")
  return sum(a != b for a, b in zip(str1, str2))


parser = argparse.ArgumentParser(description="you have run minimap-based workflow to detect SNPs from each query genome, against a single reference genome - using core_snp_pipeline_call_ref_snp_by_minimap.py")
# Add the arguments
parser.add_argument("--ref", dest = "ref_chromosome_fasta", required=True, help="reference chromosome fasta file - used in the minimap and samtools step")
parser.add_argument("--dir", dest = "snp_call_dir", required=True, help="directory where you have created <ACCESSION>.on_ref_coord.snp and <ACCESSION>.on_ref_coord.single_cov_aln files per each query genome")
parser.add_argument("--qlist", dest = "query_genome_list_file", required=True, help="list of query genome accessions to analyze in this batch")
parser.add_argument("--outdir", dest = "output_dir", required=True, help="where to write output files (fasta and statistics)")
parser.add_argument("--prefix", dest = "output_prefix", required=True, help="prefix for all output files (excluding directory path part)")
parser.add_argument("--stop_at_core", dest = "stop_at_core", required=False, help="optional:  use this flag to stop and exit when you get to count the number of core sites", action="store_true", default=False)
parser.add_argument("--stop_at_snp", dest = "stop_at_snp", required=False, help="optional:  use this flag to stop and exit when you get to count the core SNP sites", action="store_true", default=False)
parser.add_argument("--core_thresh", dest = "core_thresh_pct", required=False, type=float, default = 99, help="percent threshold for prevalence of a site to become a core site (default = 99)")
parser.add_argument("--full_fasta", dest = "write_full_fasta", required=False, default = False, action="store_true", help="give this flag to write full (every ref position) fasta")
parser.add_argument("--dist", dest = "write_distance", required=False, default = False, action = "store_true", help="give this flag to write pairwise distance table")

# Parse the arguments
args = parser.parse_args()
# Store the arguments in variables
ref_chromosome_fasta = args.ref_chromosome_fasta
snp_call_dir = args.snp_call_dir
query_genome_list_file = args.query_genome_list_file
output_dir = args.output_dir
output_prefix = args.output_prefix
stop_at_core = args.stop_at_core
stop_at_snp = args.stop_at_snp
core_thresh_pct = args.core_thresh_pct
core_thresh = core_thresh_pct/100
write_full_fasta = args.write_full_fasta
write_distance = args.write_distance



# Workflow in this script
# 1. define the core-aligned regions - based on the reference chromosome coordinate
# 2. define the sites with SNP polymorphism 
# 3. write fasta output == alleles at the polymorphic sites
# 4. write fasta output == alleles at all reference positions (gap all regions that were not single-copy-coverage-aligned)
# 5. write distance matrix based on cgSNP diff.


if not os.path.isdir(output_dir):
    os.makedirs(output_dir)


# 0. inputs
list_input_genome = []
n_input_genome = 0
fr = open(query_genome_list_file, 'r')
for line in fr:
    list_input_genome.append(line.strip())
fr.close()
n_input_genome = len(list_input_genome)
print("# number of input query genomes = " + str(n_input_genome))


# 1. define the core-aligned regions - based on the reference chromosome coordinate
## read reference chromosome sequence; length and bases
ref_chr_id = ''
ref_chr_length = 0
ref_chr_seq = ''

fr = open(ref_chromosome_fasta, 'r')
line = fr.readline()
while line != '':
    if line.strip().startswith('>'):
        ref_chr_id = line.strip()[1:].split(' ')[0]
        ref_chr_seq_l = []
        line = fr.readline()
        while line != '':
            ref_chr_seq_l.append(line.strip())
            line = fr.readline()
            if line.strip().startswith('>'):
                break
        ref_chr_seq = ''.join(ref_chr_seq_l)
        break
    else:
        line = fr.readline()
fr.close()
ref_chr_length = len(ref_chr_seq)
print("# reference chromosoem length = " + str(ref_chr_length))

## define core sites
site_cov_count_vector = [0]*ref_chr_length
print("# count alignment coverage per reference site")
for qidx in range(n_input_genome):
    q_genome = list_input_genome[qidx]
    q_genome_snp_file = os.path.join(snp_call_dir, q_genome + '.on_ref_coord.snp')
    q_genome_cov_file = os.path.join(snp_call_dir, q_genome + '.on_ref_coord.single_cov_aln')
    print("# ... " + str(qidx + 1) + "/" + str(n_input_genome))
    if not os.path.isfile(q_genome_cov_file):
        print("# ... ! file not found: " + q_genome_cov_file + " ... skip this query.  This will have a detrimental effect in core site recovery")
        continue
    fr = open(q_genome_cov_file, 'r')
    for line in fr:
        ls = line.strip().split("\t")
        aln_start_zbi = int(ls[1]) - 1
        aln_end_obi = int(ls[2])
        for site_zbi in range(aln_start_zbi, aln_end_obi):
            site_cov_count_vector[site_zbi] += 1
    fr.close()
    
site_is_core = [False]*ref_chr_length
n_core_site = 0
n_core_site_95 = 0
n_core_site_99 = 0
n_core_site_100 = 0
for site_zbi in range(ref_chr_length):
    site_freq = float(site_cov_count_vector[site_zbi])/float(n_input_genome)
    if site_freq >= core_thresh:
        site_is_core[site_zbi] = True
        n_core_site += 1
    if site_freq >= 1:
        n_core_site_95 += 1
        n_core_site_99 += 1
        n_core_site_100 += 1
    elif site_freq >= 0.99:
        n_core_site_95 += 1
        n_core_site_99 += 1
    elif site_freq >= 0.95:
        n_core_site_95 += 1
print("# number of core sites = " + str(n_core_site) + " | out of total reference chromosome length " + str(ref_chr_length))
print("# ... num core sites if thresh = 100:  " + str(n_core_site_100))
print("# ... num core sites if thresh = 99:  " + str(n_core_site_99))
print("# ... num core sites if thresh = 95:  " + str(n_core_site_95))

if stop_at_core:
    exit()



# 2. define the sites with SNP polymorphism 
site_polymorphism_freq = [0]*ref_chr_length
print("# count polymorphism frequency per reference site")
for qidx in range(n_input_genome):
    q_genome = list_input_genome[qidx]
    q_genome_snp_file = os.path.join(snp_call_dir, q_genome + '.on_ref_coord.snp')
    q_genome_cov_file = os.path.join(snp_call_dir, q_genome + '.on_ref_coord.single_cov_aln')
    print("# ... " + str(qidx + 1) + "/" + str(n_input_genome))
    if not os.path.isfile(q_genome_snp_file):
        print("# ... ! file not found: " + q_genome_snp_file + " ... skip this query")
        continue
    fr = open(q_genome_snp_file, 'r')
    for line in fr:
        ls = line.strip().split("\t")
        site_zbi = int(ls[1]) - 1
        site_polymorphism_freq[site_zbi] += 1
    fr.close()

site_is_coresnp = [False]*ref_chr_length
n_coresnp_site = 0
for site_zbi in range(ref_chr_length):
    if not site_is_core[site_zbi]:
        continue
    if site_polymorphism_freq[site_zbi] > 0:
        site_is_coresnp[site_zbi] = True
        n_coresnp_site += 1
print("# number of core SNP sites = " + str(n_coresnp_site) + " | out of total reference chromosome length " + str(ref_chr_length) + " | total core sites = " + str(n_core_site))

if stop_at_snp:
    exit()



# 3. write fasta output == only writing alleles at the polymorphic sites
output_snp_site_fastas = os.path.join(output_dir, output_prefix + ".snp_site.fasta")
fw = open(output_snp_site_fastas, 'w')
print("# writing core SNP site fasta")
for qidx in range(n_input_genome):
    q_genome = list_input_genome[qidx]
    q_genome_snp_file = os.path.join(snp_call_dir, q_genome + '.on_ref_coord.snp')
    q_genome_cov_file = os.path.join(snp_call_dir, q_genome + '.on_ref_coord.single_cov_aln')
    print("# ... " + str(qidx + 1) + "/" + str(n_input_genome))
    if not os.path.isfile(q_genome_snp_file):
        print("# ... ! file not found: " + q_genome_snp_file + " ... skip this query")
        continue
    
    genome_cov_vector = [False]*ref_chr_length
    #genome_allele_vector = ['-']*n_coresnp_site
    genome_allele_vector = []
    
    fr = open(q_genome_cov_file, 'r')
    for line in fr:
        ls = line.strip().split("\t")
        aln_start_zbi = int(ls[1]) - 1
        aln_end_obi = int(ls[2])
        for site_zbi in range(aln_start_zbi, aln_end_obi):
            genome_cov_vector[site_zbi] = True
    fr.close()
    # i am checking this "single-copy" coverage because I don't want to call a SNP that is detected by multiply aligned region

    dict_sitezbi_snpallele = {}
    fr = open(q_genome_snp_file, 'r')
    for line in fr:
        ls = line.strip().split("\t")
        site_zbi = int(ls[1]) - 1
        dict_sitezbi_snpallele[site_zbi] = ls[2]
    fr.close()

    for site_zbi in range(ref_chr_length):
        if not site_is_coresnp[site_zbi]:
            continue
        allele = '-'
        if genome_cov_vector[site_zbi]:
            allele = ref_chr_seq[site_zbi]
            if site_zbi in dict_sitezbi_snpallele:
                allele = dict_sitezbi_snpallele[site_zbi]
        genome_allele_vector.append(allele)
    
    fw.write(">" + q_genome + "\n" + ''.join(genome_allele_vector) + "\n")
fw.close()




# 4. write fasta output == alleles at all reference positions (gap all regions that were not single-copy-coverage-aligned)
if write_full_fasta:
    print("writing full length fasta?? sorry but I cannot do this for you now.")
    quit()

if not write_distance:
    print("you did not order pairwise distnace calculation ")
    quit()

# 5. write distance matrix based on cgSNP diff.
list_allele_seqs = []
fr = open(output_snp_site_fastas, 'r')
for line in fr:
    if line.strip().startswith('>'):
        continue
    list_allele_seqs.append(line.strip())
fr.close()

print("# calculating pairwise SNP distance")
output_distance_table = os.path.join(output_dir, output_prefix + ".core_snp.distance.tsv")
fw = open(output_distance_table, 'w')
fw.write("genome_a\tgenome_b\tcore_snp_dist\n")
for a_index in range(n_input_genome - 1):
    print("# ... pairwise ... from anchor " + str(a_index + 1) + "/" + str(n_input_genome) + " ... doing " + str(n_input_genome - a_index - 1) + " calculations")
    a_allele_str = list_allele_seqs[a_index]
    for b_index in range(a_index + 1, n_input_genome):
        b_allele_str = list_allele_seqs[b_index]
        coresnp_dist = count_diffs(a_allele_str, b_allele_str)
        fw.write(list_input_genome[a_index] + "\t" + list_input_genome[b_index] + "\t" + str(coresnp_dist) + "\n")
fw.close()



