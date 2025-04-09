import argparse
import os
import sys
import subprocess
from datetime import datetime
from genofaecium_utilities import run_identification_by_fastani, read_out_the_identity_call, read_out_the_mlst_call



# Time log (start)
start_time = datetime.now()
print("GenoFaecium started at:", start_time.strftime("%Y-%m-%d %H:%M:%S"))


# Add the script's directory to sys.path so local modules can be imported
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, script_dir)

# Some fixed reference files are in the data directory
data_dir = os.path.join(script_dir, 'data')
amr_rule_amr_list_file = os.path.join(data_dir, "amr_rule_amr_list.tsv")
amr_rule_amr_genotype_file = os.path.join(data_dir, "amr_rule_amr_genotype.tsv")
vf_rule_vf_list_file = os.path.join(data_dir, "vf_rule_vf_list.tsv")
vf_rule_vf_genotype_file = os.path.join(data_dir, "vf_rule_vf_genotype.tsv")


# Then import some functions from other scripts
from genofaecium_utilities import check_genome_fasta_sanity
from genofaecium_run_external import run_shell_command_tool_check, run_tool_check_in_env, run_tool_job_in_env


""" GenoFaecium workflows
(1) annotate 
-   use this script (genofaecium.py)
-   perform identity confirmation + MLST + AMR + VF
(2) core_snp
-   use another script (genofaecium_cgsnp.py)
-   given set of genome sequences, 
-   perform SNP variant calling per each input genome based on alignment against a specific reference genome (ST-specific)
-   extract core sites and core SNPs alignment per isolate set
-   phylogenetic inference
"""


""" Parse input arguments 
1. fasta file
2. GenoFaecium directory path
    GenoFaecium/
        conda_packages/
            abricate-1.0.1  amrfinder-3.12.8  mlst-2.23.0  prokka-1.14.6
        dependency_binary/
            fastANI minimap2
        genofaecium_db_pre_compiled:
            abricate_db  amrfinder_db  mlst_db  species_id_ref  species_id_ref.acc_2_spp.map
3. output directory
4. threads to use
"""
parser = argparse.ArgumentParser(description="Extract epidemiologically important information from Enterococcus faecalis genome assembly.")
parser.add_argument("--fasta", dest = "input_genome_fasta", required=True, type=str, help="Path to the input genome assembly fasta file")
parser.add_argument("--out", dest = "output_prefix", required=True, type=str, help="Path to the output prefix (PREFIX.result.txt and PREFIX.files/ will be generated)")
parser.add_argument("--tool_dir", dest = "tool_basedir", required=False, type=str, default='_NS', help="(default = auto-detection) Path to the base directory of GenoFaecium installation (e.g., '/home/user/GenoFaecium');  Under this directory path, you are expected to have 'conda_packages' subdirectory, 'dependency_binary' subdirectory, and 'genofaecium_db_pre_compiled' subdirectory")
parser.add_argument("--threads", dest = "threads_str", required=False, type=str, default="2", help="(default = 2) Number of threads to use")
parser.add_argument("--sample", dest = "sample_name", required=False, type=str, default = "_NS", help="(default = input fasta file name minus .fasta) Sample name to be written in the first column of the output file")
args = parser.parse_args()


input_genome_fasta = args.input_genome_fasta
tool_basedir = args.tool_basedir
output_prefix = args.output_prefix
sample_name_str = args.sample_name
threads_str = args.threads_str

if tool_basedir == '_NS':
    tool_basedir = script_dir


# Check if the necessary contents exist in the tool base directory 
# - Tools
basepack_install_dir = os.path.join(tool_basedir, "conda_packages", "genofaecium_base")
prokka_install_dir = os.path.join(tool_basedir, "conda_packages", "prokka-1.14.6")
amrfinder_install_dir = os.path.join(tool_basedir, "conda_packages", "amrfinder-3.12.8")
abricate_install_dir = os.path.join(tool_basedir, "conda_packages", "abricate-1.0.1")
mlst_install_dir = os.path.join(tool_basedir, "conda_packages", "mlst-2.23.0")
# - Database files
abricate_db_dir = os.path.join(tool_basedir, "genofaecium_db_pre_compiled", "abricate_db")
abricate_db_vfdb_dir = os.path.join(tool_basedir, "genofaecium_db_pre_compiled", "abricate_db", "vfdb")
amrfinder_db_dir = os.path.join(tool_basedir, "genofaecium_db_pre_compiled", "amrfinder_db")
amrfinder_db_vspec_dir = os.path.join(tool_basedir, "genofaecium_db_pre_compiled", "amrfinder_db", "2024-01-31.1")
mlst_db_dir = os.path.join(tool_basedir, "genofaecium_db_pre_compiled", "mlst_db")
mlst_db_efaecium_dir = os.path.join(mlst_db_dir, "pubmlst", "efaecium")
mlst_db_blastdb = os.path.join(mlst_db_dir, "blast", "mlst.fa")
mlst_db_pubmlstdir = os.path.join(mlst_db_dir, "pubmlst")
species_id_db_dir = os.path.join(tool_basedir, "genofaecium_db_pre_compiled", "species_id_ref")
species_id_mapping_file = os.path.join(tool_basedir, "genofaecium_db_pre_compiled", "species_id_ref.acc_2_spp.map")


# Check if dependencies are called and then collect the version information from each tool.
# 1) fastANI in the 'genofaecium_base' env
fastani_exist, fastani_version = run_tool_check_in_env(basepack_install_dir, "fastANI -h")
if fastani_exist:
    fastani_version = 'secret'
# 2) minimap2 in the 'genofaecium_base' env
minimap_exist, minimap_version = run_tool_check_in_env(basepack_install_dir, "minimap2 --version")
# 3) samtools in the 'genofaecium_base' env
samtools_exist, samtools_version = run_tool_check_in_env(basepack_install_dir, "samtools --version")
# 4) bcftools in the 'genofaecium_base' env
bcftools_exist, bcftools_version = run_tool_check_in_env(basepack_install_dir, "bcftools --version")

# 5) prokka in its own environment
prokka_exist, prokka_version = run_tool_check_in_env(prokka_install_dir, "prokka --version")
# 6) amrfinder in its own environment
amrfinder_exist, amrfinder_version = run_tool_check_in_env(amrfinder_install_dir, "amrfinder --version")
# 7) abricate in its own environment
abricate_exist, abricate_version = run_tool_check_in_env(abricate_install_dir, "abricate --version")
# 8) mlst in its own environment
mlst_exist, mlst_version = run_tool_check_in_env(mlst_install_dir, "mlst --version")


# Pass/Fail + What is wrong
dependencies_work = True
what_are_wrong = []
if not os.path.isdir(abricate_db_vfdb_dir):
    dependencies_work = False
    what_are_wrong.append("Abricate VFDB is not found at the expected path: " + abricate_db_vfdb_dir)
if not os.path.isdir(amrfinder_db_vspec_dir):
    dependencies_work = False
    what_are_wrong.append("AMRFinder specific version DB is not found at the expected path: " + amrfinder_db_vspec_dir)
if not os.path.isfile(mlst_db_blastdb):
    dependencies_work = False
    what_are_wrong.append("mlst DB for blastn is not found at the expected path: " + mlst_db_blastdb)
if not os.path.isdir(mlst_db_efaecium_dir):
    dependencies_work = False
    what_are_wrong.append("mlst DB for efaecium is not found at the expected path: " + mlst_db_efaecium_dir)
if not os.path.isdir(species_id_db_dir):
    dependencies_work = False
    what_are_wrong.append("custom DB for enterococcus species ID is not found at the expected path: " + species_id_db_dir)
if not os.path.isfile(species_id_mapping_file):
    dependencies_work = False
    what_are_wrong.append("custom mapping file for enterococcus species ID is not found at the expected path: " + species_id_mapping_file)
if not os.path.isfile(amr_rule_amr_list_file):
    dependencies_work = False
    what_are_wrong.append("the package's mandatory data file for amr reporting rule amr list file not found at the expected path: " + amr_rule_amr_list_file)
if not os.path.isfile(amr_rule_amr_genotype_file):
    dependencies_work = False
    what_are_wrong.append("the package's mandatory data file for amr reporting rule amr-genotype mapping file not found at the expected path: " + amr_rule_amr_genotype_file)
if not os.path.isfile(vf_rule_vf_list_file):
    dependencies_work = False
    what_are_wrong.append("the package's mandatory data file for VF reporting rule VF list file not found at the expected path: " + vf_rule_vf_list_file)
if not os.path.isfile(vf_rule_vf_genotype_file):
    dependencies_work = False
    what_are_wrong.append("the package's mandatory data file for VF reporting rule VF-genotype mapping file not found at the expected path: " + vf_rule_vf_genotype_file)
if fastani_exist:
    print("fastANI ... " + fastani_version + " (OK)")
else:
    #dependencies_work = False
    what_are_wrong.append("fastANI not found at the expected environment: " + basepack_install_dir)
if minimap_exist:
    print("minimap2 ... " + minimap_version + " (OK)")
else:
    dependencies_work = False
    what_are_wrong.append("minimap2 not found at the expected environment: " + basepack_install_dir)
if samtools_exist:
    print("samtools ... " + samtools_version + " (OK)")
else:
    dependencies_work = False
    what_are_wrong.append("samtools not found at the expected environment: " + basepack_install_dir)
if bcftools_exist:
    print("bcftools ... " + bcftools_version + " (OK)")
else:
    dependencies_work = False
    what_are_wrong.append("bcftools not found at the expected environment: " + basepack_install_dir)
if prokka_exist:
    print("prokka ... " + prokka_version + " (OK)")
else:
    dependencies_work = False
    what_are_wrong.append("prokka not found at the expected environment: " + prokka_install_dir)
if amrfinder_exist:
    print("amrfinder ... " + amrfinder_version + " (OK)")
else:
    dependencies_work = False
    what_are_wrong.append("amrfinder not found at the expected environment: " + amrfinder_install_dir)
if abricate_exist:
    print("prokka ... " + abricate_version + " (OK)")
else:
    dependencies_work = False
    what_are_wrong.append("abricate not found at the expected environment: " + abricate_install_dir)
if mlst_exist:
    print("prokka ... " + mlst_version + " (OK)")
else:
    dependencies_work = False
    what_are_wrong.append("mlst not found at the expected environment: " + mlst_install_dir)

if not dependencies_work:
    print("Dependencies are not working.")
    print("----------------------------------")
    print("\n".join(what_are_wrong))
    quit()


# check if the input fasta file makes sense
[input_genome_exist, input_genome_nseq, input_genome_size] = check_genome_fasta_sanity(input_genome_fasta)
if input_genome_exist:
    print("Your input genome:")
    print("- n. seq = " + str(input_genome_nseq))
    print("- size = " + str(input_genome_size) + " bp")
else:
    print("Your input genome not found at: " + input_genome_fasta)
    quit()


# fix the sample name
sample_name = sample_name_str
input_fasta_basename = os.path.basename(input_genome_fasta)
sample_file_basename = input_fasta_basename[0:input_fasta_basename.rfind('.f')]
if sample_name_str == "_NS":
    sample_name = sample_file_basename

'''
# set outputs
prepare output directory
    mkdir output_dir
'''
output_basedir = os.path.dirname(output_prefix)
if not os.path.isdir(output_basedir):
    os.makedirs(output_basedir)
output_tmp_files_dir = output_prefix + '.files'
if not os.path.isdir(output_tmp_files_dir):
    os.makedirs(output_tmp_files_dir)

output_result_row_txt = output_prefix + '.result_row.txt'
output_result_col_txt = output_prefix + '.result_col.txt'
output_tmp_fastani = os.path.join(output_tmp_files_dir, sample_file_basename + '.fastANI')
output_tmp_prokka_dir = os.path.join(output_tmp_files_dir, sample_file_basename + '.prokka_dir')
output_tmp_amrfinder = os.path.join(output_tmp_files_dir, sample_file_basename + '.amrfinder')
output_tmp_abricate_vf = os.path.join(output_tmp_files_dir, sample_file_basename + '.abricate_vf')
output_tmp_mlst = os.path.join(output_tmp_files_dir, sample_file_basename + '.mlst')


reporting_fields = []
reporting_contents = []


# 1. species ref fastANI 1.33
#   create fofn of ref genome paths -- contain all paths under the pre-compiled DB directory's 'species_id_ref' directory's .fasta sub file paths
#   run fastANI
#   output: test_genome/GCA_047261185.1.faecalis.fasta.output/GCA_047261185.1.faecalis.species_ani
dict_ref_acc_species = run_identification_by_fastani(basepack_install_dir, input_genome_fasta, species_id_db_dir, species_id_mapping_file, output_tmp_fastani, threads_str)

# 2. mlst 2.23.0
#   run mlst with --scheme efaecium
print("[GenoFaecium] Run mlst")
mlst_command = f"mlst --blastdb {mlst_db_blastdb} --datadir {mlst_db_pubmlstdir} --scheme efaecium --threads {threads_str} {input_genome_fasta} > {output_tmp_mlst}"
print("[GenoFaecium] Run mlst")
mlst_run_ok, mlst_error_message = run_tool_job_in_env(mlst_install_dir, mlst_command)

# Decide species ID
''' Consider that we need these values in the output columns
column 1    sample
column 2    n contig
column 3    size bp
column 4    species ['Enterococcus faecium', ..., 'Unidentified']
column 5    mlst    ['ST80', ..., 'new'] or 'NA' if species != 'Enterococcus faecium'
'''
species_id_value = read_out_the_identity_call(output_tmp_fastani, dict_ref_acc_species)
mlst_st_value = read_out_the_mlst_call(output_tmp_mlst)

# Record species ID
reporting_contents.append(sample_name)  #[0]
reporting_contents.append(str(input_genome_nseq))  #[1]
reporting_contents.append(str(input_genome_size))  #[2]
reporting_contents.append(species_id_value)  #[3]
reporting_contents.append(mlst_st_value)  #[4]

reporting_fields.append("sample")  #[0]
reporting_fields.append("n.contig")  #[1]
reporting_fields.append("genome.size")  #[2]
reporting_fields.append("species")  #[3]
reporting_fields.append("MLST")  #[4]


# 3. gene prediction by prokka 1.14.6
# run prokka
prokka_command = f"prokka --outdir {output_tmp_prokka_dir} --force --prefix {sample_file_basename} --cpus {threads_str} --noanno --quiet {input_genome_fasta}"
print("[GenoFaecium] Run prokka")
prokka_run_ok, prokka_error_message = run_tool_job_in_env(prokka_install_dir, prokka_command)
print("[GenoFaecium] Clean prokka unnecessary outputs")
os.remove(os.path.join(output_tmp_prokka_dir, sample_file_basename + ".err"))
os.remove(os.path.join(output_tmp_prokka_dir, sample_file_basename + ".gbk"))
os.remove(os.path.join(output_tmp_prokka_dir, sample_file_basename + ".sqn"))
os.remove(os.path.join(output_tmp_prokka_dir, sample_file_basename + ".tbl"))
os.remove(os.path.join(output_tmp_prokka_dir, sample_file_basename + ".tsv"))
os.remove(os.path.join(output_tmp_prokka_dir, sample_file_basename + ".txt"))
prokka_fna = os.path.join(output_tmp_prokka_dir, sample_file_basename + ".fna")
prokka_faa = os.path.join(output_tmp_prokka_dir, sample_file_basename + ".faa")
prokka_gff = os.path.join(output_tmp_prokka_dir, sample_file_basename + ".gff")


# 4. amrfinder 3.12.8
# run amrfinder
amrfinder_command = f"amrfinder -n {prokka_fna} --organism Enterococcus_faecium --plus -o {output_tmp_amrfinder} --threads {threads_str} -p {prokka_faa} -g {prokka_gff} -a prokka --database {amrfinder_db_vspec_dir}"
print("[GenoFaecium] Run amrfinder")
amrfinder_run_ok, amrfinder_error_message = run_tool_job_in_env(amrfinder_install_dir, amrfinder_command)

# prepare the reporting rule 
list_amr_field = []
dict_amr_field_index = {}
buff_amr_field_index = -1
fr = open(amr_rule_amr_list_file, 'r')
line = fr.readline()
for line in fr:
    amr_field = line.strip().split("\t")[1]
    if amr_field not in dict_amr_field_index:
        buff_amr_field_index += 1
        dict_amr_field_index[amr_field] = buff_amr_field_index
        list_amr_field.append(amr_field)
fr.close()
n_amr_field = len(list_amr_field)

dict_genotype_amr_field_index = {}
fr = open(amr_rule_amr_genotype_file, 'r')
line = fr.readline()
for line in fr:
    ls = line.strip().split("\t")
    amr_field = ls[0]
    genotype = ls[1]
    if amr_field not in dict_amr_field_index:
        continue
    amr_field_index = dict_amr_field_index[amr_field]
    dict_genotype_amr_field_index[genotype] = amr_field_index
fr.close()

# Find VRE-relevant AMR genotypes
list_amr_field_n_discovered_genotype = [0]*n_amr_field
list_amr_field_list_discovered_genotype = []
for amr_field_index in range(n_amr_field):
    reporting_fields.append(list_amr_field[amr_field_index])
    list_amr_field_list_discovered_genotype.append([])
fr = open(output_tmp_amrfinder, 'r')
line = fr.readline()
for line in fr:
    ls = line.strip().split("\t")
    genotype = ls[5]
    if genotype in dict_genotype_amr_field_index:
        amr_field_index = dict_genotype_amr_field_index[genotype]
        list_amr_field_n_discovered_genotype[amr_field_index] += 1
        list_amr_field_list_discovered_genotype[amr_field_index].append(genotype)
fr.close()
for amr_field_index in range(n_amr_field):
    n_genotype = list_amr_field_n_discovered_genotype[amr_field_index]
    if n_genotype == 0:
        reporting_contents.append('.')
    else:
        reporting_contents.append(','.join(list_amr_field_list_discovered_genotype[amr_field_index]))



# 5. abricate 1.0.1 with VFDB 
# run abricate - vfdb
abricate_command = f"abricate --db vfdb --threads {threads_str} --datadir {abricate_db_dir} -nopath --noheader {input_genome_fasta} > {output_tmp_abricate_vf}"
print("[GenoFaecium] Run abricate")
abricate_run_ok, abricate_error_message = run_tool_job_in_env(abricate_install_dir, abricate_command)

# prepare the reporting rule 
list_vf_field = []
dict_vf_field_index = {}
buff_vf_field_index = -1
fr = open(vf_rule_vf_list_file, 'r')
line = fr.readline()
for line in fr:
    vf_field = line.strip().split("\t")[1]
    if vf_field not in dict_vf_field_index:
        buff_vf_field_index += 1
        dict_vf_field_index[vf_field] = buff_vf_field_index
        list_vf_field.append(vf_field)
fr.close()
n_vf_field = len(list_vf_field)

dict_genotype_vf_field_index = {}
fr = open(vf_rule_vf_genotype_file, 'r')
line = fr.readline()
for line in fr:
    ls = line.strip().split("\t")
    vf_field = ls[0]
    genotype = ls[1]
    if vf_field not in dict_vf_field_index:
        continue
    vf_field_index = dict_vf_field_index[vf_field]
    dict_genotype_vf_field_index[genotype] = vf_field_index
fr.close()

# Find VRE-relevant VF genotypes
list_vf_field_n_discovered_genotype = [0]*n_vf_field
list_vf_field_list_discovered_genotype = []
for vf_field_index in range(n_vf_field):
    reporting_fields.append(list_vf_field[vf_field_index])
    list_vf_field_list_discovered_genotype.append([])
reporting_fields.append("Other VFs")
list_vf_field_n_discovered_genotype.append(0)
list_vf_field_list_discovered_genotype.append([])
fr = open(output_tmp_abricate_vf, 'r')
for line in fr:
    ls = line.strip().split("\t")
    genotype = ls[5]
    if genotype in dict_genotype_vf_field_index:
        vf_field_index = dict_genotype_vf_field_index[genotype]
        list_vf_field_n_discovered_genotype[vf_field_index] += 1
        list_vf_field_list_discovered_genotype[vf_field_index].append(genotype)
    else:
        vf_field_index = n_vf_field
        list_vf_field_n_discovered_genotype[vf_field_index] += 1
        list_vf_field_list_discovered_genotype[vf_field_index].append(genotype)
fr.close()
for vf_field_index in range(n_vf_field + 1):
    n_genotype = list_vf_field_n_discovered_genotype[vf_field_index]
    if n_genotype == 0:
        reporting_contents.append('.')
    else:
        reporting_contents.append(','.join(list_vf_field_list_discovered_genotype[vf_field_index]))




# Write output
fw = open(output_result_row_txt, 'w')
fw.write(sample_name)
for result_slot_index in range(1, len(reporting_contents)):
    fw.write("\t" + reporting_contents[result_slot_index])
fw.write("\n")
fw.close()

n_total_field = len(reporting_fields)
fw = open(output_result_col_txt, 'w')
fw.write("field\t" + sample_file_basename + "\n")
for result_slot_index in range(n_total_field):
    fw.write(reporting_fields[result_slot_index] + "\t" + reporting_contents[result_slot_index] + "\n")
fw.close()


# Time log (end)
end_time = datetime.now()
print("GenoFaecium ended at:", end_time.strftime("%Y-%m-%d %H:%M:%S"))

# Calculate the difference in seconds
elapsed_time = (end_time - start_time).total_seconds()
print(f"Total execution time: {elapsed_time} seconds")


print("        ".join(reporting_contents))

