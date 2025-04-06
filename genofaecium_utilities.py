import os
import subprocess



def check_genome_fasta_sanity(fasta_file):
    genome_exist = False
    genome_nseq = 0
    genome_size = 0
    if os.path.isfile(fasta_file):
        genome_exist = True
        fr = open(fasta_file, 'r')
        line = fr.readline()
        while line != '':
            if line.strip().startswith('>'):
                seqid = line.strip()[1:].split(' ')[0]
                line = fr.readline()
                seqlines = []
                while line != '':
                    seqlines.append(line.strip())
                    line = fr.readline()
                    if line.strip().startswith('>'):
                        break
                seq = ''.join(seqlines)
                genome_nseq += 1
                genome_size += len(seq)
            else:
                line = fr.readline()
        fr.close()
    return [genome_exist, genome_nseq, genome_size]



def run_identification_by_fastani(fastani_executable, input_genome_fasta, species_id_db_dir, species_id_mapping_file, output_tmp_fastani, threads_str):
    # species_id_db_dir = os.path.join(tool_basedir, "genofaecium_db_pre_compiled", "species_id_ref")()
    # species_id_mapping_file = os.path.join(tool_basedir, "genofaecium_db_pre_compiled", "species_id_ref.acc_2_spp.map")()
    # fastANI --rl species_id_ref.fofn -q test_genome/GCA_047261185.1.faecalis.fasta -t 4 -o test_genome/GCA_047261185.1.faecalis.fasta.output/GCA_047261185.1.faecalis.species_ani
    print("[GenoFaecium] [run_identification_by_fastani] create reference genome fofn")
    dict_acc_species = {}
    tmp_fofn = output_tmp_fastani + ".tmp_fofn"
    fw = open(tmp_fofn, 'w')
    fr = open(species_id_mapping_file, 'r')
    line = fr.readline()
    for line in fr:
        accession = line.strip().split("\t")[0]
        species = line.strip().split("\t")[1]
        fasta_path = os.path.join(species_id_db_dir, accession + ".fasta")
        dict_acc_species[accession] = species
        fw.write(fasta_path + "\n")
    fr.close()
    fw.close()
    
    print("[GenoFaecium] [run_identification_by_fastani] run fastANI")
    run_command = [fastani_executable, '--rl', tmp_fofn, '-q', input_genome_fasta, '-t', threads_str, '-o', output_tmp_fastani]
    subprocess.run(run_command)
    os.remove(tmp_fofn)
    
    return dict_acc_species



def read_out_the_identity_call(fastani_output_file, dict_ref_acc_species):
    species_id_value = 'Unidentified'
    ''' the fastANI output file looks like this:
    test_genome_temp/GCA_047613505.1.faecium.fasta  /home/kihyunee/GenoFaecium/genofaecium_db_pre_compiled/species_id_ref/GCA_001544255.1.fasta     99.6669 803     840
    test_genome_temp/GCA_047613505.1.faecium.fasta  /home/kihyunee/GenoFaecium/genofaecium_db_pre_compiled/species_id_ref/GCA_015751045.1.fasta     94.4409 742     840
    test_genome_temp/GCA_047613505.1.faecium.fasta  /home/kihyunee/GenoFaecium/genofaecium_db_pre_compiled/species_id_ref/GCA_001544215.1.fasta     80.4868 390     840
    test_genome_temp/GCA_047613505.1.faecium.fasta  /home/kihyunee/GenoFaecium/genofaecium_db_pre_compiled/species_id_ref/GCA_000271405.2.fasta     80.195  377     840
    test_genome_temp/GCA_047613505.1.faecium.fasta  /home/kihyunee/GenoFaecium/genofaecium_db_pre_compiled/species_id_ref/GCA_000393935.1.fasta     80.148  357     840
    
    find the refernece genome accession at column 2, identity at column 3
    '''
    best_hit_ani = 0
    best_hit_acc = ''
    best_hit_sp = ''
    fr = open(fastani_output_file, 'r')
    for line in fr:
        ls = line.strip().split("\t")
        ref_genome_path = ls[1]
        ani = float(ls[2])
        if ani > best_hit_ani:
            hit_filename = os.path.basename(ref_genome_path)
            acc = hit_filename[:hit_filename.rfind('.fasta')]
            sp = dict_ref_acc_species[acc]
            best_hit_sp = sp
            best_hit_acc = acc
            best_hit_ani = ani
    fr.close()
    if best_hit_ani >= 95:
        species_id_value = best_hit_sp
    return species_id_value


def read_out_the_mlst_call(mlst_output_file):
    mlst_value = 'NA'
    ''' mlst output is like this:
    test_genome_temp/GCA_047613505.1.faecium.fasta  efaecium        32      atpA(3) ddl(3)  gdh(1)  purK(2) gyd(1)  pstS(1) adk(1)
    '''
    fr = open(mlst_output_file, 'r')
    result_line = fr.readline()
    fr.close()
    result_split = result_line.strip().split("\t")
    scheme_used = result_split[1]
    st_called = result_split[2]
    if scheme_used == 'efaecium':
        if st_called.isdigit():
            mlst_value = st_called
    return mlst_value
