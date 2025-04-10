import argparse
import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord




def parse_seqs_from_fasta(fasta_file):
    list_seqid = []
    dict_sequence = {}
    fr = open(fasta_file, 'r')
    line = fr.readline()
    while line != '':
        if line.strip().startswith('>'):
            seqid = line.strip()[1:].split(' ')[0]
            seq_l = []
            line = fr.readline()
            while line != '':
                seq_l.append(line.strip())
                line = fr.readline()
                if line.strip().startswith('>'):
                    break
            list_seqid.append(seqid)
            dict_sequence[seqid] = ''.join(seq_l)
        else:
            line = fr.readline()
    fr.close()
    return [list_seqid, dict_sequence]


def translate_sequences(input_file, output_file, gcode_to_use, stop_at_stop):
    # Read input fasta file
    sequences = SeqIO.parse(input_file, "fasta")
    
    output_records = []
    
    # Translate each sequence and create SeqRecord objects
    for sequence in sequences:
        seq_id = sequence.id
        seq_description = sequence.description
        seq = sequence.seq
        protein_seq = seq.translate(table=int(gcode_to_use), to_stop=stop_at_stop)
        
        # Create SeqRecord object for the translated protein sequence
        protein_record = SeqRecord(protein_seq, id=seq_id, description=seq_description)
        
        # Append the translated protein sequence record to the output list
        output_records.append(protein_record)
    
    # Write the translated protein sequences to the output fasta file
    SeqIO.write(output_records, output_file, "fasta")



dict_valid_ntbase = {}
dict_valid_ntbase['A'] = True;  dict_valid_ntbase['a'] = True
dict_valid_ntbase['T'] = True;  dict_valid_ntbase['t'] = True
dict_valid_ntbase['C'] = True;  dict_valid_ntbase['c'] = True
dict_valid_ntbase['G'] = True;  dict_valid_ntbase['g'] = True
dict_valid_ntbase['U'] = True;  dict_valid_ntbase['u'] = True


# input structure
parser = argparse.ArgumentParser(description = "Input specifications ")
parser.add_argument("--in", dest="input_cds_nt_fasta", required=True, type=str, help="input CDS nt fasta file")
parser.add_argument("--out", dest="output_codon_align_fasta", required=True, type=str, help="output codon alignment fasta file")
parser.add_argument("--muscle", dest="syspath_to_muscle", required=False, default="muscle", type=str, help="path to executable muscle version >= 5;  default = 'muscle' (installed system-wide)")
parser.add_argument("--gcode", dest="genetic_code_number", required=False, default="11", type=str, help="genetic code number as defined by NCBI; for example, bacterial code is '11', standard code is '1', etc. (defaults = 11)")
parser.add_argument("--threads", dest="num_threads_str", required=False, default="2", type=str, help="number of threads to use (default = 2)")


args = parser.parse_args()
input_cds_nt_fasta = args.input_cds_nt_fasta
output_codon_align_fasta = args.output_codon_align_fasta
syspath_to_muscle = args.syspath_to_muscle
genetic_code_number = args.genetic_code_number
num_threads_str = args.num_threads_str



# intermediate files
intermed_prot_fasta = output_codon_align_fasta + ".itm_tmp.prot"
intermed_prot_align = output_codon_align_fasta + ".itm_tmp.prot.aln"


# translate the input CDS nt sequences --> protein sequences
# fasta --> fasta
translate_sequences(input_cds_nt_fasta, intermed_prot_fasta, genetic_code_number, False)


# check which has premature stop codons 


# align the protein sequences
aligner_command = []
aligner_command.append(syspath_to_muscle)
aligner_command.append("-super5")
aligner_command.append(intermed_prot_fasta)
aligner_command.append("-output")
aligner_command.append(intermed_prot_align)
aligner_command.append("-amino")
aligner_command.append("-threads")
aligner_command.append(num_threads_str)

print("[python subprocess] RUN: " + " ".join(aligner_command))
subprocess.run(aligner_command)


# project the gaps made in the amino acid alignments into the nucleotide sequences
[list_raw_nt_id, dict_raw_nt_seq] = parse_seqs_from_fasta(input_cds_nt_fasta)
[list_aligned_aa_id, dict_aligned_aa_seq] = parse_seqs_from_fasta(intermed_prot_align)


# number of sequences should be the same
if len(list_raw_nt_id) != len(list_aligned_aa_id):
    print("Why is there missing sequence in alignment?")
    print("  ...  n. sequences in input = " + str(len(list_raw_nt_id)))
    print("  ...  n. sequences in intermediate aa alignment = " + str(len(list_aligned_aa_id)))
    quit()


length_of_aa_aln = 0
length_of_final_nt_aln = 0
fw = open(output_codon_align_fasta, 'w')
for input_seqid in list_raw_nt_id:
    raw_nt_seq = dict_raw_nt_seq[input_seqid]
    aligned_aa_seq = dict_aligned_aa_seq[input_seqid]
    if len(aligned_aa_seq) > length_of_aa_aln:
        length_of_aa_aln = len(aligned_aa_seq)
        length_of_final_nt_aln = 3*len(aligned_aa_seq)
    
    # record the zero-based index of the <gap insertion sites>, and the length of each gap, both in <aa unit> (not in <nt unit> !!)
    dict_gapstart_zbi_to_gaplength_inaa = {}
    num_aa_in_aln = 0
    buff_physical_aazbi = 0
    buff_curr_gap_startzbi = -1
    buff_curr_gap_length = 0
    for aln_aazbi in range(len(aligned_aa_seq)):
        aaresidue = aligned_aa_seq[aln_aazbi]
        if aaresidue != '-':
            buff_physical_aazbi += 1
            num_aa_in_aln += 1
            # if there was a gap content in the buffer, register it now and initizlize the buffers.
            if buff_curr_gap_length > 0:
                dict_gapstart_zbi_to_gaplength_inaa[buff_curr_gap_startzbi] = buff_curr_gap_length
                buff_curr_gap_startzbi = -1
                buff_curr_gap_length = 0
        else:
            # open or extend a gap
            if buff_curr_gap_length == 0:
                # open a gap
                buff_curr_gap_startzbi = buff_physical_aazbi
                buff_curr_gap_length = 1
            else:
                # extend a gap
                buff_curr_gap_length += 1
    # reached the last site, if there is still a gap info in the buffer, register it
    if buff_curr_gap_length > 0:
        dict_gapstart_zbi_to_gaplength_inaa[buff_curr_gap_startzbi] = buff_curr_gap_length

    # length of the raw nt sequence should be 3 x the number of aa in the aligned aa sequence
    if num_aa_in_aln*3 != len(raw_nt_seq):
        print(" I am uncomfortable by the fact that the number of aa residuesin your aligned sequence (" + str(num_aa_in_aln) + " x 3 is not the number of nt in the CDS (" + str(len(raw_nt_seq)) + ")")
        print("   ...  for " + input_seqid)
        print("   ...  AA: " + aligned_aa_seq)
        print("   ...  nt: " + raw_nt_seq)
        print("")


    dict_gapstart_zbi_to_gaplength_innt = {}
    for aazbi in dict_gapstart_zbi_to_gaplength_inaa:
        ntzbi = aazbi*3
        ntlength = dict_gapstart_zbi_to_gaplength_inaa[aazbi]*3
        dict_gapstart_zbi_to_gaplength_innt[ntzbi] = ntlength

    fw.write(">" + input_seqid + "\n")
    n_writen_sites = 0
    for ntzbi in range(len(raw_nt_seq)):
        ntchar = raw_nt_seq[ntzbi]
        if ntzbi in dict_gapstart_zbi_to_gaplength_innt:
            gaplength = dict_gapstart_zbi_to_gaplength_innt[ntzbi]
            fw.write(''.join(['-']*gaplength))
            n_writen_sites += gaplength
        fw.write(ntchar)
        n_writen_sites += 1
    if n_writen_sites < length_of_final_nt_aln:
        fw.write(''.join(['-']*(length_of_final_nt_aln - n_writen_sites)))
    fw.write("\n")
fw.close()



os.remove(intermed_prot_fasta)
os.remove(intermed_prot_align)

