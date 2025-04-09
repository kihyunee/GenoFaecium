# 1. GenoFaecium
Genomic analysis framework for Enterococcus faecium clinical isolates.


# 2. Installation
## 2-1. Linux
This installation guide is based on Linux OS  (such as windows subsystem linux).

I tested this installation process in Windows Ubuntu 24.04.1 LTS.


## 2-2. Miniconda
If you don't have conda in your system, install a miniconda according to [Miniconda Quickstart Installation Instruction](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions).

You may try the following commands (copied from the aforementioned link).

```
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
```

Once miniconda installation is finished, do the following commands.

```
source ~/miniconda3/bin/activate
conda init --all
```

Finally, close and re-open the linux once, before you go further.


## 2-3. Clone the GenoFaecium directory in your system.
Decide the location to create the genofaecium directory.

I will use my home directory (`/home/kihyun/`) which is `~/`

_**If you are not using the home directory, please replace all `~` in the paths used in the following installation scripts to your corresponding base directory path.**_

```
# Go to the location you decided
cd ~

# Clone the repository
git clone https://github.com/kihyunee/GenoFaecium.git
```


## 2-4. Download the dependencies from zenodo file URL and unpack it.

We will download a tar file containing some conda packages and database files that GenoFaecium depends on. Inside the tar:
* conda_packages/abricate-1.0.1/abricate-1.0.1.tar.gz
* conda_packages/amrfinder-3.12.8/amrfinder-3.12.8.tar.gz
* conda_packages/genofaecium_base/genofaecium_base.tar.gz
* conda_packages/mlst-2.23.0/mlst-2.23.0.tar.gz
* conda_packages/prokka-1.14.6/prokka-1.14.6.tar.gz
* genofaecium_db_pre_compiled.tar.gz
* test_genome/GCA_047261185.1.faecalis.fasta
* test_genome/GCA_047613505.1.faecium.fasta

This tar file **must be unpacked in the `GenoFaecium` directory** that you just created in the step **2-3** using git clone.

```
# Go to the GenoFaecium directory and download the tar file from zenodo URL
cd GenoFecium
wget https://zenodo.org/records/15183043/files/genofaecium_dependencies.tar --no-check-certificate
```

You have to first unpack the `tar` and then subsequently unpack each of the contained `.tar.gz` files. Use the following commands.

```
# Unpack the master tar first
tar -xf genofaecium_dependencies.tar

# Unpack the pre-compiled database directory
tar -xzf genofaecium_db_pre_compiled.tar.gz

# Unpack each conda package inside the ./conda_package/ directory
cd conda_packages/abricate-1.0.1/
tar -xzf abricate-1.0.1.tar.gz
cd ../amrfinder-3.12.8/
tar -xzf amrfinder-3.12.8.tar.gz
cd ../genofaecium_base/
tar -xzf genofaecium_base.tar.gz
cd ../mlst-2.23.0/
tar -xzf mlst-2.23.0.tar.gz
cd ../prokka-1.14.6/
tar -xzf prokka-1.14.6.tar.gz
cd ../../
```


# 3. Execute GenoFaecium

```
python genofaecium.py -h

usage: genofaecium.py [-h] --fasta INPUT_GENOME_FASTA --out OUTPUT_PREFIX --tool_dir TOOL_BASEDIR [--threads THREADS_STR] [--sample SAMPLE_NAME]

options:
  -h, --help            show this help message and exit
  --fasta INPUT_GENOME_FASTA
                        Path to the input genome assembly fasta file
  --out OUTPUT_PREFIX   Path to the output prefix (PREFIX.result.txt and PREFIX.files/ will be generated)
  --tool_dir TOOL_BASEDIR
                        Path to the base directory of GenoFaecium installation (e.g., '/home/user/GenoFaecium'); Under this directory path, you are
                        expected to have 'conda_packages' subdirectory, 'dependency_binary' subdirectory, and 'genofaecium_db_pre_compiled' subdirectory
  --threads THREADS_STR
                        Number of threads to use (default = 2)
  --sample SAMPLE_NAME  sample name to be written in the first column of the output file; default = input fasta file name minus fasta
```

For example, to you will find an E. faecium genome fasta file at `test_genome/GCA_047613505.1.faecium.fasta` in the tool's directory.

Analyze that genome sequence:

```
python genofaecium.py --fasta test_input/GCA_047785525.1.fasta --out test_output/GCA_047785525.1 --tool_dir /home/kihyunee/GenoFaecium
```

The output files 
## prefix.result_col.txt
- Is a two-column tab-delimited table of field names (column 1) and report values (column 2).
- The order or fields is consistent across all output files.
- Example output content:
```
field   GCA_047785525.1
sample  GCA_047785525.1
n.contig        296
genome.size     3079105
species Enterococcus faecium
MLST    80
Aminoglycoside  aac(6')-I,aac(6')-Ie/aph(2'')-Ia,aph(3')-IIIa,ant(6)-Ia
Penicillin pbp5_M485A   pbp5_M485A
Penicillin pbp5_M485T   .
Beta-lactamase  .
vanA    vanA
vanB    .
vanD    .
vancomycin-others       .
Oxazolidinone 23S mutation      .
Oxazolidinone Cfr       .
Oxazolidinone OptrA     .
Oxazolidinone PoxtA     .
Tetracycline    tet(S),tet(L)
Daptomycin      liaR_W73C,liaS_T120A
Rifamycin       .
Fosfomycin      .
Streptogramin   .
Phenicol        .
Quinolone       gyrA_S83I,parC_S80I
Macrolides      msr(C),erm(B)
MSCRAMM acm     acm
Biofilm sgrA    sgrA
MSCRAMM fss3    fss3
MSCRAMM ecbA    .
MSCRAMM scm     scm
Biofilm esp     .
Other VFs       .
```

## preifx.result_row.txt
- Is a single-row text file, tab-delimited list of report values ordered as the rows presented in the above column-format output file.
- This is made for you to easily concatenate the reports generated for many isolates into a single table using `cat`.
- Example output content:
```
GCA_047785525.1 296     3079105 Enterococcus faecium    80      aac(6')-I,aac(6')-Ie/aph(2'')-Ia,aph(3')-IIIa,ant(6)-Ia pbp5_M485A      .       .       vanA.       .       .       .       .       .       .       tet(S),tet(L)   liaR_W73C,liaS_T120A    .       .       .       .       gyrA_S83I,parC_S80I     msr(C),erm(B)       acm     sgrA    fss3    .       scm     .       .
```

## prefix.files/ 
- Is a directory containing intermediate files
