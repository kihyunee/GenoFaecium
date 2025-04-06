# 1. GenoFaecium
Genomic analysis framework for Enterococcus faecium clinical isolates 


# 2. Installation
## 2.(1) Linux
The installation guide is based on Linux OS.
I tested this installation process in Windows Ubuntu 24.04.1 LTS.


## 2.(2) Miniconda
If you don't have conda in your system, install a miniconda according to [Miniconda Quickstart Installation Instruction](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions).
Here are the commands you can use (copied from the aforementioned link).

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


## 2.(3) Some of the basic dependencies you may miss in fresh ubuntu & required in FastANI or Minimap
```
sudo apt update
sudo apt install unzip
sudo apt install libgomp1
```


## Next steps
- Clone the GenoFaecium directory onto your system
- Install FastANI
- Install Prokka
- Install AMRFinderPlus
- Install Abricate
- Install mlst
- Install Minimap




## 2.(4) Clone the GenoFaecium base directory
Decide the location to create the genofaecium directory.
I will use '~/' ( = '/home/kihyunee/').

```
# Go to the location you decided
cd /home/kihyunee/

# Clone the repository
git clone https://github.com/kihyunee/GenoFaecium.git

# Enter the GenoFaecium directory
cd /home/kihyunee/GenoFaecium
mkdir dependency_binary
```


## 2.(5) FastANI version 1.34
We use [fastANI](https://github.com/ParBLiSS/FastANI) in the identity confirmation step.
We will simply download the executable binary file that the developers have distributed and put it in the GenoFaecium/install directory.

```
# Go into the directory where you want to have FastANI and Minimap.
cd /home/kihyunee/GenoFaecium/dependency_binary

# Download the zipped binary file from the FastANI developer's site
wget https://github.com/ParBLiSS/FastANI/releases/download/v1.34/fastANI-linux64-v1.34.zip

# Extract fastANI
unzip fastANI-linux64-v1.34.zip

# You will have a file called 'fastANI' now. Add this executable 'fastANI' file to your path
# add the current directory (~/GenoFaecium/dependency_binary) to your PATH
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc
source ~/.bashrc
```


## 2.(6) Minimap2 version 2.28
We use [Minimap2](https://github.com/lh3/minimap2) as the core engine of coreSNP pipeline.
You can download executable file like this, according to [Minimap2 readme](https://github.com/lh3/minimap2?tab=readme-ov-file#install)
```
cd /home/kihyunee/GenoFaecium/dependency_binary
curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf -
cp minimap2-2.28_x64-linux/minimap2 .
```

## Rest of the installations will use conda
From here on, we install the tools using conda, in conda environment made for each specific tool.
However, **you don't have to activate any of these environments before executing the Genofaecium script** because the script takes care of environments switching.

## 2.(7) Prokka
Prokka will be used to predict genes in the input genome sequences.
To ensure that you can always use this version of prokka,
I added the conda package file of prokka version 1.14.6 in a zenodo repository.
The package file `prokka-1.14.6.tar.gz` was created by installing and packaging the Prokka according to its [source](https://github.com/tseemann/prokka) using following commands:
**You dont't have to do this because this just tells you how the package file was prepared.**
```
conda create -n prokka-1.14.6 prokka=1.14.6 -c conda-forge -c bioconda
conda pack -n prokka-1.14.6 -o prokka-1.14.6.tar.gz
```

**What you should do is to download [prokka-1.14.6.tar.gz](https://zenodo.org/records/15104799/files/prokka-1.14.6.tar.gz), 
and unpack prokka-1.14.6.tar.gz in your GenoFaecium directory:**
```
mkdir /home/kihyunee/GenoFaecium/conda_packages
cd /home/kihyunee/GenoFaecium/conda_packages

mkdir /home/kihyunee/GenoFaecium/conda_packages/prokka-1.14.6
cd /home/kihyunee/GenoFaecium/conda_packages/prokka-1.14.6
wget https://zenodo.org/records/15104799/files/prokka-1.14.6.tar.gz --no-check-certificate
tar -xzf prokka-1.14.6.tar.gz
```


## 2.(8) AMRFinderPlus
AMRFinderPlus will be used to find resistance genes in the input genome sequences.
To ensure that you can always use this version of amrfinder,
I added the conda package file of amrfinder version 3.12.8 in a zenodo repository.
The package file `amrfinder-3.12.8.tar.gz` was created by installing and packaing the AMRFinderPlus according to its [source](https://github.com/ncbi/amr/wiki) using following commands:
**You dont't have to do this because this just tells you how the package file was prepared.**
```
conda create -y -c conda-forge -c bioconda -n amrfinder-3.12.8 --strict-channel-priority ncbi-amrfinderplus=3.12.8
conda pack -n amrfinder-3.12.8 -o amrfinder-3.12.8.tar.gz
```

**What you should do is to download [amrfinder-3.12.8.tar.gz](https://zenodo.org/records/15104799/files/amrfinder-3.12.8.tar.gz), 
and unpack amrfinder-3.12.8.tar.gz in your GenoFaecium directory:**
```
mkdir /home/kihyunee/GenoFaecium/conda_packages/amrfinder-3.12.8
cd /home/kihyunee/GenoFaecium/conda_packages/amrfinder-3.12.8
wget https://zenodo.org/records/15104799/files/amrfinder-3.12.8.tar.gz --no-check-certificate
tar -xzf amrfinder-3.12.8.tar.gz
```


## 2.(9) Abricate
Abricate will be used to find virulence factor genes in the input genome sequences.
To ensure that you can always use this version of abricate,
I added the conda package file of abricate version 1.0.1 in a zenodo repository.
The package file `abricate-1.0.1.tar.gz` was created by installing and packaing the Abricate according to its [source](https://github.com/tseemann/abricate) using following commands:
**You dont't have to do this because this just tells you how the package file was prepared.**
```
conda create -n abricate-1.0.1 -c conda-forge -c bioconda -c defaults abricate=1.0.1
conda pack -n abricate-1.0.1 -o abricate-1.0.1.tar.gz
```

**What you should do is to download [abricate-1.0.1.tar.gz](https://zenodo.org/records/15104799/files/abricate-1.0.1.tar.gz?download=1), 
and unpack abricate-1.0.1.tar.gz in your GenoFaecium directory:**
```
mkdir /home/kihyunee/GenoFaecium/conda_packages/abricate-1.0.1
cd /home/kihyunee/GenoFaecium/conda_packages/abricate-1.0.1
wget https://zenodo.org/records/15104799/files/abricate-1.0.1.tar.gz --no-check-certificate
tar -xzf abricate-1.0.1.tar.gz
```



## 2.(10) Mlst
Seemann's mlst will be used to assign MLST to the input genome sequences.
To ensure that you can always use this version of abricate,
I added the conda package file of mlst version 2.23.0 in a zenodo repository.
The package file `mlst-2.23.0.tar.gz` was created by installing and packaing the mlst according to its [source](https://github.com/tseemann/mlst) using following commands:
**You dont't have to do this because this just tells you how the package file was prepared.**
```
conda create -n mlst-2.23.0 mlst=2.23.0 -c conda-forge -c bioconda -c defaults
conda pack -n mlst-2.23.0 -o mlst-2.23.0.tar.gz
```

**What you should do is to download [mlst-2.23.0.tar.gz](https://zenodo.org/records/15104799/files/mlst-2.23.0.tar.gz?download=1), 
and unpack mlst-2.23.0.tar.gz in your GenoFaecium directory:**
```
mkdir /home/kihyunee/GenoFaecium/conda_packages/mlst-2.23.0
cd /home/kihyunee/GenoFaecium/conda_packages/mlst-2.23.0
wget https://zenodo.org/records/15104799/files/mlst-2.23.0.tar.gz --no-check-certificate
tar -xzf mlst-2.23.0.tar.gz
```

When the GenoFaecium activates mlst-2.23.0, it will be done by: 
```
source /home/kihyunee/GenoFaecium/conda_packages/mlst-2.23.0/bin/activate
```

When the GeGenoFaecium deactivates mlst-2.23.0, it will be done by:
```
source ~/.bashrc
```


# 3. Pre-compiled database for GenoFaecium
Go to the GenoFaecium directory (the directory you cloned above).
Download the [genofaecium_db_pre_compiled.tar.gz](https://zenodo.org/records/15104799/files/genofaecium_db_pre_compiled.tar.gz) file containing pre-compiled databases and extract it.
```
cd /home/kihyunee/GenoFaecium
wget -O genofaecium_db_pre_compiled.tar.gz https://zenodo.org/records/15104799/files/genofaecium_db_pre_compiled.tar.gz --no-check-certificate
tar -xzf genofaecium_db_pre_compiled.tar.gz
```


# 4. Execute GenoFaecium

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

For example, to extract results from E. faecium genome fasta file at `test_input/GCA_`.
```
python genofaecium.py --fasta test_input/GCA_047785525.1.fasta --out test_output/GCA_047785525.1 --tool_dir /home/kihyunee/GenoFaecium
```

The output files 
## <prefix>.result_col.txt
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

## <preifx>.result_row.txt
- Is a single-row text file, tab-delimited list of report values ordered as the rows presented in the above column-format output file.
- This is made for you to easily concatenate the reports generated for many isolates into a single table using `cat`.
- Example output content:
```
GCA_047785525.1 296     3079105 Enterococcus faecium    80      aac(6')-I,aac(6')-Ie/aph(2'')-Ia,aph(3')-IIIa,ant(6)-Ia pbp5_M485A      .       .       vanA.       .       .       .       .       .       .       tet(S),tet(L)   liaR_W73C,liaS_T120A    .       .       .       .       gyrA_S83I,parC_S80I     msr(C),erm(B)       acm     sgrA    fss3    .       scm     .       .
```

## <prefix>.files/ 
- Is a directory containing intermediate files
