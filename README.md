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
```
Tool            Environment
AMRFinderPlus   amrfinder
Abricate
Prokka
mlst
```

## 2.(7) Prokka
Prokka will be used to predict genes in the input genome sequences.
To ensure that you can always use this version of prokka,
I added the conda package file of prokka version 1.14.6 in a zenodo repository.
The package file `prokka-1.14.6.tar.gz` was created by installing and packaging the Prokka according to its [source](https://github.com/tseemann/prokka) using following commands:
```
conda create -n prokka-1.14.6 prokka=1.14.6 -c conda-forge -c bioconda
conda pack -n prokka-1.14.6 -o prokka-1.14.6.tar.gz
```

What you should do is to download [prokka-1.14.6.tar.gz](https://zenodo.org/records/15104799/files/prokka-1.14.6.tar.gz), 
and unpack prokka-1.14.6.tar.gz in your GenoFaecium directory
```
mkdir /home/kihyunee/GenoFaecium/conda_packages
cd /home/kihyunee/GenoFaecium/conda_packages

mkdir /home/kihyunee/GenoFaecium/conda_packages/prokka-1.14.6
cd /home/kihyunee/GenoFaecium/conda_packages/prokka-1.14.6
wget https://zenodo.org/records/15104799/files/prokka-1.14.6.tar.gz --no-check-certificate
tar -xzf prokka-1.14.6.tar.gz
```

When the GenoFaecium activates prokka-1.14.6, it will be done by: 
```
source /home/kihyunee/GenoFaecium/conda_packages/prokka-1.14.6/bin/activate
```
When the GenoFaecium deactivates prokka-1.14.6, it will be done by:
```
source ~/.bashrc
```


## 2.(8) AMRFinderPlus
AMRFinderPlus will be used to find resistance genes in the input genome sequences.
To ensure that you can always use this version of amrfinder,
I added the conda package file of amrfinder version 3.12.8 in a zenodo repository.
The package file `amrfinder-3.12.8.tar.gz` was created by installing and packaing the AMRFinderPlus according to its [source](https://github.com/ncbi/amr/wiki) using following commands:
```
conda create -y -c conda-forge -c bioconda -n amrfinder-3.12.8 --strict-channel-priority ncbi-amrfinderplus=3.12.8
conda pack -n amrfinder-3.12.8 -o amrfinder-3.12.8.tar.gz
```

What you should do is to download [amrfinder-3.12.8.tar.gz](https://zenodo.org/records/15104799/files/amrfinder-3.12.8.tar.gz), 
and unpack amrfinder-3.12.8.tar.gz in your GenoFaecium directory
```
mkdir /home/kihyunee/GenoFaecium/conda_packages/amrfinder-3.12.8
cd /home/kihyunee/GenoFaecium/conda_packages/amrfinder-3.12.8
wget https://zenodo.org/records/15104799/files/amrfinder-3.12.8.tar.gz --no-check-certificate
tar -xzf amrfinder-3.12.8.tar.gz
```

When the GenoFaecium activates amrfinder-3.12.8, it will be done by: 
```
source /home/kihyunee/GenoFaecium/conda_packages/amrfinder-3.12.8/bin/activate
```

When the GeGenoFaecium deactivates amrfinder-3.12.8, it will be done by:
```
source ~/.bashrc
```


## 2.(9) Abricate
Abricate will be used to find virulence factor genes in the input genome sequences.
To ensure that you can always use this version of abricate,
I added the conda package file of abricate version 1.0.1 in a zenodo repository.
The package file `abricate-1.0.1.tar.gz` was created by installing and packaing the Abricate according to its [source](https://github.com/tseemann/abricate) using following commands:
```
conda create -n abricate-1.0.1 -c conda-forge -c bioconda -c defaults abricate=1.0.1
conda pack -n abricate-1.0.1 -o abricate-1.0.1.tar.gz
```

What you should do is to download [abricate-1.0.1.tar.gz](https://zenodo.org/records/15104799/files/abricate-1.0.1.tar.gz?download=1), 
and unpack abricate-1.0.1.tar.gz in your GenoFaecium directory
```
mkdir /home/kihyunee/GenoFaecium/conda_packages/abricate-1.0.1
cd /home/kihyunee/GenoFaecium/conda_packages/abricate-1.0.1
wget https://zenodo.org/records/15104799/files/abricate-1.0.1.tar.gz --no-check-certificate
tar -xzf abricate-1.0.1.tar.gz
```

When the GenoFaecium activates abricate-1.0.1, it will be done by: 
```
source /home/kihyunee/GenoFaecium/conda_packages/abricate-1.0.1/bin/activate
```

When the GeGenoFaecium deactivates abricate-1.0.1, it will be done by:
```
source ~/.bashrc
```


## 2.(10) Mlst
Seemann's mlst will be used to assign MLST to the input genome sequences.
To ensure that you can always use this version of abricate,
I added the conda package file of mlst version 2.23.0 in a zenodo repository.
The package file `mlst-2.23.0.tar.gz` was created by installing and packaing the mlst according to its [source](https://github.com/tseemann/mlst) using following commands:
```
conda create -n mlst-2.23.0 mlst=2.23.0 -c conda-forge -c bioconda -c defaults
conda pack -n mlst-2.23.0 -o mlst-2.23.0.tar.gz
```

What you should do is to download [mlst-2.23.0.tar.gz](https://zenodo.org/records/15104799/files/mlst-2.23.0.tar.gz?download=1), 
and unpack mlst-2.23.0.tar.gz in your GenoFaecium directory
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




