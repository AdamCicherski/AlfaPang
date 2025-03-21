# AlfaPang
AlfaPang constructs variation graphs, leveraging its alignment-free and reference-free approach, based solely on intrinsic sequence properties. This design allows AlfaPang's runtime and memory usage to scale linearly with the size of input sequences, enabling it to handle significantly larger genome sets compared to other methods. 

## Instalation
To install AlfaPang, download the repo, then navigate to the project root and run following commands:

```bash
    git submodule update --init --recursive  
    mkdir build  
    cd build  
    cmake ..  
    cmake --build .
```
## Usage

To run the program, use the following command format:

```bash
./AlfaPang <input_fasta> <output_gfa> <k>
```
**Note**: The value of k must be an odd integer.


## Example
Download example data with following command:

```bash
wget https://zenodo.org/records/7937947/files/ecoli50.fa.gz
gzip -d ecoli50.fa.gz 
```
Then, try::
```bash
./AlfaPang ecoli50.fa ecoli50_ap.gfa 47
```
Refining AlfaPang graphs is highly recommended. We suggest using [smoothxg](https://github.com/pangenome/smoothxg) and [gfaffix](https://github.com/marschall-lab/GFAffix).
Both tools with dependencies can be easily installed following [pggb](https://github.com/pangenome/pggb) repository.
Basic commands for graph refainment: 

```bash
smoothxg -g ecoli50_ap.gfa -r 50 -V -o ecoli50_ap_smooth.gfa
gfaffix ecoli50_ap_smoth.gfa -o ecoli50_final.gfa
```

## Parameter k choice
We suggest choosing parameter $k$ based on the fraction of rare $k$-mers (those occurring only once in the $k$-mer spectrum). In our tests, values of $k$ yielding around 5% rare $k$-mers result in a reasonable graph structure.  

For this purpose, we provide the script `kmer_fractions.sh`, which uses the disk-based $k$-mer counter [KMC](https://github.com/refresh-bio/KMC). The script produces a `.tsv` file with the fraction of rare $k$-mers calculated for a given $k$ range.
```bash
./AlfaPang/scripts/kmer_fractions.sh -i <input_fasta_file> -o <output_dir_name> -k <min_k_value> -K <max_k_value> -s <step> 
```
**Note**: If KMC executable path in your system is different than `${HOME}/kmc/bin/kmc` modify script variable `KMC_PATH`. 


