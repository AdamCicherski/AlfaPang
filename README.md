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





