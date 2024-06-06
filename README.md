# AlfaPang
AlfaPang constructs variation graphs, leveraging its alignment-free and reference-free approach, based solely on intrinsic sequence properties. This design allows AlfaPang's runtime and memory usage to scale linearly with the size of input sequences, enabling it to handle significantly larger genome sets compared to other methods. 

## Instalation
To install AlfaPang, download the repo, then go to the project root and run following commands:

```bash  
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
