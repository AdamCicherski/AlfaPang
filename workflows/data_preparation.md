### E. coli 50, 100, 200 and 400
E. coli data used in these experiments was downloaded from [Zenodo](https://zenodo.org/)

For 50 ecoli sequences run:
```bash
wget https://zenodo.org/records/7937947/files/ecoli50.fa.gz
gzip -d ecoli50.fa.gz
```
For the remaining datasets, run:
```bash
wget https://zenodo.org/records/7937947/files/ecoli500.fa.gz
gzip -d ecoli500.fa.gz
python AlfaPang/scripts/preprocess_data.py ecoli500.fa AlfaPang/data/ecoli100_names.txt ecoli100.fa 
python AlfaPang/scripts/preprocess_data.py ecoli500.fa AlfaPang/data/ecoli200_names.txt ecoli200.fa
python AlfaPang/scripts/preprocess_data.py ecoli500.fa AlfaPang/data/ecoli400_names.txt ecoli400.fa
```
### E. coli 800, 1600 and 3412

E. coli data used in these experiments was downloaded from [GenBank](`https://www.ncbi.nlm.nih.gov/assembly`).

- **Date:** 18.06.2024
- **Release Type:** RefSeq
- **Query title:** `Search txid562[Organism] AND (latest[filter] AND "complete genome"[filter] AND all[filter] NOT partial[filter]) AND latest[filter] AND all[filter] NOT anomalous[filter]`
- **Search results count:** 4534
- **Filtered out:** 1122 entries that do not have the requested Release Type or are suppressed.
- **Remaining entries:** 3412 assemblies
- **Further filtering:** Removed plasmid sequences and sampled the first 800 and 1600 sequences.

Since June 2024, Genome and Assembly have been replaced with the NCBI datasets. We recommend downloading data using the datasets CLI tool:
```bash
mkdir ecoli3412_archives
while IFS= read -r line; do
datasets download genome accession "$line" --filename "ecoli3412_archives/${line}.zip"
done < AlfaPang/data/ecoli_3412_accession.txt
```

