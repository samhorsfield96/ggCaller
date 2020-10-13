# ggCaller: a gene caller for Bifrost graphs

Traverses Bifrost graphs constructed from bacterial genomes to identify putative protein coding sequences, known as open reading frames (ORFs).

ggCaller uses pyGFA to convert Bifrost graphs to Networkx graph objects to enable graph traversal.

##Installation

### Installation via conda (recommended)

Install through [bioconda](http://bioconda.github.io/):

```conda install ggcaller```

If conda is not installed, first install [miniconda](https://docs.conda.io/en/latest/miniconda.html), then add the correct channels:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Installation from source

Requirements:
-   ```python3```
-   ```biopython```
-   ```networkx```
-   ```SeQan3```
-   ```cmake```
-   ```openMP```
-   ```pthreads```
-   a C++17 compiler (e.g. gcc >=7.3)

```
git clone --recursive https://github.com/samhorsfield96/ggCaller 
python setup.py install
```

## Input files for ggCaller

ggCaller takes a list of fasta files, a Bifrost GFA file generated by ```Bifrost build```, and a TSV file generated by ```Bifrost query```. See the Bifrost repository for installation (https://github.com/pmelsted/bifrost).

### Generate Bifrost GFA

Using ```Bifrost build```, generate a GFA file with an associated colours file using a list of reference genomes. Here k-mer size is set to 31 bp (```-k 31```).

```Bifrost build -k 31 -c -r references.txt -o graph```

### Generate TSV file

The TSV file should be generated by querying the constituent unitigs of the GFA file against the GFA itself.
The accompanying script 'gfa_to_fasta.py' can be used to convert a GFA file to a FASTA file for use in ```Bifrost query```.

```python gfa_to_fasta.py graph.gfa```

The resulting FASTA file will take the same name as the graph, and can then be used in ```Bifrost query```. NOTE: exact k-mer matching should be used (```-e 1.0```).

```Bifrost query -e 1.0 -g graph.gfa -f graph.bfg_colors -q graph.fasta -o colours```

This will produce a TSV file named ```colours.tsv``` which can be used in ggCaller.

## Usage

To run ggCaller from the command line, pass arguments specifying input/output and parameters for traversal.
Arguments taken by ggcaller are:
- ```--graph``` Input GFA
- ```--colours``` Input TSV
- ```--source``` List of fasta files used to generate Input GFA (one file path per line)
- ```--kmer``` k-mer size used in ```Bifrost build```
- ```--path``` maximum path length in base-pairs (default: 10000 bp)
- ```--orf``` minimum ORF length to return in base-pairs (default: 90 bp)
- ```--write-idx``` Write FMIndexes to file, in same directory as fasta files (default: true)
- ```--threads``` Number of threads for FMIndexing (default: 1)
- ```--out``` output file in FASTA format (default: 'calls.fasta')

In the below example, ggCaller is run against ```graph.gfa```, ```colours.tsv``` and ```references.txt```, querying a graph constructed using k=31. Maximum path length is set to 5000 bp, and minimum ORF length is set to 150 bp. A FASTA file ```calls.fasta``` is then generated in the same directory.

```ggcaller --graph graph.gfa --colours colours.tsv --source references.txt --kmer 31 --path 5000 --orf 150 --out calls.fasta```

Test data is available in the ```data``` directory.

## Interpreting output FASTA

The output ```calls.fasta``` contains ORF sequences, with headers which contain a unique identifier (```Gene_ID```), the strand of the ORF (```Strand```) and a colours array (```Colours```), which describes the presence/absence of the ORF in the source genomes in the same order as the matrix in ```colours.tsv```.

```
>[Gene_ID: 1] [Strand: +] [Colours: ['0', '0', '1', '1']]
TTGTCATTTTATTGGTTTATTTTCTTAGCTTTGTTAGAGAGACAGAACTTGAACGTTCTTCTAG
>[Gene_ID: 2] [Strand: +] [Colours: ['1', '1', '1', '1']]
ATGGGATTAATGCTGTCTTTATGGGTGTTGGTGGCAGTTTTGATGTATTATCAGGACACATTAAACGAGCTCCATTATGGATGCAAAAATTGA
>[Gene_ID: 3] [Strand: +] [Colours: ['0', '0', '0', '1']]
TTGTTCAAAGGTGGTGTTACGATTTCAAGAACTCCTCTCAGTTCTGAGGACACGGTAATGATTGATGCGATAG
```

## Citation

If you use this code, please cite:

Bifrost: 
Holley G., Melsted, P. Bifrost – Highly parallel construction and indexing of colored and compacted de Bruijn graphs. (2019) bioRxiv 695338. doi: https://doi.org/10.1101/695338

pyGFA: 
https://github.com/AlgoLab/pygfa

Networkx: 
Hagberg AA, Schult DA, Swart PJ. Exploring Network Structure, Dynamics, and Function using NetworkX. In: Varoquaux G, Vaught T, Millman J, editors. Proceedings of the 7th Python in Science Conference. Pasadena, CA USA; 2008. p. 11–5

SeQan3: 
Reinert, K. et al. The SeqAn C++ template library for efficient sequence analysis: A resource for programmers. (2017) Journal of biotechnology, 261, 157-168 doi: https://doi.org/10.1016/j.jbiotec.2017.07.017




