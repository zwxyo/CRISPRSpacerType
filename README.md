## **CRISPRSpacerType**
***
This is a tool for typing _Cronobacter_ based on CRISPR spacer.  
If you want to use this tool, please cite:  
- Mitrofanov A, Alkhnbashi O S, Shmakov S A, et al. CRISPRidentify: identification of CRISPR arrays using machine learning approach[J]. Nucleic acids research, 2021, 49(4): e20-e20.
- Zhang T, Jia Y, Li H, et al. CRISPRCasStack: A stacking strategy-based ensemble learning framework for accurate identification of Cas proteins[J]. Briefings in bioinformatics, 2022, 23(5): bbac335.
- Abby S S, Néron B, Ménager H, et al. MacSyFinder: a program to mine genomes for molecular systems with an application to CRISPR-Cas systems[J]. PloS one, 2014, 9(10): e110726.
- Néron B, Denise R, Coluzzi C, et al. MacSyFinder v2: Improved modelling and search engine to identify molecular systems in genomes[J]. Peer Community Journal, 2023, 3.

## Installation
***
Please make sure you have [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installed before proceeding.

```
conda create -n crisprspacertype
```
```
conda activate crisprspacertype
```
```
wget https://github.com/zwxyo/CRISPRSpacerType/releases/download/v0.1/crisprspacertype_channel.tar.gz
```
```
tar -xvzf crisprspacertype_channel.tar.gz
```
```
conda install -c ./crisprspacertype_channel crisprspacertype
```


## Usage

### Activate environment
Make sure that the dependencies in the CRISPRSpacerType.yaml file are all installed !

```
conda activate crisprspacertype
```

### See more options
You can see more parameters using `-h` or `--help`.

```
CRISPRSpacerType -h
CRISPRSpacerType --help
```



### Input: Whole Genome Sequences (WGS)

- WGS inputs must be in **FASTA (`.fna`) format**. 
- Each file should start with a header line beginning with "`>`", followed by the sequence identifier (e.g., `accession number`), strain name, and optional description. 
- The next lines should contain the nucleotide sequence.  

**Example format:** 
> \>NC_009778.1 Cronobacter sakazakii ATCC BAA-894, complete sequence
> ACGCTTCCCGGAAGACACGCCAAACGGTGTTTATGCAACAATGCGCGTAGTTTCCATAAAAGAAATGAAGAGAATTGGCG...

**CRISPR typing:** 

CRISPR classfication associated genes' BLAST use local database.
```
CRISPRSpacerType --ct --cas --bt s
```

CRISPR classfication associated genes' BLAST use custom database building by `makeblastdb`.

The `--db` parameter specifies the path to (`.fna`) files that will be used to build the BLAST database.  
```
CRISPRSpacerType --ct --cas --bt m --db <path>
```

**CRISPR identification:** Identifying CRISPR with [CRISPRidentify](https://github.com/BackofenLab/CRISPRidentify) and [CRISPRCasStack](https://github.com/yrjia1015/CRISPRCasStack)

```
CRISPRSpacerType -i input_file.fna -o output_folder
```

**Cas identificaton:** Identifying Cas with [MacSyFidner](https://github.com/gem-pasteur/macsyfinder) and [Prodigal](https://github.com/hyattpd/Prodigal)
```
CRISPRSpacerType --cas
```

**ST:** Multi-Locus Sequence Typing using [MLST](https://github.com/tseemann/mlst)
```
CRISPRSpacerType --mlst
```

### Input: PCR Sequences

- For each strain, integrate all PCR sequencing results into a single `.fna` file.  
- The first line must start with "`>`", followed by the CRISPR type and ID (e.g., `CRISPR1_210407132`), strain name.
- The next lines should contain the nucleotide sequence.  

**Example format:** 
> \>CRISPR1_210407132 Cronobacter sakazakii strain 7G
> CCTGCGCCTTGTCTCTTTTTTACCGCTTGAAAATCAGTAAGTTAATCGTTCTTTAACAATGTGAGATTGTGAACCAAATG...

**CRISPR typing:**
```
CRISPRSpacerType --pcr --ct
```

## 