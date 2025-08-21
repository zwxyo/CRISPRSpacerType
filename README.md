## **CRISPRSpacerType**
***
This is a tool for typing _Cronobacter_ based on CRISPR spacer.  
If you want to use this tool, please cite:  
- Mitrofanov A, Alkhnbashi O S, Shmakov S A, et al. CRISPRidentify: identification of CRISPR arrays using machine learning approach[J]. Nucleic acids research, 2021, 49(4): e20-e20.
- Zhang T, Jia Y, Li H, et al. CRISPRCasStack: A stacking strategy-based ensemble learning framework for accurate identification of Cas proteins[J]. Briefings in bioinformatics, 2022, 23(5): bbac335.
- Abby S S, Néron B, Ménager H, et al. MacSyFinder: a program to mine genomes for molecular systems with an application to CRISPR-Cas systems[J]. PloS one, 2014, 9(10): e110726.
- Néron B, Denise R, Coluzzi C, et al. MacSyFinder v2: Improved modelling and search engine to identify molecular systems in genomes[J]. Peer Community Journal, 2023, 3.

## Installation via Conda
***
Please make sure you have [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) installed before proceeding.

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


## How to use

### Activate environemt
Make sure that the dependencies in the CRISPRSpacerType.yaml file are all installed !

```
conda activate crisprspacertype
```

### See more options

```
CRISPRSpacerType -h
```

or

```
CRISPRSpacerType --help
```



### Input: Whole Genome Sequences (WGS)
- 


CT: CRISPR typing number

```
CRISPRSpacerType --ct --cas --bt s
```

```
CRISPRSpacerType --ct --cas --bt m --bi <path> --db <path>
```

CRISPR identification: Identifying CRISPR with CRISPRCasStack

```
CRISPRSpacerType -i input_file.fna -o output_folder
```

Cas identificaton: Identifying Cas with MacSyFidner and prodigal
```
CRISPRSpacerType --cas
```

ST: Multi-Locus Sequence Typing using MLST
```
CRISPRSpacerType --mlst
```

### Input: PCR Sequences
CT: CRISPR typing number
```
CRISPRSpacerType --pcr --ct
```