## **CRISPRSpacerType**
***
This is a tool for typing Cronobacter based on CRISPR spacer.  
If you want to use this tool, please cite:  
- Mitrofanov A, Alkhnbashi O S, Shmakov S A, et al. CRISPRidentify: identification of CRISPR arrays using machine learning approach[J]. Nucleic acids research, 2021, 49(4): e20-e20.
- Zhang T, Jia Y, Li H, et al. CRISPRCasStack: A stacking strategy-based ensemble learning framework for accurate identification of Cas proteins[J]. Briefings in bioinformatics, 2022, 23(5): bbac335.
- Abby S S, Néron B, Ménager H, et al. MacSyFinder: a program to mine genomes for molecular systems with an application to CRISPR-Cas systems[J]. PloS one, 2014, 9(10): e110726.
- Néron B, Denise R, Coluzzi C, et al. MacSyFinder v2: Improved modelling and search engine to identify molecular systems in genomes[J]. Peer Community Journal, 2023, 3.

## Installation
***
Configure the **conda environment** before installing to the system.

`conda env create -f CRISPRSpacerType.yml`
- python 3.7.10
- blast 2.12.0
- prodigal 2.6.3
- mlst 2.23.0
- macsyfinder 2.1.3

### CentOS
` `

### Ubuntu
` `

## How to use

### See more options
`CRISPRSpacerType -h`
or
`CRISPRSpacerType --help`

### Activate environemt
Make sure that the dependencies in the CRISPRSpacerType.yml file are all installed !

`conda activate CRISPRSpacerType`

### CRISPR typing number

`CRISPRSpacerType --ct`

`CRISPRSpacerType --ct --cas`  

`CRISPRSpacerType --ct --cas --bt s`

`CRISPRSpacerType --ct --bt m --bi  --db`

### CRISPR identification
Identifying CRISPR with CRISPRCasStack

`CRISPRSpacerType --csi input_file.fna --cso output_folder`

### Cas identificaton
Identifying Cas with MacSyFidner and prodigal

`CRISPRSpacerType --cas`

Due to GitHub's file size constraints, if you 
