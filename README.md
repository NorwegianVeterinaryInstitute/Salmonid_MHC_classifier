# Salmonid_MHC_classifier
IPD-MHC classifier for Salmonids


## Dependencies
Muscle, emboss (transep and water tools) which can installed with Conda as below
```bash
$ conda install muscle emboss
```

Uses os, sys, subprocess, Bio and ete3 packages in python.
First three packages are available within python while the last two can be installed with Conda as below
```bash
$ conda install -c conda-forge biopython ete3 
```

## Initiation
```bash
$ python tool.py [input.fa] [Db_name] [report.txt]
```  

input.fa: Sequences in fasta format  
Db_name: SASA-DAA, SASA-DAB, SASA-UAB  
report.txt: Name for the output report file  

## Execution

## Output
