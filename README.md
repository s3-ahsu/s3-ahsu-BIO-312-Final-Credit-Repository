# s3-ahsu-BIO-312-Final Credit Repository
This repository contains all the commands that would need to be run to re-do the analysis for labs 3,4,5,6,and 8. 

# Contents

1. Lab 3: Finding Homologs with BLAST  
2. Lab 4: Gene Family Sequence Alignment 
3. Lab 5: Gene Family Phylogeny using IQ-TREE 
4. Lab 6: Gene-and-Species Tree Reconciliation
5. Lab 8: Protein Domain Prediction

# Finding Homologs with BLAST 
This lab aimed to identify homologs of my gene family against known BLAST databases. Only high-scoring matches were needed, so the e-value was required to be less than 1e-30.

To start the lab, create a new folder for the gene family:
```
mkdir ~/lab03-$MYGIT/NP_036387.2 
```
Once you created the folder, use the cd command to make sure it is there:
```
cd ~/lab03-$MYGIT/NP_036387.2
```
To make sure you are in the current folder, use the pwd command. 
```
pwd
```
To download the globin protein from Homo sapiens as query sequence:
```
ncbi-acc-download -F fasta -m protein "NP_036387.2" 
```
Next, we need to perform a BLAST search using the query protein:
```
blastp -db ../allprotein.fas -query NP_036387.2.fa -outfmt 0 -max_hsps 1 -out globins.blastp.typical.out
```
To view the BLAST output file, use the less command.
```
less globins.blastp.typical.out
```
We can create the same analysis but with a more detailed and easier-to-process output by requesting for tabuluar output:
```
blastp -db ../allprotein.fas -query NP_036387.2.fa  -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out globins.blastp.detail.out 
```
We can take a look at the output by using the less command again.
```
less -S globins.blastp.detail.out
```
Instead of counting by hand, we can compute the total human hits by using the grep command:
```
grep -c H.sapiens globins.blastp.detail.out
```







