# BIO 312 Final Credit Repository-s3-ahsu
This repository contains all the commands that would need to be run to re-do the analysis for labs 3-6, 8.

# Contents

1. Lab 3: Finding Homologs with BLAST  
2. Lab 4: Gene Family Sequence Alignment 
3. Lab 5: Gene Family Phylogeny using IQ-TREE 
4. Lab 6: Reconciliaing a Gene and Species Tree
5. Lab 8: Protein Domain Prediction

# Finding Homologs with BLAST 
This lab aimed to identify homologs of my gene family against known BLAST databases. Only high-scoring matches were needed, so the e-value was required to be less than 1e-30.


To start the lab, create a new folder for the gene family:
```
mkdir ~/lab03-$MYGIT/NP_036387.2 
```
Once you created the folder, use the cd command to go there:
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
Next, we would like to filter the BLAST output file for high-scoring homologs only. To do this, we need to set the e-value less than 1e-30. 
```
awk '{if ($6< 1e-30)print $1 }' globins.blastp.detail.out > globins.blastp.detail.filtered.out
```
To count the total number of hits in the BLAST results after the filter easier, we can use the wc command. 
```
wc -l globins.blastp.detail.filtered.out
```
Lastly, we want to find how many homologs are in each species.
```
grep -o -E "^[A-Z]\.[a-z]+" globins.blastp.detail.filtered.out  | sort | uniq -c
```
It is desirable to work between 20 and 85 homologs. If it is not within that range, you may need to change the e-value threshold that will either decrease or increase the number of hits. 


# Gene Family Sequence Alignment 
In this lab, we aligned the sequences that we obtained from last lab and uncover more information about them, including percent identity, length of the alignment, and homologous positions. 


To start the lab, create a new folder for the gene family. 
```
mkdir ~/lab04-$MYGIT/NP_036387.2 
```
Once you created the folder, use the cd command to go there:
```
cd ~/lab04-$MYGIT/NP_036387.2 
```
To make sure you are in the current folder, use the pwd command. 
```
pwd
```
 We need to obtain the sequences that are in the BLAST output file from last lab before we can align them:
 ```
seqkit grep --pattern-file~/lab03-$MYGIT/NP_036387.2/NP_036387.2.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homologs.fas
```
Next, we can make a multiple sequence alignment using the muscle command.
```
muscle -align ~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homologs.fas -output ~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homologs.al.fas
```
To view the alignment that was made, use the alv command.
```
alv -kli  ~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homologs.al.fas | less -RS 
```
We can look at the majority of the alignment using the alv command again.
```
alv -kli --majority ~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homologs.al.fas | less -RS
```
We can also print the alignment to a large pdf file to see better.
```
Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R
~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homologs.al.fas
```
To calculate the width (length) of the alignment, alignbuddy can be used.
```
alignbuddy  -al  ~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homologs.al.fas 
```
If we want to calculate the length of the alignment after removing any column with gaps, here is how to do it.
```
alignbuddy -trm all  ~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homologs.al.fas | alignbuddy  -al 
```
If we want to calculate the length of the alignment after removing invariant (completely conserved) positions, here is how to do it.
```
alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homologs.al.fas | alignbuddy  -al
```
One way to calculate average percent identity is by using t_coffee. 
```
t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homologs.al.fas -output sim
```
Another way to calculate average percent identity is by using alignbuddy.
```
alignbuddy -pi ~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homologs.al.fas | awk ' (NR>2)
{ for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } '
```

# Gene Family Phylogeny using IQ-TREE
This lab demonstrated how to construct a phylogenetic tree for the homologs we found from sequence data in lab 4. 


To start the lab, create a new folder for the gene family. 
```
mkdir ~/lab05-$MYGIT/NP_036387.2
```
Once you created the folder, use the cd command to go there:
```
cd ~/lab05-$MYGIT/NP_036387.2
```
We want to remove any sequences that contains a duplicate label tag and put a copy in lab 5 directory. 
```
sed 's/ /_/g'  ~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homolo
gs.al.fas | seqkit grep -v -r -p "dupelabel" >  ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.fas
```
We can find the maximum tree estimate. This calculates the optimal amino acid substitution model and amino acid frequencies. A tree search can also be performed that estimates the tree length.
```
iqtree -s ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.a
l.fas -bb 1000 -nt 2
```
To display the iqtree file that includes the ASCII graphics version of the tree, we can use the nw_display program. The .treefile is newick formatted. 
```
nw_display ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.fas.treefile
```
We can also look at the unrooted tree with a graphical display. To make it easier to read the genes, we can make the size of the text labels smaller and set the label lengths to 15.
```
Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R  
~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homol
ogsf.al.fas.treefile ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.fas.treefile.pdf 0.4 15
```
Next, we can use gotree to re-root the unrooted tree to midpoint tree rooting. 
```
gotree reroot midpoint -i ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile
```

# Reconciliaing a Gene and Species Tree
This lab ulitized programs such as Notung to perform a gene tree-species tree reconciliation, in which duplication events and loss events will be estimated. In addition, we were able to determine orthologous and paralogous relationships. 

To start the lab, create a new folder for the gene family:
```
mkdir ~/lab06-$MYGIT/NP_036387.2
```
Once you created the folder, use the cd command to go there:
```
cd ~/lab06-$MYGIT/NP_036387.2
```
To make a copy of 














