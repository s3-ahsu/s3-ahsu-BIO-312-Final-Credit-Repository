# BIO 312 Final Credit Repository-s3-ahsu
This repository contains all the commands that would need to be run to re-do the analysis for labs 3-6, 8.

# Contents

1. Lab 3: Finding Homologs with BLAST  
2. Lab 4: Gene Family Sequence Alignment 
3. Lab 5: Gene Family Phylogeny using IQ-TREE 
4. Lab 6: Reconciling a Gene and Species Tree
5. Lab 8: Protein Domain Prediction

# Finding Homologs with BLAST 
This lab aimed to identify the homologs of my gene family against known BLAST databases. Only high-scoring matches were needed, so the e-value was required to be less than 1e-30.


To start the lab, create a new folder for the gene family:
```
mkdir ~/lab03-$MYGIT/NP_036387.2 
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
We can create the same analysis but with a more detailed and easier-to-process output by requesting for tabular output:
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
To count the total number of hits in the BLAST results after the filter, we can use the wc command. 
```
wc -l globins.blastp.detail.filtered.out
```
We can also count the total number of hits before the filter.
```
wc -l globins.blastp.detail.out
```
Lastly, we want to find how many homologs are in each species.
```
grep -o -E "^[A-Z]\.[a-z]+" globins.blastp.detail.filtered.out  | sort | uniq -c
```
| Species      |   Count   |
|--------------|-----------|
| C.carcharias |     2      |
| C.mydas      |     2     |
| D.rerio      |     2      |
| E.caballus   |      2     |
| F.catus      |      2     |
| G.aculeatus  |      2     |
| G.gallus     |      2     |
| H.sapiens    |      2     |
| S.salar      |      4    |
| S.townsendi  |      2     |
| X.laevis     |      4     |

It is desirable to work between 20 and 85 homologs. If it is not within that range, you may need to change the e-value threshold that will either decrease or increase the number of hits. 
Since my gene family had a total of 26 homologs, the threshold was not modified. 


# Gene Family Sequence Alignment 
In this lab, we aligned the sequences that we obtained from last lab and discovered more information about them, including percent identity, length of the alignment, and homologous positions. 


To start the lab, create a new folder for the gene family. 
```
mkdir ~/lab04-$MYGIT/NP_036387.2 
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
We want to remove any sequences that contains a duplicate label tag and put a copy in the lab 5 directory. 
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
Next, we can use gotree to re-root the unrooted tree using midpoint tree rooting. 
```
gotree reroot midpoint -i ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile
```
We can now take a look at the new rooted tree at the command line.
```
nw_order -c n ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile  | nw_display -
```
However, for better visual purposes, we should make a graphic version of the tree. 
```
nw_order -c n ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s  >  ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile.svg -
```
Even better, we can also convert the svg image to a pdf.  
```
convert  ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homol
ogsf.al.mid.treefile.svg  ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile.pdf
```
Most programs, including the nw_display that we used, would display the tree as a phylogram. 
```
nw_order -c n ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile | nw_t
opology - | nw_display -s  -w 1000 > ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.midCl.treefile.svg -
```
However, this isn't as convenient in cases where the branch lengths are very short, which makes it difficult to view the individual clades. Therefore, we can convert the phylogram to a cladogram for better visualization.
```
convert ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.midCl.treefile.svg ~/lab0
5-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.midCl.treefile.pdf
```

# Reconciling a Gene and Species Tree
This lab ulitized programs such as Notung to perform a gene tree-species tree reconciliation, in which duplication events and loss events will be estimated. In addition, we were able to determine orthologous and paralogous relationships. 

To start the lab, create a new folder for the gene family:
```
mkdir ~/lab06-$MYGIT/NP_036387.2
```
Make a copy of the gene tree from the lab 5 gene family folder into lab 6. 
```
cp ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile 
~/lab06-$MYGIT/NP_036387.2/
```
To check if the gene tree is in the lab 6 folder, use ls.
```
ls ~/lab06-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile 
/home/bio312-user/lab06-s3-ahsu/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile
```
Next, we need to reconcile the gene tree and species tree in notung. 
```
java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/NP_036387.2/
```
We can also use nw_display to provide names for internal lineages.
```
nw_display ~/lab05-$MYGIT/species.tre
```
As an alternative, we can see how notung assigned node names to the internal nodes.
```
grep NOTUNG-SPECIES-TREE ~/lab06-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile.rec.ntg | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -
```
Next, we can generate a RECPhyloXML object and view the reconciliation via thirdkind.
```
 python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile.rec.ntg --include.species
```
We can then create a graphic version of the reconciliation by using thirdkind.
```
thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile.rec.ntg.xml -o  ~/lab06-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile.rec.svg
```
Lastly, for better viewing, we can convert the graphic to a pdf. 
```
convert  -density 150 ~/lab06-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile.rec.svg ~/lab06-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.mid.treefile.rec.pdf
```

# Protein Domain Prediction
This lab focused on the identification of Pfam domains within the protein sequences by the using RPS-BLAST.

First, make a new directory for the gene family sequences and change into that directory. 
```
mkdir ~/lab08-$MYGIT/NP_036387.2 && cd ~/lab08-$MYGIT/NP_036387.2
```
Then, make a copy of the raw unaligned sequence that removes any stop codons (which are represented by asterisks). We will use the sed command to substitute these asterisks with nothing and direct the output to our gene family folder in lab 8. 
```
sed 's/*//' ~/lab04-$MYGIT/NP_036387.2/NP_036387.2.homologs.fas > ~/lab08-$MYGIT/NP_036387.2/NP_036387.2.homologs.fas
```
We will now run rps blast by the commmand below:
```
rpsblast -query ~/lab08-$MYGIT/NP_036387.2/NP_036387.2.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/NP_036387.2/NP_036387.2.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```
To copy the final gene tree from Lab 5 to Lab 8, we use the cp command.
```
cp ~/lab05-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.fas.treefile ~/lab08-$MYGIT/NP_036387.2
```
We will run R script that plots the pfam domain predictions from rps-blast next to their protein on the phylogeny. Here's an easy way to do it from the command line without opening the console.
```
Rscript --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/NP_036387.2/NP_036387.2.homologsf.al.fas.treefile ~/lab08-$MYGIT/NP_036387.2/NP_036387.2.rps-blast.out ~/lab08-$MYGIT/NP_036387.2/NP_036387.2.tree.rps.pdf
```
You should see the predicted domains on the tips of the tree and a legend with all the names of the domains. If some proteins are missing domains, then you may need to fix it manually. 

We can use the mlr command to make better sense of the domains.
```
mlr --inidx --ifs "\t" --opprint cat ~/lab08-$MYGIT/NP_036387.2/NP_036387.2.rps-blast.out | tail -n +2 | less -S
```
Instead of counting by hand, we can see how many proteins have no annotations using the command below:
```
cut -f 1 ~/lab08-$MYGIT/NP_036387.2/NP_036387.2.rps-blast.out | sort | uniq -c
```
We can see which Pfam domain annotation is most commonly found by using the command below:
```
cut -f 6 ~/lab08-$MYGIT/NP_036387.2/NP_036387.2.rps-blast.out | sort | uniq -c
```
Here's an easier way to see which protein has the longest annotated protein domain:
```
awk '{a=$4-$3;print $1,'\t',a;}' ~/lab08-$MYGIT/NP_036387.2/NP_036387.2.rps-blast.out | sort -k2nr
```
We can also find the protein with a domain annotation that has the best e-value, which is the lowest one. 
```
cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/NP_036387.2/NP_036387.2.rps-blast.out
```

















