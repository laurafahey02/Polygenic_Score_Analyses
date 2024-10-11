# remove lines beginning with #
sed '/^#/d' Homo_sapiens.GRCh37.87.gtf > grch37.87.gtf

# print first and third column "groups"
awk -F";" '{print $1, $3}' grch37.87.gtf > gtf_13.txt

# Subset to only contain rows for genes 
awk '$3 == "gene"' gtf_13.txt > gtf_genes13.txt

# Subset to only contain columns chr, start, stop ...
awk '{print $1, $4, $5, $12}' gtf_genes13.txt | sed 's/"//g' > gtf_genes_corord.txt
 
# convert tasb to column in geneset
tr "\t" "\n" < ../pathway_files/REACTOME_MRNA_SPLICING_MINOR_PATHWAY.txt > rna_splice_minor.txt
# remove top two lines
sed -i "1d" rna_splice_minor.txt
sed -i "1d" rna_splice_minor.txt
