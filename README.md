# Cell-cell-crosstalk-prediction
Using pseudo-bulk RNA-seq data to predict ligand-receptor relationships between cell populations


### How to run the program
First step is to compile the source code to executive file by C++ compliler 
```sh
g++ LR_generator.cpp -o LR_generator
```
Before running the program, you have to prepare three files, two tables of expressed gene list for cell type 1 and 2 and one table of human ligand-receptor list.
Usage: ./LR_generator No_of_CT1_clusters No_of_CT2_clusters CT1_table_file CT2_table_file L-R_pair_table

```sh
./LR_generator 10 2 macrophage_UQ_10CMs.tsv neutrophil_UQ_2CMs.tsv human_LR_pairs.txt
```

### About the parameters and input files
No_of_CT1_clusters (n): Number of clusters for cell type 1
No_of_CT2_clusters (m): Number of clusters for cell type 2

Table of expressed ligand or receptor gene in each cluster for a cell type. Tf the average expression level of a gene is greater than upper quantile in a cluster then define it as a expressed gene.
CT1_table_file: A talbe with 3 + n columns. The first column lists the gene name in target species. The second colum lists the gene name in human. The third column lists the gene type (1 for ligand and 2 for receptor). Then following by n columns of expression  
CT2_table_file: A talbe with 3 + m columns. The first column lists the gene name in target species. The second colum lists the gene name in human. The third column lists the gene type (1 for ligand and 2 for receptor).

L-R_pair_table: list of ligand-receptor pairs collected from public database, CellTalkDB (http://tcm.zju.edu.cn/celltalkdb/).  
