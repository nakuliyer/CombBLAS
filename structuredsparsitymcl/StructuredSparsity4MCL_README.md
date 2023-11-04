# Structured Sparsity for Marcov Cluster Algorithm

## Introduction
The Markov cluster algorithm (MCL) performs flow-based clustering on a graph, often constructed to represent similarities among data points. MCL (and its distributed version HipMCL) has high memory footprint and high communication cost when running at large concurrencies, which makes it hard to apply to extremely large datasets. We are investigating how to utilize the structure of the similarity matrix to accelerate the computation and reduce the memory footprint.

In this work, we utilize the cluster information from the similarity matrix and reshuffle the rows and columns so that vertices belonging to the same cluster are close to each other. This will "push" the non-zero elements towards the diagonal. We can further split the matrix into two parts: 1) a diagonal block that contains many more non-zero elements (``dense'') and 2) an off-diagonal part that contains fewer non-zero elements.

TODO: add Citation

## Pipelines
Record the ideas and related files.

### First, where is the dataset? Format?
[dataset portal](https://portal.nersc.gov/project/m1982/HipMCL/viruses/)
We are look into 3 datasets, virus, archaea, eukarya.
We use index.triple format as input.

#### convert index.triple format to mtx format.
why? because CombBLAS loads mtx format faster than index.triple format.
**mtx format is always one based in this research.**

`
python convertindextriple2mtx.py InputIndexTripleFile 
`


output is `$InputIndexTripleFile.mtx`
Example input and output 

```
(base) yuxihong@perlmutter:login02:/pscratch/sd/y/yuxihong/hipmcl_dataset> head vir_vs_vir_30_50length.indexed.triples
isolate_virus|2546825548|2547000401|	isolate_virus|2546825788|2547020166|	1
isolate_virus|2546825548|2547000402|	isolate_virus|2546825548|2547000403|	0.4709
isolate_virus|2546825548|2547000403|	isolate_virus|2546825548|2547000402|	0.4709
isolate_virus|2546825548|2547000404|	isolate_virus|2546825571|2547004679|	0.3448
isolate_virus|2546825548|2547000404|	isolate_virus|2546825578|2547005500|	0.3556
isolate_virus|2546825548|2547000404|	isolate_virus|2546825579|2547005599|	0.3678
isolate_virus|2546825548|2547000404|	isolate_virus|2546825605|2547008174|	0.3333
isolate_virus|2546825548|2547000404|	isolate_virus|2546825689|2547014668|	0.3769
isolate_virus|2546825548|2547000404|	isolate_virus|2546825839|2547024240|	0.3571
isolate_virus|2546825548|2547000404|	isolate_virus|2546825937|2547031834|	0.3333
```

```
(base) yuxihong@perlmutter:login02:/pscratch/sd/y/yuxihong/hipmcl_dataset> head vir_vs_vir_30_50length.indexed.triples.mtx
%%MatrixMarket matrix coordinate real general
219715 219715 4583048
1	2	1
3	4	0.4709
4	3	0.4709
5	6	0.3448
5	7	0.3556
5	8	0.3678
5	9	0.3333
5	10	0.3769
```

check the maximum value and minmum value of row and col of output. 
- min row `awk 'NR > 2 && ($1 < min || NR == 3) {min=$1} END {print min}' filename` : should be 1
- max row `awk 'NR > 2 && ($1 > max || NR == 3) {max=$1} END {print max}' filename` : should be 219715
- min row `awk 'NR > 2 && ($2 < min || NR == 3) {min=$2} END {print min}' filename` : should be 1
- max row `awk 'NR > 2 && ($2 > max || NR == 3) {max=$2} END {print max}' filename` : should be 219715

#### 



### HipMCL: The baseline
```
srun --constraint=cpu --ntasks=4 -N 1 -c 32 --time=00:05:00 --qos debug \
/global/homes/y/yuxihong/graphclustering/CombBLAS/build/Applications/mcl \
-M /pscratch/sd/y/yuxihong/hipmcl_dataset/vir_vs_vir_30_50length.indexed.triples.mtx -I 2 -per-process-mem 64
```



### Graph Partition Using METIS 
> Aydin Oct 18, 2023 through email: As for technical stuff, I was thinking of using graph partitioning to find some of these initial clusters you are defining your blocks from. I noticed that the virus data is perfectly partitionable for example and this might be the case for other datasets as well.

To generate the initial block forms, we can use METIS to partition the graph.

TODO: What's the time complexity of the METIS partition algorithms?

#### install METIS 
git repo: https://github.com/KarypisLab/METIS

#### graph partioner binary
graph partioner:
- gpmetis [options] graphfile nparts

#### convert index.triple format to METIS graph format
code in: `convert2metisinput.py`

Launch command: 
`python convert2metisinput.py /pscratch/sd/y/yuxihong/hipmcl_dataset/vir_vs_vir_30_50length.indexed.triples`
output `/pscratch/sd/y/yuxihong/hipmcl_dataset/vir_vs_vir_30_50length.indexed.triples.graph`
#### run partioner 
partition in 10 parts: 
`gpmetis vir_vs_vir_30_50length.indexed.triples.graph 10`

#### shuffle the input matrix 

## some utility function

### Check input matrix is symmetric.

