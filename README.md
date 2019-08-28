[![Build Status](https://travis-ci.org/jambler24/GenGraph.svg?branch=master)](https://travis-ci.org/jambler24/GenGraph)

# GenGraph
A repository for the GenGraph toolkit for the creation and manipulation of graph genomes

This project aims to make genome graphs simple to create, and provide tools for their manipulation and intergration into common workflows.

Please refer to the [wiki](https://github.com/jambler24/GenGraph/wiki) for installation and usage instructions. 

## Quickstart 

    python3 ./gengraphTool.py make_genome_graph --seq_file <sequence_file.txt> --out_file_name <filename> --recreate_check

A Docker file is also available that will contain everything you need for your first run. 



## Sequence file

The sequence file is a tab delimited file with 4 columns:

seq_name    |	aln_name	|   seq_path	|   annotation_path
------------ | ------------- | ------------- | -------------
H37Rv |	seq0 |	/Users/panix/Desktop/genomes/H37Rv/sequence.fasta |	NA
Beijing |	seq1 |	/Users/panix/Desktop/genomes/Beijing-NITR203/sequence.fasta |	NA
H37Ra |	seq2 |	/Users/panix/Desktop/genomes/H37Ra/sequence.fasta |	NA


The seq_name column is used as the node and edge 'ids' attribute value. It needs to be unique. 

