# Microbial-Mats-West-Antarctica
We employed a metagenomic approach to analyze 14 microbial mats from meltwater streams of western Antarctica, covering the Maritime, Peninsula, and Dry Valleys regions.

# Metagenomic Analysis Pipeline

This pipeline describes the bioinformatics workflow used for processing and analyzing metagenomic data from Antarctic microbial mats.

## Table of Contents
- [Prerequisites](#prerequisites)
- [1. Quality Control and Preprocessing](#1-quality-control-and-preprocessing)
- [2. Taxonomic Classification](#2-taxonomic-classification)
- [3. Diversity Analysis](#3-diversity-analysis)
- [4. Assembly](#4-assembly)
- [5. Gene Prediction](#5-gene-prediction)
- [6. Functional Annotation](#6-functional-annotation)
- [7. rRNA Analysis](#7-rrna-analysis)
- [8. Phylogenetic Reconstruction](#8-phylogenetic-reconstruction)
- [References](#references)

## Prerequisites

### Required Software
```bash
# Install required tools
conda install -c bioconda fastp kaiju megahit augustus prodigal mmseqs2 barrnap cd-hit mafft blast raxml
pip install eggnog-mapper
```

### Required Databases
```bash
# Download NCBI nr database for Kaiju (2021-02 version)
wget https://kaiju.binf.ku.dk/database/kaiju_db_nr_euk_2021-02-24.tgz
tar -xzf kaiju_db_nr_euk_2021-02-24.tgz

# Download eggNOG database
download_eggnog_data.py -y
```

## 1. Quality Control and Preprocessing

Process raw metagenomic FASTQ files for quality control and adapter removal [fastp](https://github.com/OpenGene/fastp).

```bash
#!/bin/bash

# Input files
RAW_R1="raw.1.fq.gz"
RAW_R2="raw.2.fq.gz"

# Output files
QC_R1="qc.1.fq.gz"
QC_R2="qc.2.fq.gz"

# Quality control with fastp
# - Remove adapters automatically
# - Filter low-quality reads
# - Remove duplicated reads  
# - Trim first 10 bp from both paired-end reads
fastp \
    -i ${RAW_R1} \
    -o ${QC_R1} \
    -I ${RAW_R2} \
    -O ${QC_R2} \
    --trim_front1 10 \
    --trim_front2 10 \
    --dedup \
    --dup_calc_accuracy 6 \
    --html fastp_report.html \
    --json fastp_report.json
```

## 2. Taxonomic Classification

Assign taxonomic classification using [Kaiju](https://github.com/bioinformatics-centre/kaiju) against NCBI nr database.

```bash
#!/bin/bash

# Kaiju database files
NODES_DMP="nodes.dmp"
KAIJU_DB="kaiju_db_nr_euk.fmi"
NAMES_DMP="names.dmp"

# Taxonomic classification with Kaiju
# Uses NCBI nr 2021-02 database with default parameters
kaiju \
    -t ${NODES_DMP} \
    -f ${KAIJU_DB} \
    -i ${QC_R1} \
    -j ${QC_R2} \
    -o kaiju.out \
    -z 4  # Number of parallel threads

# Generate taxonomic summary table at genus level
# Includes taxonomic hierarchy from superkingdom to genus
kaiju2table \
    -t ${NODES_DMP} \
    -n ${NAMES_DMP} \
    -r genus \
    -l superkingdom,phylum,class,order,family,genus \
    -o kaiju_summary.tsv \
    kaiju.out
```

## 3. Diversity Analysis

Profile alpha diversity using single marker genes and detect eukaryotic sequences with [SingleM](https://github.com/wwood/singlem), [EukDetect](https://github.com/allind/EukDetect), and [Metaxa2](https://microbiology.se/software/metaxa2/).

```bash
#!/bin/bash

# SingleM for alpha diversity profiling
# Generates OTU tables from single marker genes
singlem pipe \
    --forward ${QC_R1} \
    --reverse ${QC_R2} \
    --otu-table otu.table.tsv \
    --threads 4

# Eukaryote detection with EukDetect
# Detects eukaryotic sequences in metagenomic datasets
eukdetect \
    --mode mags \
    --reads ${QC_R1} ${QC_R2} \
    --output eukdetect_output

# Alternative eukaryote detection with Metaxa2
metaxa2 \
    -1 ${QC_R1} \
    -2 ${QC_R2} \
    -o metaxa2_output \
    --plus T \
    --graphical F
```

## 4. Assembly

Assemble quality-controlled reads into contigs using [MEGAHIT](https://github.com/voutcn/megahit) and [metaSPAdes](https://github.com/ablab/spades).

```bash
#!/bin/bash

# Assembly with MegaHit
# Uses default parameters for metagenomic assembly
megahit \
    -1 ${QC_R1} \
    -2 ${QC_R2} \
    -o mat.sample \
    --min-contig-len 1000 \
    -t 4  # Number of threads

# Filter contigs >= 1000 bp
CONTIGS="mat.sample/final.contigs.fa"

# Assembly with metaSPAdes for better bin recovery
metaspades.py \
    --only-assembler \
    -1 forward.paired-end.reads \
    -2 reverse.paired-end.reads \
    -o output.dir
```

## 5. Gene Prediction

Predict genes from assembled contigs based on taxonomic classification. For euk [Augustus](https://github.com/Gaius-Augustus/Augustus), prok [Prodigal](https://github.com/hyattpd/Prodigal); and ORFs getorf [EMBOSS](https://emboss.sourceforge.net/).

```bash
#!/bin/bash

# Separate contigs by taxonomic assignment
python separate_contigs.py  # Custom script to separate contigs

# Gene prediction for eukaryotic contigs using Augustus
augustus \
    --species=generic \
    --gff3=on \
    --outfile=eukaryotic_genes.gff3 \
    eukaryotic_contigs.fa

# Gene prediction for prokaryotic contigs using Prodigal  
prodigal \
    -i prokaryotic_contigs.fa \
    -a prokaryotic_genes.faa \
    -d prokaryotic_genes.fna \
    -o prokaryotic_genes.gff \
    -f gff \
    -p meta  # Metagenomic mode

# ORF prediction for unclassified contigs using EMBOSS getorf
getorf \
    -sequence unclassified_contigs.fa \
    -outseq unclassified_orfs.faa \
    -minsize 300 \
    -find 1  # Find all ORFs
```

## 6. Functional Annotation

Cluster predicted genes with [MMseqs2](https://github.com/soedinglab/MMseqs2) and perform functional annotation with [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper).

```bash
#!/bin/bash

# Combine all predicted genes
cat eukaryotic_genes.faa prokaryotic_genes.faa unclassified_orfs.faa > all_genes.faa

# Cluster genes to reduce redundancy
# 40% identity and 80% coverage threshold
mmseqs createdb all_genes.faa all_genes_db
mmseqs cluster all_genes_db clustered_db tmp \
    --min-seq-id 0.4 \
    --coverage-mode 1 \
    -c 0.8 \
    --cov-mode 0

mmseqs createseqfiledb all_genes_db clustered_db clustered_seqs
mmseqs result2flat all_genes_db all_genes_db clustered_seqs clustered_genes.faa

# Functional annotation with eggNOG-mapper
emapper.py \
    -i clustered_genes.faa \
    -o functional_annotation \
    --itype proteins \
    --data_dir /path/to/eggnog/data \
    --cpu 4 \
    --override

# Calculate GPM (Genes Per Million) values
# GPM = (counts of gene X / total number of paired reads) * 1,000,000
python calculate_gpm.py  # Custom script for GPM calculation
```

## 7. rRNA Analysis

Predict and analyze ribosomal RNA sequences [barrnap](https://github.com/tseemann/barrnap).

```bash
#!/bin/bash

# Predict rRNA sequences using Barrnap
barrnap \
    --kingdom bac \
    --outseq bacterial_rRNA.fasta \
    ${CONTIGS}

barrnap \
    --kingdom euk \
    --outseq eukaryotic_rRNA.fasta \
    ${CONTIGS}

# Annotate rRNA sequences against SILVA database
# Note: ACT Silva requires specific setup in 
ACT: Alignment, Classification and Tree Service
https://www.arb-silva.de/aligner

# Cluster eukaryotic rRNA sequences
cd-hit-est \
    -i eukaryotic_rRNA.fasta \
    -o eukaryotic_rRNA_clustered.fasta \
    -c 0.97 \
    -n 10 \
    -d 0 \
    -M 16000 \
    -T 4

# Align clustered sequences
mafft \
    --auto \
    eukaryotic_rRNA_clustered.fasta > eukaryotic_rRNA_aligned.fasta

# BLASTn annotation against nr database  
blastn \
    -query eukaryotic_rRNA_aligned.fasta \
    -db nt \
    -out rRNA_blast_results.txt \
    -outfmt 6 \
    -evalue 1e-5 \
    -num_threads 4
```

## 8. Phylogenetic Reconstruction

Perform phylogenetic analysis for *Adineta vaga*.

```bash
#!/bin/bash

# Extract Adineta vaga sequences (manual curation required)
grep -A1 "Adineta_vaga" eukaryotic_rRNA_aligned.fasta > adineta_sequences.fasta

# Add outgroup sequence (Seison nebaliae)
cat adineta_sequences.fasta seison_nebaliae.fasta > phylo_input.fasta

# Multiple sequence alignment with MAFFT
mafft \
    --auto \
    phylo_input.fasta > phylo_aligned.fasta

# Phylogenetic reconstruction with RAxML
raxml \
    -s phylo_aligned.fasta \
    -n adineta_tree \
    -m GTRGAMMAI \
    -p 12345 \
    -x 12345 \
    -# 1000 \
    -f a  # Rapid bootstrap analysis

# Visualize tree in FigTree (GUI application)
# Load RAxML_bipartitions.adineta_tree in FigTree
```

## Custom Scripts

### separate_contigs.py
```python
#!/usr/bin/env python3
"""
Separate contigs based on Kaiju taxonomic classification
"""

import sys
from Bio import SeqIO

def separate_contigs(contig_file, kaiju_output, output_prefix):
    """Separate contigs into eukaryotic, prokaryotic, and unassigned"""
    
    # Parse Kaiju output
    classifications = {}
    with open(kaiju_output, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[0] == 'C':  # Classified
                contig_id = parts[1]
                taxonomy = parts[2]
                classifications[contig_id] = taxonomy
    
    # Separate contigs
    euk_contigs = []
    prok_contigs = []
    unassigned_contigs = []
    
    for record in SeqIO.parse(contig_file, "fasta"):
        contig_id = record.id
        if contig_id in classifications:
            taxonomy = classifications[contig_id]
            if 'Eukaryota' in taxonomy:
                euk_contigs.append(record)
            else:
                prok_contigs.append(record)
        else:
            unassigned_contigs.append(record)
    
    # Write separated files
    SeqIO.write(euk_contigs, f"{output_prefix}_eukaryotic.fa", "fasta")
    SeqIO.write(prok_contigs, f"{output_prefix}_prokaryotic.fa", "fasta")  
    SeqIO.write(unassigned_contigs, f"{output_prefix}_unassigned.fa", "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python separate_contigs.py <contigs.fa> <kaiju.out> <output_prefix>")
        sys.exit(1)
    
    separate_contigs(sys.argv[1], sys.argv[2], sys.argv[3])
```

### calculate_gpm.py  
```python
#!/usr/bin/env python3
"""
Calculate GPM (Genes Per Million) values
GPM = (counts of gene X / total number of paired reads) * 1,000,000
"""

import pandas as pd
import sys

def calculate_gpm(gene_counts_file, total_reads, output_file):
    """Calculate GPM values for gene counts"""
    
    # Read gene counts
    df = pd.read_csv(gene_counts_file, sep='\t')
    
    # Calculate GPM
    df['GPM'] = (df['gene_counts'] / total_reads) * 1000000
    
    # Save results
    df.to_csv(output_file, sep='\t', index=False)
    print(f"GPM values calculated and saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python calculate_gpm.py <gene_counts.tsv> <total_reads> <output.tsv>")
        sys.exit(1)
    
    calculate_gpm(sys.argv[1], int(sys.argv[2]), sys.argv[3])
```

## References

Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of Molecular Biology, 215(3), 403–410. https://doi.org/10.1016/S0022-2836(05)80360-2

Bengtsson-Palme, J., Hartmann, M., Eriksson, K. M., Pal, C., Thorell, K., Larsson, D. G. J., & Nilsson, R. H. (2015). metaxa2: Improved identification and taxonomic classification of small and large subunit rRNA in metagenomic data. Molecular Ecology Resources, 15(6), 1403–1414. https://doi.org/10.1111/1755-0998.12399

Buchfink, B., Reuter, K., & Drost, H.-G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods, 18(4), 366–368. https://doi.org/10.1038/s41592-021-01101-x

Cantalapiedra, C. P., Hernández-Plaza, A., Letunic, I., Bork, P., & Huerta-Cepas, J. (2021). eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. Molecular Biology and Evolution, 38(12), 5825–5829. https://doi.org/10.1093/molbev/msab293

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

Huerta-Cepas, J., Szklarczyk, D., Heller, D., Hernández-Plaza, A., Forslund, S. K., Cook, H., Mende, D. R., Letunic, I., Rattei, T., Jensen, L. J., von Mering, C., & Bork, P. (2019). eggNOG 5.0: A hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Nucleic Acids Research, 47(D1), D309–D314. https://doi.org/10.1093/nar/gky1085

Hyatt, D., Chen, G.-L., LoCascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: Prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics, 11(1), 119. https://doi.org/10.1186/1471-2105-11-119

Li, D., Liu, C.-M., Luo, R., Sadakane, K., & Lam, T.-W. (2015). MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674–1676. https://doi.org/10.1093/bioinformatics/btv033

Li, W., & Godzik, A. (2006). Cd-hit: A fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics, 22(13), 1658–1659. https://doi.org/10.1093/bioinformatics/btl158

Lind, A. L., & Pollard, K. S. (2021). Accurate and sensitive detection of microbial eukaryotes from whole metagenome shotgun sequencing. Microbiome, 9(1), 58. https://doi.org/10.1186/s40168-021-01015-y

Menzel, P., Ng, K. L., & Krogh, A. (2016). Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nature Communications, 7(1), 11257. https://doi.org/10.1038/ncomms11257

Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J., & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. Nucleic Acids Research, 41(D1), D590–D596. https://doi.org/10.1093/nar/gks1219

Rice, P., Longden, I., & Bleasby, A. (2000). EMBOSS: The European Molecular Biology Open Software Suite. Trends in Genetics, 16(6), 276–277. https://doi.org/10.1016/S0168-9525(00)02024-2

Sayers, E. W., Bolton, E. E., Brister, J. R., Canese, K., Chan, J., Comeau, D. C., Connor, R., Funk, K., Kelly, C., Kim, S., Madej, T., Marchler-Bauer, A., Lanczycki, C., Lathrop, S., Lu, Z., Thibaud-Nissen, F., Murphy, T., Phan, L., Skripchenko, Y., … Sherry, S. T. (2022). Database resources of the national center for biotechnology information. Nucleic Acids Research, 50(D1), D20–D26. https://doi.org/10.1093/nar/gkab1112

Seemann, T. (2018). BAsic Rapid Ribosomal RNA Predictor. Barrnap [Computer software]. https://github.com/tseemann/barrnap

Stamatakis, A. (2014). RAxML version 8: A tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30(9), 1312–1313. https://doi.org/10.1093/bioinformatics/btu033

Stanke, M., Schöffmann, O., Morgenstern, B., & Waack, S. (2006). Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources. BMC Bioinformatics, 7(1), 62. https://doi.org/10.1186/1471-2105-7-62
