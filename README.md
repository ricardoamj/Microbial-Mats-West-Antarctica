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
- [9. Recovery bins](#9-recovery-bins)
- [Expect Results](#expect_results)
- [References](#references)


## Prerequisites

We run our analysis in server Ubuntu 24.04.3 LTS, 32 core Intel(R) Xeon(R) CPU E5-2640 v2 @ 2.00GHz, 256 GB RAM

- Miniconda (conda 25.7.0)
- FastQC v0.12.1
- fastp 0.24.0
- Kaiju 1.10.1
- SingleM v0.19.0
- EukDetect v1.0.1
- Metaxa2 v2.2.3
- MEGAHIT v1.2.9
- SPAdes genome assembler v3.15.4 [metaSPAdes mode]
- AUGUSTUS (3.5.0)
- PRODIGAL v2.6.3
- getorf EMBOSS:6.6.0.0
- MMseqs2 Version: 1046260d43f8d721041dec43a1763ecc450a6ea9
- emapper-2.0.1
- barrnap 0.9
- CD-HIT version 4.8.1
- MAFFT v7.525
- blast 2.16.0
- RAxML version 8.2.12
- MaxBin 2.2.7
- concoct 1.1.0
- MetaBAT2 version 2:2.18
- CheckM v1.2.4
- DAS Tool 1.1.7


### Required Software
Install conda

If you do not have the conda package manager installed already, follow the instructions to install [Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

```bash
# Install required tools
mamba create -n meta_pipeline -c bioconda fastp kaiju megahit augustus prodigal mmseqs2 barrnap cd-hit blast raxml spades eggnog-mapper fastqc metaeuk seqtk  
# Clean shell
mamba deactivate
conda deactivate
source ~/.bashrc
mamba activate meta_pipeline

# Install mafft from source
wget https://mafft.cbrc.jp/alignment/software/mafft-7.525-with-extensions-src.tgz
# follow the indication for non-root user [here](https://mafft.cbrc.jp/alignment/software/installation_without_root.html)

# Install SingleM in our environment to resolve meta_pipeline environment conflicts
mamba create -c conda-forge -c bioconda --override-channels --name singlem singlem'>='0.19.0
# Again clean shell
mamba deactivate
conda deactivate
source ~/.bashrc
mamba activate singlem

# Install EukDetect
# Please redirect to [EukDetect install](https://github.com/allind/EukDetect?tab=readme-ov-file#installation) follow the indications
# We create a env in miniconda 
# Run in folder install EukDetect
python setup.py install
# Modify default_configfile.yml to run your datasets

# Install Metaxa2 from https://microbiology.se/software/metaxa2/

# Create binning environment in conda 
mamba create -c conda-forge -c bioconda --override-channels --name binning spades concoct maxbin2 metabat2 quast
conda activate binning
pip install drep
mamba install hmmer prodigal pplacer
pip3 install numpy
pip3 install matplotlib
pip3 install pysam
pip3 install checkm-genome
```

### Required Databases
```bash
# Download NCBI nr database for Kaiju
# We recomended download pre-built from https://bioinformatics-centre.github.io/kaiju/downloads.html
# for Subset of NCBI BLAST nr database containing Archaea, bacteria, viruses; additionally including fungi and microbial eukaryotes

wget https://kaiju-idx.s3.eu-central-1.amazonaws.com/2021/kaiju_db_nr_euk_2021-02-24.tgz

# uncompress tar gz in where a you need it
tar -xzf kaiju_db_nr_euk_2021-02-24.tgz

# Download eggNOG database
download_eggnog_data.py -y

# Storage Requirements for eggnog_data

    # ~40 GB for the eggNOG annotation databases (eggnog.db and eggnog.taxa.db)
    # ~9 GB for Diamond database of eggNOG sequences (required if using -m diamond, which is the default search mode).
    # ~11 GB for MMseqs2 database of eggNOG sequences (~86 GB if MMseqs2 index is created) (required if using -m mmseqs).
    # ~3 GB for PFAM database (required if using --pfam_realign options for realignment of queries to PFAM domains).
    # The size of eggNOG diamond/mmseqs databases created with create_dbs.py is highly variable, depending on the size of the chosen taxonomic groups.

# For singlem you need
mamba activate singlem
singlem data --output-directory singlem_data

# For checkm you need this database
 wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz


```

## 1. Quality Control and Preprocessing

Process raw metagenomic FASTQ files for quality control [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and adapter removal [fastp](https://github.com/OpenGene/fastp).

```bash
#!/bin/bash

# Input files
RAW_R1="raw.1.fq.gz"
RAW_R2="raw.2.fq.gz"

# Output files
QC_R1="qc.1.fq.gz"
QC_R2="qc.2.fq.gz"

# Quality check with FastQC
fastqc \
    ${RAW_R1}\
    ${RAW_R2}\
    -o FastQC

# Check the FastQC Report for your raw_sequences: basic statistics, Per base sequence quality, Adapter Content

# Quality control with fastp
# - Remove adapters automatically
# - Filter low-quality reads
# - Remove duplicated reads  
# - (In this case) Trim first 10 bp from both paired-end reads
fastp \
    -i ${RAW_R1} \
    -o ${QC_R1} \
    -I ${RAW_R2} \
    -O ${QC_R2} \
    --trim_front1 10 \
    --trim_front2 10 \
    --html fastp_report.html \
    --json fastp_report.json

# Check the new report for quality control reads in fastp_report.html
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
# pre-trimming reads not recommended

kaiju \
    -t ${NODES_DMP} \
    -f ${KAIJU_DB} \
    -i ${RAW_R1} \
    -j ${RAW_R2} \
    -o kaiju.out \
    -z 16  # Number of parallel threads that you have available

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

Profile alpha diversity using single marker genes [SingleM](https://github.com/wwood/singlem), and detect eukaryotic sequences with [EukDetect](https://github.com/allind/EukDetect), auxiliary program [Metaxa2](https://microbiology.se/software/metaxa2/).

```bash
#!/bin/bash

# SingleM for alpha diversity profiling
# Generates OTU tables from single marker genes
# Recomended use raw reads
# export database variable
export SINGLEM_METAPACKAGE_PATH='/home/user/singlem_data/S5.4.0.GTDB_r226.metapackage_20250331.smpkg.zb
conda activate singlem
# run singlem 
singlem pipe \
    -1 ${RAW_R1} \
    -2 ${RAW_R2} \
    --otu-table otu.table.tsv \
    --threads 4 # Number of parallel threads that you have available

# Eukaryote detection with EukDetect
# Detects eukaryotic sequences in metagenomic datasets

conda activate eukdetect
eukdetect \
    --mode runall \
    --configfile [config file] \
    --cores [cores] # Number of parallel threads that you have available

# Alternative eukaryote detection with Metaxa2
metaxa2 \
    -1 ${QC_R1} \ #DNA FASTQ input file containing the first reads in the read pairs to investigate
    -2 ${QC_R2} \ #DNA FASTQ input file containing the second reads in the pairs to investigate
    -o metaxa2_output \ #output file
    --format p \ #Specifies the format of the input file, p paired end
    --plus T \ #Runs blast search through blast+ instead of the legacy blastall engine
    -t e \ #Profile set to use for the search e eukaryota
```

## 4. Assembly

Assemble quality-controlled reads into contigs using [MEGAHIT](https://github.com/voutcn/megahit) and [metaSPAdes](https://github.com/ablab/spades).

Key advantages of MEGAHIT is Computational efficiency: MegaHit is optimized to use memory very efficiently, allowing large datasets to be assembled.

Advantages of metaSPAdes: Low-abundance species: Better for detecting and assembling rare organisms in the community.
Complex downstream analyses: When you need high-quality contigs for functional annotation, binning, or phylogenetic analyses.


```bash
#!/bin/bash

# Activate conda environment
conda activate meta_pipeline

# Assembly with MegaHit
# Uses default parameters for metagenomic assembly
megahit \
    -1 ${QC_R1} \ #comma-separated list of fasta/q paired-end #1 files, paired with files in <pe2>
    -2 ${QC_R2} \ #comma-separated list of fasta/q paired-end #2 files, paired with files in <pe2>
    -o mat.sample \ #out dir
    --min-contig-len 1000 \ #minimum length of contigs to output
    -t 4  # number of CPU threads [# of logical processors]

# Filter contigs >= 1000 bp
CONTIGS="mat.sample/final.contigs.fa"

# Assembly with metaSPAdes for better bin recovery
metaspades.py \
    --only-assembler \ #runs only assembling (without read error correction)
    -1 forward.paired-end.reads \ #file with forward paired-end reads
    -2 reverse.paired-end.reads \ #file with reverse paired-end reads
    -o output.dir \ #output dir 
    --threads 16 \ # <int> number of threads. [default: 16]
```

## 5. Gene Prediction

Predict genes from assembled contigs based on taxonomic classification. For euk [Augustus](https://github.com/Gaius-Augustus/Augustus), prok [Prodigal](https://github.com/hyattpd/Prodigal); and ORFs getorf [EMBOSS](https://emboss.sourceforge.net/).

```bash
#!/bin/bash

# Filter contigs with less 1000 bp
# For example with seqtk
seqtk seq -L 1000 input.contigs.fasta > output.gt1k.contigs.fasta

# Separate contigs by taxonomic assignment with kaiju

kaiju \
    -t nodes.dmp \
    -f kaiju_db_nr_euk.fmi \
    -i output.gt1k.contigs.fasta \
    -o kaiju.classification \
    -e 2 -m 15 -s 70 -E 1e-06 \
        #run mode: Greedy
         # minimum match length: 15
         # seed length: 7
         # minimum blosum62 score for matches: 70
         # minimum E-value: 1e-06
         # max number of mismatches within a match: 2
    -v \ #verbose mode
    -z 16 # number of threads

# Annotate taxonomic name from kaiju.classification

 kaiju-addTaxonNames \
    -i kaiju.classification \
    -o output.taxon.names \
    -t nodes.dmp \
    -n names.dmp \
    -u \ # Unclassified are not contained in the output
    -r superkingdom,phylum,class,order,family,genus \ #Print taxon path containing only ranks specified by a comma-separated list

# Filter list contigs by Eukaryota, Bacteria or Unclass
cat output.taxon.names | grep Eukaryota | cut -f2 | sort -u > list.contigs.Eukaryota
cat output.taxon.names | egrep 'Archaea|Bacteria' | cut -f2 | sort -u > list.contigs.Prokaryota
cut -f2 kaiju.classification | grep -v -F -w -f <(cat list.contigs.Prokaryota list.contigs.Eukaryota) > list.contigs.unclassified

# Subseq contigs from lists
seqtk subseq contigs.gt1k.fasta list.contigs.Prokaryota > contigs.gt1k.prok.fasta
seqtk subseq contigs.gt1k.fasta list.contigs.Eukaryota > contigs.gt1k.euk.fasta
seqtk subseq contigs.gt1k.fasta list.contigs.unclassified > contigs.gt1k.unclass.fasta

# Gene prediction for eukaryotic contigs using Augustus
augustus \
    --species=<your favorite model> \
    --gff3=on \
    --outfile=eukaryotic_genes.gff3 \
    contigs.gt1k.euk.fasta

# Gene prediction for prokaryotic contigs using Prodigal  
prodigal \
    -i contigs.gt1k.prok.fasta \
    -a prokaryotic_genes.faa \
    -d prokaryotic_genes.fna \
    -o prokaryotic_genes.gff \
    -f gff \
    -p meta  # Metagenomic mode

# ORF prediction for unclassified contigs using EMBOSS getorf
getorf \
    -sequence contigs.gt1k.unclass.fasta \
    -outseq unclassified_orfs.faa \
    -minsize 150 \
    -find 1  # Find all ORFs
```

## 6. Functional Annotation

Cluster predicted genes with [MMseqs2](https://github.com/soedinglab/MMseqs2) and perform functional annotation with [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper).

```bash
#!/bin/bash

# Cluster prediction genes to reduce redundancy
# 40% identity and 80% coverage threshold
mmseqs createdb predicted_genes.faa all_genes_db
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
    -evalue 1e-10 \
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
## 9. Recovery bins

We will obtain the bins for each of the metagenomes using [Maxbin2](https://sourceforge.net/projects/maxbin2/), [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/src/master/), and [CONCOCT](https://github.com/BinPro/CONCOCT).

Then we will use the [DAS tool](https://github.com/cmks/DAS_Tool)  to optimize redundant bins in each assembly.

On the other hand, we will review the quality of the bins using [checkm](https://github.com/Ecogenomics/CheckM). 

And subsequently, if warranted (Optional), we will dereplicate bins that come from the same type of microorganism (at the species or strain level) with [dRep](https://github.com/MrOlm/drep).


```bash
#!/bin/bash

# Map reads to contigs for coverage calculation
bowtie2-build ${CONTIGS} contigs_index
bowtie2 -x contigs_index -1 ${QC_R1} -2 ${QC_R2} -S alignment.sam
samtools view -bS alignment.sam | samtools sort -o alignment_sorted.bam
samtools index alignment_sorted.bam

# Generate depth file for MetaBAT2
jgi_summarize_bam_contig_depths --outputDepth depth.txt alignment_sorted.bam

# Binning with MetaBAT2
metabat2 -i ${CONTIGS} -a depth.txt -o metabat2_bins/bin

# Binning with MaxBin2
run_MaxBin.pl -contig ${CONTIGS} -reads ${QC_R1} -reads2 ${QC_R2} -out maxbin2_bins/bin

# Binning with CONCOCT
# [Agregar comandos de CONCOCT]

# Optimize bins with DAS Tool
DAS_Tool -i metabat2_bins.txt,maxbin2_bins.txt -l MetaBAT2,MaxBin2 -c ${CONTIGS} -o das_tool_output

# Quality assessment with CheckM
checkm lineage_wf das_tool_output checkm_output

```

## Expected Outputs
- Quality reports (FastQC, fastp)
- Taxonomic profiles (Kaiju tables)
- Diversity metrics (SingleM OTU tables)
- Assembled contigs and scaffolds
- Predicted genes and functional annotations
- High-quality genome bins
- Phylogenetic trees


## References

Alneberg, J., Bjarnason, B. S., de Bruijn, I., Schirmer, M., Quick, J., Ijaz, U. Z., Lahti, L., Loman, N. J., Andersson, A. F., & Quince, C. (2014). Binning metagenomic contigs by coverage and composition. Nature Methods, 11(11), 1144–1146. https://doi.org/10.1038/nmeth.3103

Altschul, S. F., Gish, W., Miller, W., Myers, E. W., & Lipman, D. J. (1990). Basic local alignment search tool. Journal of Molecular Biology, 215(3), 403–410. https://doi.org/10.1016/S0022-2836(05)80360-2

Bengtsson-Palme, J., Hartmann, M., Eriksson, K. M., Pal, C., Thorell, K., Larsson, D. G. J., & Nilsson, R. H. (2015). metaxa2: Improved identification and taxonomic classification of small and large subunit rRNA in metagenomic data. Molecular Ecology Resources, 15(6), 1403–1414. https://doi.org/10.1111/1755-0998.12399

Buchfink, B., Reuter, K., & Drost, H.-G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods, 18(4), 366–368. https://doi.org/10.1038/s41592-021-01101-x

Cantalapiedra, C. P., Hernández-Plaza, A., Letunic, I., Bork, P., & Huerta-Cepas, J. (2021). eggNOG-mapper v2: Functional Annotation, Orthology Assignments, and Domain Prediction at the Metagenomic Scale. Molecular Biology and Evolution, 38(12), 5825–5829. https://doi.org/10.1093/molbev/msab293

Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890. https://doi.org/10.1093/bioinformatics/bty560

Huerta-Cepas, J., Szklarczyk, D., Heller, D., Hernández-Plaza, A., Forslund, S. K., Cook, H., Mende, D. R., Letunic, I., Rattei, T., Jensen, L. J., von Mering, C., & Bork, P. (2019). eggNOG 5.0: A hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Nucleic Acids Research, 47(D1), D309–D314. https://doi.org/10.1093/nar/gky1085

Hyatt, D., Chen, G.-L., LoCascio, P. F., Land, M. L., Larimer, F. W., & Hauser, L. J. (2010). Prodigal: Prokaryotic gene recognition and translation initiation site identification. BMC Bioinformatics, 11(1), 119. https://doi.org/10.1186/1471-2105-11-119

Kang, D. D., Li, F., Kirton, E., Thomas, A., Egan, R., An, H., & Wang, Z. (2019). MetaBAT 2: An adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. PeerJ, 7, e7359. https://doi.org/10.7717/peerj.7359

Li, D., Liu, C.-M., Luo, R., Sadakane, K., & Lam, T.-W. (2015). MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. Bioinformatics, 31(10), 1674–1676. https://doi.org/10.1093/bioinformatics/btv033

Li, W., & Godzik, A. (2006). Cd-hit: A fast program for clustering and comparing large sets of protein or nucleotide sequences. Bioinformatics, 22(13), 1658–1659. https://doi.org/10.1093/bioinformatics/btl158

Lind, A. L., & Pollard, K. S. (2021). Accurate and sensitive detection of microbial eukaryotes from whole metagenome shotgun sequencing. Microbiome, 9(1), 58. https://doi.org/10.1186/s40168-021-01015-y

Menzel, P., Ng, K. L., & Krogh, A. (2016). Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nature Communications, 7(1), 11257. https://doi.org/10.1038/ncomms11257

Olm, M. R., Brown, C. T., Brooks, B., & Banfield, J. F. (2017). dRep: A tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication. The ISME Journal, 11(12), 2864–2868. https://doi.org/10.1038/ismej.2017.126

Parks, D. H., Imelfort, M., Skennerton, C. T., Hugenholtz, P., & Tyson, G. W. (2015). CheckM: Assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome Research, 25(7), 1043–1055. https://doi.org/10.1101/gr.186072.114

Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J., & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. Nucleic Acids Research, 41(D1), D590–D596. https://doi.org/10.1093/nar/gks1219

Rice, P., Longden, I., & Bleasby, A. (2000). EMBOSS: The European Molecular Biology Open Software Suite. Trends in Genetics, 16(6), 276–277. https://doi.org/10.1016/S0168-9525(00)02024-2

Sayers, E. W., Bolton, E. E., Brister, J. R., Canese, K., Chan, J., Comeau, D. C., Connor, R., Funk, K., Kelly, C., Kim, S., Madej, T., Marchler-Bauer, A., Lanczycki, C., Lathrop, S., Lu, Z., Thibaud-Nissen, F., Murphy, T., Phan, L., Skripchenko, Y., … Sherry, S. T. (2022). Database resources of the national center for biotechnology information. Nucleic Acids Research, 50(D1), D20–D26. https://doi.org/10.1093/nar/gkab1112

Seemann, T. (2018). BAsic Rapid Ribosomal RNA Predictor. Barrnap [Computer software]. https://github.com/tseemann/barrnap

Sieber, C. M. K., Probst, A. J., Sharrar, A., Thomas, B. C., Hess, M., Tringe, S. G., & Banfield, J. F. (2018). Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. Nature Microbiology, 3(7), 836–843. https://doi.org/10.1038/s41564-018-0171-1

Stamatakis, A. (2014). RAxML version 8: A tool for phylogenetic analysis and post-analysis of large phylogenies. Bioinformatics, 30(9), 1312–1313. https://doi.org/10.1093/bioinformatics/btu033

Wu, Y.-W., Simmons, B. A., & Singer, S. W. (2016). MaxBin 2.0: An automated binning algorithm to recover genomes from multiple metagenomic datasets. Bioinformatics, 32(4), 605–607. https://doi.org/10.1093/bioinformatics/btv638


Stanke, M., Schöffmann, O., Morgenstern, B., & Waack, S. (2006). Gene prediction in eukaryotes with a generalized hidden Markov model that uses hints from external sources. BMC Bioinformatics, 7(1), 62. https://doi.org/10.1186/1471-2105-7-62
