# Overview

This repository is dedicated to comparing the performance of the meta-assembler [**Aletsch**](https://github.com/Shao-Group/aletsch) against two leading meta-assemblers, [TransMeta](https://github.com/yutingsdu/TransMeta), [PsiCLASS](https://github.com/splicebox/PsiCLASS), as well as two meta-assembly pipelines based on single-sample assemblers, [StringTie2-Merge](https://ccb.jhu.edu/software/stringtie/index.shtml) and [Scallop2](https://github.com/Shao-Group/scallop2) combined with [TACO](https://tacorna.github.io). Here we provide instructions for downloading necessary tools, preparing datasets, executing the tools/pipelines, scoring Aletsch's output transcripts, and reproducing the results presented in the Aletsch paper.

# Step 1: Download and Link Tools

Our experiments involve the following tools:

| Tool                                                         | Version |                Description                |
| ------------------------------------------------------------ | :-----: | :---------------------------------------: |
| [Aletsch](https://github.com/Shao-Group/aletsch)             | v1.1.0  |              Meta Assembler               |
| [Transmeta](https://github.com/yutingsdu/TransMeta)          |  v.1.0  |              Meta Assembler               |
| [PsiCLASS](https://github.com/splicebox/PsiCLASS)            | v1.0.3  |              Meta Assembler               |
| [StringTie2](https://ccb.jhu.edu/software/stringtie/index.shtml) | v2.2.1  | Single-sample Assembler + GTF merged tool |
| [Scallop2](https://github.com/Shao-Group/scallop2)           | v1.1.2  |          Single-sample Assembler          |
| [TACO](https://tacorna.github.io)                            | v0.7.3  |              GTF merged tool              |
| [STAR](https://github.com/alexdobin/STAR/tree/master)        | v2.7.11 |              RNA-seq Aligner              |
| [GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml#gffcompare_dl) | v0.11.2 |      Evaluate assembled transcripts       |
| [gtfcuff](https://github.com/Kingsford-Group/rnaseqtools)    |    -    |               RNA-seq tool                |

#### Step 1.1: Download Tools

* Access the homepages of the respective tools using the links provided above.

- Follow the download and compilation instructions on each tool's homepage.

#### Step 1.2: Link or Copy Executables

- For tools with available executable files, link or copy them to the `programs` directory. This includes `aletsch`, `scallop2`, `taco_run`, `stringtie`, `STAR`, `gffcompare`, and `gtfcuff`.
- For tools without standalone executables (TransMeta and PsiCLASS), link the entire directory to `programs`.

Ensure the tools are accessible via the following paths:

```
your/path/to/programs/aletsch
your/path/to/programs/TransMeta/TransMeta
your/path/to/programs/psiclass/psiclass
your/path/to/programs/stringtie
your/path/to/programs/scallop2
your/path/to/programs/taco_run
your/path/to/programs/STAR
your/path/to/programs/gffcomapre
your/path/to/programs/gtfcuff
```

You may need to rename some of the executable files to match the paths listed above.

# Step 2: Download Datasets and Align

We evaluate the performance of the five methods using eight datasets, as outlined below. Each dataset is identified by its unique prefix (used in this repository) and accession ID for reference.

| Name in paper |      prefix in test(Ensembl)      |      Protocol       |                         Accession ID                         |
| :-----------: | :-------------------------------: | :-----------------: | :----------------------------------------------------------: |
|     BK-H1     |       **encode10_ensembl**        | Illumina paired-end |           Refer to<br /> `data/encode10.sra.list`            |
|     BK-H2     |      **PRJNA575230_ensembl**      | Illumina paired-end |          Refer to<br /> `data/PRJNA575230.sra.list`          |
|     SC-H1     |  **smartseq3_ensembl_human_100**  |      Smartseq3      | Random 100 cells from HEK293T of<br /> [E-MTAB-8735](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735) |
|     SC-H2     | **smartseq3_ensembl_xpress_run2** |  Smartseq3-Xpress   | All 1066 cells from PBMCs_run2 of<br /> [E-MTAB-11452](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-11452/sdrf) |
|     BK-H3     |      **PRJNA489891_ensembl**      | Illumina paired-end |          Refer to<br /> `data/PRJNA489891.sra.list`          |
|     SC-H3     |  **smartseq3_ensembl_human_92**   |      Smartseq3      | Random 92 cells from HEK293T of<br /> [E-MTAB-8735](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735) |
|     BK-M1     |      **PRJEB18790_ensembl**       | Illumina paired-end |          Refer to<br /> `data/PRJEB18790.sra.list`           |
|     SC-M1     |    **smartseq3_ensembl_mouse**    |      Smartseq3      | All 369 cells from Mouse-Fibroblast of<br /> [E-MTAB-8735](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735) |

Use [STAR](https://github.com/alexdobin/STAR/tree/master) for read alignments for each sample/cell. For every dataset, compile a list of all BAM file paths as required by the different meta-assemblers. Example for Aletsch: `data/encode10_ensembl.star.list`. Your lists should follow this format:

```
your/path/to/data/encode10_ensembl.star.list
your/path/to/data/PRJNA575230_ensembl.star.list
...
```

Annotations for evaluation are available from [Ensembl](http://useast.ensembl.org/Homo_sapiens/Info/Index) and [RefSeq](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/).

# Step 3: Run the Methods

Execute the provided scripts in the `results` directory to run the assemblers:

```
./run.aletsch.sh 
./run.transmeta.sh
./run.psiclass.sh
./run.stringtie2.sh
./run.st2merge.sh
./run.scallop.sh
./run.sc2taco.sh
```

Each script accepts different list arguments for the assembler. Example usage:

```
Usage: ./run.aletsch.sh {encode10|xx230|sc-human|sc-xpress|xx891|mouse790|sc-mouse} dataset_index
For example: ./run.aletsch.sh xx230 3
```

Note: Ensure `run.stringtie2.sh` is completed before `run.st2merge.sh`, and `run.scallop.sh` before `run.sc2taco.sh`. Scripts may be adjusted for parallel processing with specified CPU core usage. BAM files for specific chromosomes are prepared using [samtools](http://www.htslib.org/doc/samtools.html), e.g., for extracting chr1:

```
samtools view -b -F 4 -@ 10 -o "$chr1_bam" "$bam_file" 1
```

Alternatively, assemble with all chromosomes and later extract chr1 results from the GTF file.

# Step 4: Test with Aletsch's Pretrained Model and Compare

Navigate to the `random-forest` directory and run `prepare.sh` as done in Step 3, e.g., `./prepare.sh xx891 2`. This will collect Aletsch's feature files in `your/path/to/random-forest/aletsch`.

Download Aletsch's pretrained model from [Zenodo](https://doi.org/10.5281/zenodo.10602529), typically "Aletsch-ChrAll.ensembl.joblib". Execute the following command to test and compare the results:

```
cd your/path/to/random-forest
./prepare.sh xx891 2 # Example for dataset "PRJNA489891_ensembl"
python3 test.py -i xx891
```

Ensure you have the necessary Python packages (`numpy`, `pandas`, `sklearn`) installed. The `test.py` script can be modified to use a different pre-trained model, with outputs stored in `your/path/to/random-forest/logs`. An example log file, `Aletsch-Chr1.ensembl.sc-human100.other_chrs.log`, demonstrates the result of Aletsch-Chr1 on dataset SC-H1, including metrics like adjusted precision and pAUC for comparison.
