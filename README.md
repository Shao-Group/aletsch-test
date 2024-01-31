# Overview

This repository compares the performance of meta-assembler [**Aletsch**](https://github.com/Shao-Group/aletsch) with two leading meta-assemblers, [TransMeta](https://github.com/yutingsdu/TransMeta), [PsiCLASS](https://github.com/splicebox/PsiCLASS), and two additional meta-assembly pipelines based on single-sample assemblers, [StringTie2-Merge](https://ccb.jhu.edu/software/stringtie/index.shtml) and [Scallop2](https://github.com/Shao-Group/scallop2)+[TACO](https://tacorna.github.io). Here we provide instructions to download and link all necessary tools, prepare input lists for datasets, run all the five tools/pipelines, evaluate output transcripts and reproduce results in Aletsch paper:

1. Download necessary tools (in `program`)
2. Download datasets and prepare bam lists(in `data`)
3. Run the methods to produce results(in `results`)
4. Test Aletsch's pre-trained model and compare results(in `random-forest`)

# Step 1: Download and Link Tools

Our experiments involve the following tools. Users need to separately download all neccessary tools and link them to the folder `programs` before running experiments.

| Tool                                                         | Version |                Description                |
| ------------------------------------------------------------ | :-----: | :---------------------------------------: |
| [Aletsch](https://github.com/Shao-Group/aletsch)             | v1.1.0  |              Meta Assembler               |
| [Transmeta](https://github.com/yutingsdu/TransMeta)          |  v.1.0  |              Meta Assembler               |
| [PsiCLASS](https://github.com/splicebox/PsiCLASS)            | v1.0.3  |              Meta Assembler               |
| [StringTie2]((https://ccb.jhu.edu/software/stringtie/index.shtml)) | v2.2.1  | Single-sample Assembler + GTF merged tool |
| [Scallop2](https://github.com/Shao-Group/scallop2)           | v1.1.2  |          Single-sample Assembler          |
| [TACO](https://tacorna.github.io)                            | v0.7.3  |              GTF merged tool              |
| [STAR](https://github.com/alexdobin/STAR/tree/master)        | v2.7.11 |              RNA-seq Aligner              |
| [GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml#gffcompare_dl) | v0.11.2 |      Evaluate assembled transcripts       |
| [gtfcuff](https://github.com/Kingsford-Group/rnaseqtools)    |    -    |               RNA-seq tool                |

**Step 1.1**: Click above tools and those hyperlinks will navigate users to the homepage of tools. Then please follow the instructions provided in tools' homepages to download and/or compile above tools.

**Step 1.2**: Please link or copy the executable files to `programs` directory if they are avaliable (aletsch, scallop2, taco_run, stringtie2, STAR, gffcompare, gtfcuff). Otherwise please link the directory here (TransMeta and PsiCLASS).

Please make sure tools can be called by the following paths:

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

You may need to rename some executable files to make sure that all tools can be called successfully by above paths. 

# Step 2: Download Datasets and Align

We compare the five methods on eight datasets:

| Name in paper |      Protocol       |                         Accession ID                         |      |
| :-----------: | :-----------------: | :----------------------------------------------------------: | :--: |
|     BK-H1     | Illumina paired-end |              Refer to `data/encode10.sra.list`               |      |
|     BK-H2     | Illumina paired-end |             Refer to `data/PRJNA575230.sra.list`             |      |
|     SC-H1     |      Smartseq3      | Random 100 cells from HEK293T of  [E-MTAB-8735](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735) |      |
|     SC-H2     |  Smartseq3-Xpress   | All 1066 cells from run2 of [E-MTAB-11452](http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-11452/) |      |
|     BK-H3     | Illumina paired-end |             Refer to `data/PRJNA489891.sra.list`             |      |
|     SC-H3     |      Smartseq3      | Random 92 cells from HEK293T of [E-MTAB-8735](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735) |      |
|     BK-M1     | Illumina paired-end |             Refer to `data/PRJEB18790.sra.list`              |      |
|     SC-M1     |      Smartseq3      | All 369 cells from Mouse-Fibroblast of [E-MTAB-8735](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8735) |      |

For each sample/cell, we use [STAR](https://github.com/alexdobin/STAR/tree/master) to generate read alignments. For each dataset, generate a list about all paths of  bam files, according to requirements of different meta-assemblers. We provide an example of bam list fed into Aletsch: `data/encode10_ensembl.star.list`. Please prepare all your lists similar to the following:

```
your/path/to/data/encode10_ensembl.star.list
your/path/to/data/PRJNA575230_ensembl.star.list
your/path/to/data/smartseq3_ensembl_human.star.list
your/path/to/data/smartseq3_ensembl_xpress_run2.star.list
...
```

Besides, we also need the annotation files for evaluation: [Ensembl](http://useast.ensembl.org/Homo_sapiens/Info/Index) and [RefSeq](https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/).

# Step 3: Run the Methods

Once the datasets and programs are all available, use the following scripts in `results`
to run the assemblers on the datasets:

```
./run.aletsch.sh 
./run.transmeta.sh
./run.psiclass.sh
./run.stringtie2.sh
./run.st2merge.sh
./run.scallop.sh
./run.sc2taco.sh
```

For each script, you can choose different lists fed into the assembler, for example:

```
Usage: ./run.aletsch.sh {encode10|xx230|sc-human|sc-mouse|sc-xpress} dataset_index
For example: ./run.aletsch.sh xx230 3
Dataset Choices:
  For encode10(encode10):
  0: encode10_ensembl_chr1
  1: encode10_ensembl_other_chrs
  2: encode10_ensembl
  3: encode10_refseq_chr1
  4: encode10_refseq_other_chrs
  5: encode10_refseq
  all: submit all choices in this group
  ...
```

Note that `run.stringtie2.sh` should be finished before `run.st2merge.sh` , and `run.scallop.sh` should be finished before `run.sc2taco.sh`. You can modify the scripts to specify how many CPU cores you want to use to run the jobs in parallel. 

# Step 4: Test on Aletsch's Pretrained Model and Compare

Go to directory "random-forst" and run `prepare.sh` similar to the scripts in Step3, for example: `./prepare.sh sc-human 2`. Then all Aletsch's feature files will be collected in `your/path/to/random-forest/aletsch`. Run the following command to test and compare results of Aletsch and the other tools:

```
python3 test.py
```

Note that basic package `numpy`, `pandas` and `sklearn` should be installed. The output of `test.py` are saved in `your/path/to/random-forest/logs`. All logs under Ensembl annotation have been saved.
