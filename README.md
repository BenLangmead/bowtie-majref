This collection of scripts are used to build [1KGenomes project](http://www.internationalgenome.org)-informed major-allele reference genomes for Bowtie 1 and Bowtie 2.  The current scripts build *linear* references with major alleles. We provide SNP-only and SNP-and-indel indexes for both Bowtie versions. The inclusion of indels change the genomic coordinate system. We recommend using [levioSAM](https://github.com/alshai/levioSAM) to perform accurate and scalable lift-over for alignments against the SNP-and-indel indexes.

### Downloads

#### SNP-only 
Pre-built major-allele-SNP reference indexes are available:

| Aligner  | Reference                     | Index zip                                                            |
|----------|-------------------------------|----------------------------------------------------------------------|
| Bowtie 2 | GRCh38 + major SNPs           | [https](https://genome-idx.s3.amazonaws.com/bt/grch38_1kgmaj_snvs_bt2.zip)  |
| Bowtie 2 | hg19 + major SNPs             | [https](https://genome-idx.s3.amazonaws.com/bt/hg19_1kgmaj_snvs_bt2.zip)  |

#### SNP-and-indel

Pre-built major-allele **SNP-and-indel** reference indexes are available below.
The pre-built levioSAM index (`.lft`) is also included.
We provide a [tutorial](https://github.com/alshai/levioSAM/wiki/Alignment-with-variant-aware-reference-genomes) of using levioSAM in a major-allele alignment workflow.

| Aligner  | Reference                     | Index zip                                                            |
|----------|-------------------------------|----------------------------------------------------------------------|
| Bowtie 2 | GRCh38 + major SNP-and-indels | [https](https://genome-idx.s3.amazonaws.com/bt/grch38_1kgmaj_snvindels_bt2.zip)  |
| Bowtie 2 | hg19 + major SNP-and-indels   | [https](https://genome-idx.s3.amazonaws.com/bt/hg19_1kgmaj_snvindels_bt2.zip)  |

Those links **will** also appear on the [Bowtie web page](http://bowtie-bio.sourceforge.net) and [Bowtie 2 web page](http://bowtie-bio.sourceforge.net/bowtie2) in the right-hand sidebar.

The FASTA files with major-allele SNPs inserted are also available:

| Reference                     | FASTA file                                                        | LevioSAM index |
|-------------------------------|-------------------------------------------------------------------|----------------|
| GRCh38 + major SNPs           | [https](https://genome-idx.s3.amazonaws.com/bt/grch38_1kgmaj_snvs_bt2.fa.gz) | N/A |
| GRCh38 + major SNP-and-indels | [https](https://genome-idx.s3.amazonaws.com/bt/grch38_1kgmaj_snvindels_bt2.fa.gz) | [https](https://genome-idx.s3.amazonaws.com/bt/grch38_1kgmaj_snvindels.lft) |
| hg19 + major SNPs             | [https](https://genome-idx.s3.amazonaws.com/bt/hg19_1kgmaj_snvs_bt2.fa.gz) | N/A |
| hg19 + major SNP-and-indels   | [https](https://genome-idx.s3.amazonaws.com/bt/hg19_1kgmaj_snvindels_bt2.fa.gz) | [https](https://genome-idx.s3.amazonaws.com/bt/hg19_1kgmaj_snvindels.lft) |

### Instructions

Workflow:

1. `cd` to appropriate subdir
2. `./buildXX.sh` to build index
    * Use `sbatch buildXX_marcc.sh` if on JHU MARCC cluster
3. `./testXX.sh` to test
    * You may see a limited number of warnings, usually due to VCF formatting issues
    * Use `sbatch testXX_marcc.sh` if on JHU MARCC cluster
4. `./index_bt_marcc.sh` and `./index_bt2_marcc.sh` to build Bowtie and Bowtie 2 indexes
    * Use `sbatch` if on JHU MARCC cluster
5. `./zip_bt.sh` and `./zip_bt2.sh` to zip indexes
    * Archives all include `README.md`
6. `./scp_bt.sh` and `./scp_bt2.sh` to copy over to FTP server
    * Assumes you're on MARCC or other JHU cluster with access to `gwln1`

Requirements for building major-allele FASTAs (first 3 steps above):

* [bcftools](https://samtools.github.io/bcftools/) ([conda](https://anaconda.org/bioconda/bcftools))
* [samtools](https://samtools.github.io) ([conda](https://anaconda.org/bioconda/samtools))
* [tabix](http://www.htslib.org/doc/tabix.html) ([conda](https://anaconda.org/bioconda/tabix))

Requirements for building genome indexes (4th step above):

* [bowtie](http://bowtie-bio.sourceforge.net) ([conda](https://anaconda.org/bioconda/bowtie))
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2) ([conda](https://anaconda.org/bioconda/bowtie2))

Requirement for building major-allele FASTAs including indels (2th step above):
* [levioSAM](https://github.com/alshai/levioSAM) ([conda](https://anaconda.org/bioconda/leviosam))

### Authors

* Nae-Chyun Chen
* Ben Langmead

