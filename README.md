This collection of scripts are used to build [1KGenomes project](http://www.internationalgenome.org)-informed major-allele reference genomes for Bowtie 1 and Bowtie 2.  The current scripts build *linear* references with major-allele *SNPs only*.  Since no indels are included, the coordinate system will remain the same even once major alleles are inserted.  Thus, these indexes should be completely compatible with downstream tools as long as those tools are also using the major-allele references.  We are still considering the best way to also handle major-allele indels, etc.

#### Downloads

Pre-built major-allele-SNP reference indexes are available:

| Aligner  | Reference           | Index zip                                                            |
|----------|---------------------|----------------------------------------------------------------------|
| Bowtie   | GRCh38 + major SNPs | [ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/grch38_1kgmaj_bt.zip](ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/grch38_1kgmaj_bt.zip) |
| Bowtie   | hg19 + major SNPs   | [ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19_1kgmaj_bt.zip](ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19_1kgmaj_bt.zip)     |
| Bowtie 2 | GRCh38 + major SNPs | [ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/grch38_1kgmaj_bt2.zip](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/grch38_1kgmaj_bt2.zip)   |
| Bowtie 2 | hg19 + major SNPs   | [ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19_1kgmaj_bt2.zip](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19_1kgmaj_bt2.zip) |

Those links also appear on the [Bowtie web page](http://bowtie-bio.sourceforge.net) and [Bowtie 2 web page](http://bowtie-bio.sourceforge.net/bowtie2) in the right-hand sidebar.

The FASTA files with major-allele SNPs inserted are also available:

| Reference           | FASTA file                                                        |
|---------------------|-------------------------------------------------------------------|
| GRCh38 + major SNPs | [ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/grch38_1kgmaj.fa.gz](ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/grch38_1kgmaj.fa.gz) |
| hg19 + major SNPs   | [ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19_1kgmaj.fa.gz](ftp://ftp.ccb.jhu.edu/pub/data/bowtie_indexes/hg19_1kgmaj.fa.gz)   |

#### Instructions

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
* [vcftools](https://vcftools.github.io/perl_module.html) ([conda](https://anaconda.org/bioconda/vcftools))
* [samtools](https://samtools.github.io) ([conda](https://anaconda.org/bioconda/samtools))
* [tabix](http://www.htslib.org/doc/tabix.html) ([conda](https://anaconda.org/bioconda/tabix))

Requirements for building genome indexes (4th step above):

* [bowtie](http://bowtie-bio.sourceforge.net) ([conda](https://anaconda.org/bioconda/bowtie))
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2) ([conda](https://anaconda.org/bioconda/bowtie2))

#### Authors

* Nae-Chyun Chen
* Ben Langmead

