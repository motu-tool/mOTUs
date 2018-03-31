![alt text](motu_logo.png)

mOTUs profiler
========

The mOTUs profiler is a computational tool that estimates relative abundance of known and currently unknown microbial community members using metagenomic shotgun sequencing data. The profiling method is based on ~100,000 universal single-copy marker genes that have been extracted from more than 3,000 metagenomic samples and ~25,000 reference genomes.

Check the [wiki](https://github.com/motu-tool/mOTUs_v2/wiki) for more information.

Pre-requisites
--------------

The mOTU profiler requires:
* Python 3 or Python 2.7 (or higher)
* the Burrow-Wheeler Aligner ([bwa](https://github.com/lh3/bwa))
* SAMtools ([link](http://samtools.sourceforge.net/))

In order to use the command ```snv_call``` you need:
* [metaSNV](http://metasnv.embl.de/), now available also on [bioconda](https://anaconda.org/bioconda/metasnv)

Installation
--------------
```bash
git clone https://github.com/motu-tool/mOTUs_v2.git
cd mOTUs_v2
python setup.py
export PATH=`pwd`:$PATH
```

Note: in the following examples we assume that the python script ```motus``` is in the system path.


Simple examples
--------------
Here is a simple example on how to obtain a taxonomic profiling from a raw read file:

```bash
motus profile -s metagenomic_sample.fastq > taxonomic_profiling.txt
```

You can separate the previous call as:
```bash
motus map_tax -s metagenomic_sample.fastq -o mapped_reads.sam
motus calc_mgc -i mapped_reads.sam -o mgc_ab_table.count
motus calc_motu -i mgc_ab_table.count > taxonomic_profiling.txt
rm mapped_reads.sam mgc_ab_table.count
```


The use of multiple threads (`-t`) is recommended, since bwa will finish faster. Here is an example with Paired-End reads:

```bash
motus profile -f for_sample.fastq -r rev_sample.fastq -s no_pair.fastq -t 6 > taxonomic_profiling.txt
```

You can merge taxonomy files from different samples with `mOTU merge`:

```shell
motus profile -s metagenomic_sample_1.fastq -o taxonomic_profiling_1.txt
motus profile -s metagenomic_sample_2.fastq -o taxonomic_profiling_2.txt
motus merge -i taxonomic_profiling_1.txt,taxonomic_profiling_2.txt > all_samples_profiling.txt
```

You can profile samples that have been sequenced through different runs:
```shell
motus profile -f sample1_run1_for.fastq,sample1_run2_for.fastq -r sample1_run1_rev.fastq,sample1_run2_rev.fastq -s sample1_run1_single.fastq > taxonomic_profiling.txt
```
