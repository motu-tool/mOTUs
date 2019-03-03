![alt text](motu_logo.png)

[![Build status](https://ci.appveyor.com/api/projects/status/0x4veuuoabm6018v/branch/master?svg=true)](https://ci.appveyor.com/project/AlessioMilanese/motus-v2/branch/master)
[![Install with Bioconda](https://anaconda.org/bioconda/motus/badges/installer/conda.svg)](https://anaconda.org/bioconda/motus)
[![container ready](https://quay.io/repository/biocontainers/motus/status)](https://quay.io/repository/biocontainers/motus)
[![license](https://anaconda.org/bioconda/motus/badges/license.svg)](https://github.com/motu-tool/mOTUs_v2/blob/master/LICENSE)
[![Install with Bioconda](https://img.shields.io/conda/dn/bioconda/motus.svg?style=flat)](https://anaconda.org/bioconda/motus)


mOTUs profiler
========

The mOTUs profiler is a computational tool that estimates relative abundance of known and currently unknown microbial community members using metagenomic shotgun sequencing data.

Check the [wiki](https://github.com/motu-tool/mOTUs_v2/wiki) for more information.

If you are using mOTUs2, please cite:

> Alessio Milanese, Daniel R Mende, Lucas Paoli, Guillem Salazar, Hans-Joachim Ruscheweyh, Miguelangel Cuenca,
> Pascal Hingamp, Renato Alves, Paul I Costea, Luis Pedro Coelho, Thomas S B Schmidt,
> Alexandre Almeida, Alex L Mitchell, Robert D Finn, Jaime Huerta-Cepas,
> Peer Bork, Georg Zeller, Shinichi Sunagawa.
> **Microbial abundance, activity and population genomic profiling with mOTUs2**; _Nature Communications_ (2019).
> doi: [10.1038/s41467-019-08844-4](https://doi.org/10.1038/s41467-019-08844-4)


Pre-requisites
--------------

The mOTU profiler requires:
* Python 3 (or higher)
* the Burrow-Wheeler Aligner v0.7.15 or higher ([bwa](https://github.com/lh3/bwa))
* SAMtools v1.5 or higher ([link](http://www.htslib.org/download/))

In order to use the command ```snv_call``` you need:
* [metaSNV v1.0.3](https://git.embl.de/costea/metaSNV), available also on [bioconda](https://anaconda.org/bioconda/metasnv) (we assume metaSNV.py to be in the system path)

Check [installation wiki](https://github.com/motu-tool/mOTUs_v2/wiki/Installation) to see how to install the dependencies with conda.

Installation
--------------
```bash
git clone https://github.com/motu-tool/mOTUs_v2.git
cd mOTUs_v2
python setup.py
python test.py
export PATH=`pwd`:$PATH
```

Note: in the following examples we assume that the python script ```motus``` is in the system path.


Simple examples
--------------
Here is a simple example on how to obtain a taxonomic profiling from a raw read file:

```bash
motus profile -s metagenomic_sample.fastq > taxonomy_profile.txt
```

You can separate the previous call as:
```bash
motus map_tax -s metagenomic_sample.fastq -o mapped_reads.sam
motus calc_mgc -i mapped_reads.sam -o mgc_ab_table.count
motus calc_motu -i mgc_ab_table.count > taxonomy_profile.txt
rm mapped_reads.sam mgc_ab_table.count
```


The use of multiple threads (`-t`) is recommended, since bwa will finish faster. Here is an example with Paired-End reads:

```bash
motus profile -f for_sample.fastq -r rev_sample.fastq -s no_pair.fastq -t 6 > taxonomy_profile.txt
```

You can merge taxonomy files from different samples with `mOTU merge`:

```shell
motus profile -s metagenomic_sample_1.fastq -o taxonomy_profile_1.txt
motus profile -s metagenomic_sample_2.fastq -o taxonomy_profile_2.txt
motus merge -i taxonomy_profile_1.txt,taxonomy_profile_2.txt > all_sample_profiles.txt
```

You can profile samples that have been sequenced through different runs:
```shell
motus profile -f sample1_run1_for.fastq,sample1_run2_for.fastq -r sample1_run1_rev.fastq,sample1_run2_rev.fastq -s sample1_run1_single.fastq > taxonomy_profile.txt
```

ChangeLog
--------------
**Version 2.1.0 2019-03-03 by AlessioMilanese**
* Correct error \'\t\t\' when printing -C recall
* Update bioconda environment .yaml file
* Update database (gene coordinates)

**Version 2.0.1 2018-08-23 by AlessioMilanese**
* Add -C to print the result in CAMI format (BioBoxes format 0.9.1)
* Add -K to snv_call command to keep all the directories produced by metaSNV

**Version 2.0.0 2018-06-12 by AlessioMilanese**
* Set relative abundances as default (instead of counts)
* Add -B to print the result in BIOM format
* Add test directory
* Python2 is not supported anymore
* Minor bug fixes

**Version 2.0.0-rc1 2018-05-10 by AlessioMilanese**
* First release supporting all basic functionality.
