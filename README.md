<p align="center">
  <img src="https://raw.githubusercontent.com/motu-tool/mOTUs/master/pics/motu_logo.png" alt="mOTUs logo">
</p>

[![Actions](https://img.shields.io/github/actions/workflow/status/motu-tool/mOTUs/python-app.yml?branch=mOTUs4.1&logo=github&style=flat-square&maxAge=300)](https://github.com/motu-tool/mOTUs/actions)
[![PyPI](https://img.shields.io/pypi/v/motus-tool.svg?logo=pypi&style=flat-square&maxAge=3600)](https://pypi.org/project/motus-tool)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/motus?logo=anaconda&style=flat-square&maxAge=3600)](https://anaconda.org/bioconda/motus)
[![Python Versions](https://img.shields.io/pypi/pyversions/motus-tool.svg?logo=python&style=flat-square&maxAge=3600)](https://pypi.org/project/motus-tool/#files)
[![License](https://img.shields.io/badge/license-GPL--3.0-blue.svg?style=flat-square)](https://choosealicense.com/licenses/gpl-3.0/)
[![GitHub issues](https://img.shields.io/github/issues/motu-tool/mOTUs.svg?style=flat-square&maxAge=600)](https://github.com/motu-tool/mOTUs/issues)
[![Docs](https://img.shields.io/badge/docs-motus--tool.org-informational?style=flat-square)](https://www.motus-tool.org/)
[![Database](https://img.shields.io/badge/database-motus--db.org-orange?style=flat-square)](https://www.motus-db.org/)
[![Downloads](https://img.shields.io/conda/dn/bioconda/motus?style=flat-square&label=bioconda%20downloads&color=303f9f)](https://anaconda.org/bioconda/motus)
[![Paper 2022](https://img.shields.io/badge/paper-Microbiome%202022-teal.svg?style=flat-square)](https://doi.org/10.1186/s40168-022-01410-z)
[![Paper 2024](https://img.shields.io/badge/paper-NAR%202024-teal.svg?style=flat-square)](https://doi.org/10.1093/nar/gkae1004)

---

# mOTUs profiler


The mOTU profiler is a computational tool that estimates taxonomic abundance of known and currently unknown microbial community members using metagenomic shotgun sequencing data.


The current version of the mOTUs profiler is built on top of the genomic mOTUs database ([motus-db](https://motus-db.org/)) which is constructed from 919K isolate and single cell-amplified (SAGs) genomes and 2.83M metagenome-assembled genomes (MAGs) generated from over 117K metagenomic samples spanning diverse microbiomes, which include (in addition to the human and ocean microbiome) soil, freshwater and gastrointestinal tract microbiomes of ruminants and other animals, environments we found to be greatly underrepresented by reference genomes.  

In the current version, 124,295 species-level taxonomic units (mOTUs) were constructed using sequences of 10 single-copy marker genes recovered from these genomes. 30,256 mOTUs are represented by an isolate genome, whereas 94,039 mOTUs are represented by MAGs only.


Please cite the paper(s) corresponding to the version(s) you use:

| Version | Journal | Year | DOI | Citations |
|---------|---------|------|-----|-----------|
| mOTUs v1 | _Nature Methods_ | 2013 | [10.1038/nmeth.2693](https://doi.org/10.1038/nmeth.2693) | [![](https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fbadge.dimensions.ai%2Fdetails%2Fdoi%2F10.1038%2Fnmeth.2693&search=%3Cdiv%20class%3D%22count%22%3E(%5Cd*)%3C%2Fdiv%3E&replace=%241&style=flat-square&label=cited&cacheSeconds=3600)](https://badge.dimensions.ai/details/doi/10.1038/nmeth.2693) |
| mOTUs v2 | _Nature Communications_ | 2019 | [10.1038/s41467-019-08844-4](https://doi.org/10.1038/s41467-019-08844-4) | [![](https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fbadge.dimensions.ai%2Fdetails%2Fdoi%2F10.1038%2Fs41467-019-08844-4&search=%3Cdiv%20class%3D%22count%22%3E(%5Cd*)%3C%2Fdiv%3E&replace=%241&style=flat-square&label=cited&cacheSeconds=3600)](https://badge.dimensions.ai/details/doi/10.1038/s41467-019-08844-4) |
| mOTUs v3 (profiler) | _Microbiome_ | 2022 | [10.1186/s40168-022-01410-z](https://doi.org/10.1186/s40168-022-01410-z) | [![](https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fbadge.dimensions.ai%2Fdetails%2Fdoi%2F10.1186%2Fs40168-022-01410-z&search=%3Cdiv%20class%3D%22count%22%3E(%5Cd*)%3C%2Fdiv%3E&replace=%241&style=flat-square&label=cited&cacheSeconds=3600)](https://badge.dimensions.ai/details/doi/10.1186/s40168-022-01410-z) |
| mOTUs v3 (protocol) | _Current Protocols_ | 2022 | [10.1002/cpz1.218](https://doi.org/10.1002/cpz1.218) | [![](https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fbadge.dimensions.ai%2Fdetails%2Fdoi%2F10.1002%2Fcpz1.218&search=%3Cdiv%20class%3D%22count%22%3E(%5Cd*)%3C%2Fdiv%3E&replace=%241&style=flat-square&label=cited&cacheSeconds=3600)](https://badge.dimensions.ai/details/doi/10.1002/cpz1.218) |
| mOTUs database (v4) | _Nucleic Acids Research_ | 2025 | [10.1093/nar/gkae1004](https://doi.org/10.1093/nar/gkae1004) | [![](https://img.shields.io/badge/dynamic/regex?url=https%3A%2F%2Fbadge.dimensions.ai%2Fdetails%2Fdoi%2F10.1093%2Fnar%2Fgkae1004&search=%3Cdiv%20class%3D%22count%22%3E(%5Cd*)%3C%2Fdiv%3E&replace=%241&style=flat-square&label=cited&cacheSeconds=3600)](https://badge.dimensions.ai/details/doi/10.1093/nar/gkae1004) |




---

## 📦 Installation

The mOTUs profiler, written in Python 3 (>=3.12), can be executed on a 64-bit Linux or MacOS system. However, there are external dependencies that need to be pre-installed. These dependencies can be manually installed or, more conveniently, using the conda package manager.


### Installation with Conda


<details>
<summary>Miniconda</summary>

The installation using the conda package manager is generally preferable, as it encapsulates the entire installation process into a single command once conda is installed. Execute the following command to install conda:

```bash
$ curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
$ sh Miniconda3-latest-Linux-x86_64.sh
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```


If working on a MacOS system, the download link has to be replaced by: `https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh`.

</details>


mOTUs is available as a package in [bioconda](https://bioconda.github.io/recipes/motus/README.html) and can be installed in an isolated environment:

```bash
$ conda create -n mOTUs4 motus
$ conda activate mOTUs4
```



---

## 🚀 Usage



After installation, you can test whether the tool was installed correctly by executing:


```bash
$ motus --help
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.1.0


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Usage:
        motus <command> [options]


    Commands:

        -- Taxonomic profiling

            profile       Perform taxonomic profiling (map_tax + calc_mgc + calc_motu) in a single step

            map_tax       Map reads to the marker gene database
            calc_mgc      Calculate marker gene cluster (MGC) abundance
            calc_motu     Summarize MGC abundances into a mOTU profile


        -- Tool utilities

            downloadMGDB  Download the mOTUs marker gene database
            merge         Merge multiple taxonomic profiling results into one table
            classify      Classify user genomes into mOTUs
            prep_long     Prepare long reads to be profiled by mOTUs


        -- Genome accession

            genomes       Search the mOTUs-db by keyword (taxonomic, functional)
            download      Download sequence files from mOTUs-db


    Type motus <command> to print the help menu for a specific command

```

---

### Commands

The `profile` function in mOTUs is the main function that executes `map_tax`, `calc_mgc`, and `calc_motu` in sequence. It takes short read metagenomic sequencing data as input and generates a taxonomic profile.

Helper functions include `download`, which provides users with programmatic access to the ~4 million genomes in the motus-db; `downloadMGDB`, which downloads the marker gene database of mOTUs; `merge`, which merges multiple taxonomic profiles; and `classify`, which assigns user-submitted genomes to existing mOTUs.

---

### Profile

```bash
$ motus profile
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.1.0


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The profile command in mOTUs is the main function that executes map_tax, calc_mgc,
        and calc_motu in sequence. It takes short read metagenomic sequencing data as input
        and generates a taxonomic profile.


    Usage:
       motus profile -f FILE [FILE ...] -r FILE [FILE ...] -s FILE [FILE ...] -o FILE [options]
       motus profile -f FILE [FILE ...] -r FILE [FILE ...] -o FILE [options]
       motus profile -s FILE [FILE ...] -o FILE [options]


    Input options:
        -f, --forward  FILE [FILE ...]
            Input file(s) for reads in forward orientation, fastQ/A(.gz)-formatted

        -r, --reverse  FILE [FILE ...]
            Input file(s) for reads in reverse orientation, fastQ/A(.gz)-formatted

        -s, --single  FILE [FILE ...]
            Input file(s) for unpaired reads, fastQ/A(.gz)-formatted

        -n, --sample-name  STR
            Sample name (default: 'unnamed sample')

    Output options:
        -o, --output-file  FILE
            Output file name [required]

    Algorithm options:
        -g, --marker-genes  INT
            Required number of marker genes for a mOTU to be called present:
            1=higher recall, 6=higher precision, 10=maximum (default: 3)

        -l, --alignment-length  INT
            Minimum length of the alignment (bp) (default: 75)

        -t, --threads  INT
            Number of threads (default: 1)

        -y, --counting-mode  STR
            Which scale the abundances are reported in (default: INSERT_SCALED)
            Choices: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]

        --skip-pair-check
            Skip validation that forward and reverse read headers match.
            Use when reads are unsorted or contain singletons.

        -db  PATH
            Alternative path for the mOTUs marker gene database

```

---

### Map Tax

```bash
$ motus map_tax
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.1.0


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The map_tax command takes short read metagenomic sequencing data as input and
        maps reads to the mOTUs marker gene database.


    Usage:
        motus map_tax -f FILE [FILE ...] -r FILE [FILE ...] -s FILE [FILE ...] -o FILE [options]
        motus map_tax -f FILE [FILE ...] -r FILE [FILE ...] -o FILE [options]
        motus map_tax -s FILE [FILE ...] -o FILE [options]


    Input options:
        -f, --forward  FILE [FILE ...]
            Input file(s) for reads in forward orientation, fastQ/A(.gz)-formatted

        -r, --reverse  FILE [FILE ...]
            Input file(s) for reads in reverse orientation, fastQ/A(.gz)-formatted

        -s, --single  FILE [FILE ...]
            Input file(s) for unpaired reads, fastQ/A(.gz)-formatted

    Output options:
        -o, --output-file  FILE
            Output file name [required]

    Algorithm options:
        -l, --alignment-length  INT
            Minimum length of the alignment (bp) (default: 75)

        -t, --threads  INT
            Number of threads (default: 1)

        --skip-pair-check
            Skip validation that forward and reverse read headers match.
            Use when reads are unsorted or contain singletons.

        -db  PATH
            Alternative path for the mOTUs marker gene database

```


---

### Calc MGC

```bash
$ motus calc_mgc
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.1.0


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The calc_mgc command takes a file storing the alignments of sequencing reads
        to the mOTUs marker gene database and calculates marker gene cluster abundances.


    Usage:
        motus calc_mgc -i FILE -o FILE [options]


    Input options:
        -i, --input-file  FILE
            Path to BAM file generated after running the motus map_tax command [required]

    Output options:
        -o, --output-file  FILE
            Output file name [required]

    Algorithm options:
        -l, --alignment-length  INT
            Minimum length of the alignment (bp) (default: 75)

        -db  PATH
            Alternative path for the mOTUs marker gene database

```
---

### Calc mOTU

```bash
$ motus calc_motu
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.1.0


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The calc_motu command takes a file containing marker gene cluster
        abundances and generates a taxonomic profile.


    Usage:
        motus calc_motu -i FILE -o FILE [options]


    Input options:
        -i, --input-file  FILE
            MGC abundance table generated by the calc_mgc command [required]

        -n, --sample-name  STR
            Sample name (default: 'unnamed sample')

    Output options:
        -o, --output-file  FILE
            Output file name [required]

    Algorithm options:
        -g, --marker-genes  INT
            Required number of marker genes for a mOTU to be called present:
            1=higher recall, 6=higher precision, 10=maximum (default: 3)

        -y, --counting-mode  STR
            Which scale the abundances are reported in (default: INSERT_SCALED)
            Choices: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]

        -db  PATH
            Alternative path for the mOTUs marker gene database (default: built-in location)

```

---

### merge

```bash
$ motus merge
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.1.0


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The merge command takes multiple profiles produced after running the
        profile command and combines them into a single table.


    Usage:
        motus merge -i FILE [FILE ...] -o FILE


    Input options:
        -i, --input-files  FILE [FILE ...]
            A list of mOTUs profile files or a text file containing the list of profile
            files to be merged, with one line per file [required]

    Output options:
        -o, --output-file  FILE
            Output file name [required]

    Database options:
        -db  PATH
            Alternative path for the mOTUs marker gene database

```


---


### downloadMGDB

```bash
$ motus downloadMGDB
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.1.0


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The downloadMGDB command downloads the marker gene reference database used
        by the profile and map_tax commands.


    Usage:
        motus downloadMGDB [options]


    Options:
        -f, --force
            Force download even when database is already present

        --toy
            Download the lightweight toy database (v4.1-toy) instead of the full database.
            Useful for testing and development. Note: the genomes and download commands
            are not available with the toy database.

        -db  PATH
            Alternative path for the mOTUs marker gene database

```



---

### classify

```bash
$ motus classify
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.1.0


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The classify command takes a list of genome sequence files as input and
        assigns these genomes to existing mOTUs in the database.
        Requires vsearch to be installed and on PATH.


    Usage:
        motus classify -i FILE -o FILE [options]


    Input options:
        -i, --input-file  FILE
            Text file listing genome sequence files in fastA(.gz) format to classify.
            One line per genome file [required]

    Output options:
        -o, --output-file  FILE
            Output file name [required]

    Algorithm options:
        -t, --threads  INT
            Number of threads (default: 1)

        -db  PATH
            Alternative path for the mOTUs marker gene database

```

The output file is a tab-separated table with one row per input genome:

| Column | Description |
|---|---|
| `GENOME` | Input genome filename |
| `CLOSEST_MOTU` | Best-matching mOTU by combined marker gene similarity. `no_mOTU` if no hit was found with ≥6 marker genes; `no_mOTU_<6MGs` if fewer than 6 marker genes were extracted. |
| `SIMILARITY` | Combined percent identity to the closest mOTU (0–100). `-1.0` if no mOTU could be assigned. |
| `ASSIGNED_TO_MOTU` | `True` if similarity ≥96.5% (genome is within the mOTU boundary), `False` otherwise. |
| `TAXONOMY` | GTDB taxonomy of the closest mOTU. `d__;p__;c__;o__;f__;g__;s__` if no mOTU was assigned. |
| `#MGs` | Number of marker genes extracted from the genome by fetchMGs. |


---



### prep_long

```bash
$ motus prep_long
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.1.0


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The prep_long command takes long-read sequencing data and converts it
        into the appropriate input format to be used by the profile and map_tax commands.


    Usage:
        motus prep_long -i FILE -o FILE [options]


    Input options:
        -i, --input-file  FILE
            Long-read sequencing file to convert, can be in fastQ/A(.gz) format [required]

    Output options:
        -o, --output-file  FILE
            Output file name. This converted file is ready to be used by motus profile [required]

    Algorithm options:
        -sl, --splitting-length  INT
            Target fragment length (in bp) for splitting long reads (default: 300)

        -ml, --minimum-length  INT
            Minimum read length after splitting. Shorter reads are discarded (default: 50)

           
```


---



### download

```bash
$ motus download
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.1.0


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The download command downloads listed genome files from mOTUs-db.


    Usage:
        motus download -i FILE -o PATH [options]
        motus download -i STR [STR ...] -o PATH [options]


    Input options:
        -i, --input-genomes  FILE/STR
            Can be either a list of genome identifiers separated by spaces or a text file
            listing the identifiers of genomes for download. One line per genome. The output of
            the motus genomes command can be used as input for this command [required]

    Output options:
        -o, --output-folder  PATH
            Path to output folder where the downloaded sequences will be saved [required]

        -r, --representatives
            Download only sequences from representative genomes.

    Algorithm options:
        -t, --file-type  STR
            File type to download (default: genome)
            Choices: [genome, gene_fna, gene_faa, gene_gff, antismash, pfam, eggnog, kegg, trna, rrna]

    Database options:
        -db  PATH
            Alternative path for the mOTUs marker gene database

```



---




### genomes

```bash
$ motus genomes
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.1.0


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The genomes command queries the mOTUs-db based on identifiers, functional,
        or taxonomic annotations and returns a list of genomes matching indicated query.


    Usage:
        motus genomes -i FILE -o FILE [options]
        motus genomes -i STR [STR ...] -o FILE [options]
        motus genomes -l GENOME|TAXON|PFAM|KEGG|EGGNOG -o FILE [options]


    Input options:
        -i, --input-queries  FILE/STR
            Can be either a list of search queries or a text file listing search queries
            with one line per query. Queries can be genome or mOTUs identifiers, PFAM, KEGG, EGGNOG,
            or GTDB taxonomy names. If the query does not exactly match any database entry,
            alternative queries will be suggested [required unless -l is used]

        -l, --list  STR
            List all searchable entries for a given category and write them to -o.
            Choose from [GENOME, TAXON, PFAM, KEGG, EGGNOG]. When used, -i is not required.

    Output options:
        -o, --output-file  FILE
            Output file containing a list of genome identifiers matching search queries and their
            annotations as indicated by the -d parameter. This output file can be used as input
            for the motus download command [required]

        -d, --details  STR [STR ...]
            List of annotations to report. Choose any combination of [KEGG, PFAM, EGGNOG, TAXONOMY],
            for example, -d KEGG PFAM.

        -db  PATH
            Alternative path for the mOTUs database

```
    

---



## ❓ Need Help?

Write an issue on GitHub

---

## 📋 Changelog

### v4.1.0

**Database**
- Default marker gene database updated to v4.1
- Annotation database (used by `genomes`) is now version-matched to the installed marker gene DB; v4.0 and v4.1 annotation DBs are downloaded and stored separately
- GTDB taxonomy files parsed by column name rather than position to handle format differences between v4.0 (`GTDBR220`) and v4.1 (`GTDB`)
- Toy database added (`downloadMGDB --toy`): lightweight database for testing; disables `genomes` and `download` commands

**classify**
- Output columns changed to `GENOME`, `CLOSEST_MOTU`, `SIMILARITY`, `ASSIGNED_TO_MOTU`, `TAXONOMY`, `#MGs`
- Reports one best mOTU per genome with GTDB taxonomy; `ASSIGNED_TO_MOTU` is `True` if similarity ≥96.5%
- Unclassified genomes reported as `no_mOTU` (hits found but all below threshold) or `no_mOTU_<6MGs` (fewer than 6 marker genes extracted)

**CLI**
- `-db PATH` flag added to all commands to specify a custom database parent folder
- `--skip-pair-check` added to `profile` and `map_tax` for unsorted inputs or inputs containing singletons
- `motus genomes -l GENOME|TAXON|PFAM|KEGG|EGGNOG` lists all searchable entries of a given type without requiring `-i`
- Tool checks at startup whether `bwa` is on PATH (hard error if missing) and whether `vsearch` is on PATH (warning if missing)

**Bug fixes**
- Database download: incomplete downloads (content-length mismatch) no longer leave a broken DB with a valid completion marker
- `merge`: blank lines in a profile list file no longer cause a cryptic `FileNotFoundError`
- `MergedmOTUsFile`: class-level `_singlemotusfiles` dict replaced with an instance variable, preventing state leakage between calls in the same process
- Multimapper resolution: replaced non-deterministic `random.choice` with `sorted()[0]` for reproducible MG selection within an MGC
- Edge correction: fixed inaccurate weight distribution in inverse padding

**Algorithm**
- `INSERT_NORM` and `INSERT_SCALED` calculation corrected: length-normalisation was incorrectly folded into the scaling denominator

---

### v4.0.x

- Initial public release of mOTUs4 with database v4.0 (124,295 mOTUs)
- Three-stage pipeline: `map_tax` → `calc_mgc` → `calc_motu`; `profile` runs all three in sequence
- `classify` command: genome-to-mOTU assignment via fetchMGs + vsearch
- `genomes` and `download` commands for programmatic access to mOTUs-db genome sequences
- `merge` command for combining multiple single-sample profiles into one table
- `prep_long` command: splits long reads into ~300 bp fragments for profiling
- Read name normalisation: `/1`/`/2` suffixes stripped and re-appended as BAM qname tags to track orientation through the pipeline