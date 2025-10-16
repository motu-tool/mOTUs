![alt text](https://raw.githubusercontent.com/motu-tool/mOTUs/master/pics/motu_logo.png)

[![license](https://anaconda.org/bioconda/motus/badges/license.svg)](https://github.com/motu-tool/mOTUs4-dev/blob/main/LICENSE)

---

# mOTUs profiler


The mOTU profiler is a computational tool that estimates taxonomic abundance of known and currently unknown microbial community members using metagenomic shotgun sequencing data.


The current version of the mOTUs profiler is built on top of the genomic mOTUs database ([motus-db](https://motus-db.org/)) which is constructed from 919K isolate and single cell-amplified (SAGs) genomes and 2.83M metagenome-assembled genomes (MAGs) generated from over 117K metagenomic samples spanning diverse microbiomes, which include (in addition to the human and ocean microbiome) soil, freshwater and gastrointestinal tract microbiomes of ruminants and other animals, environments we found to be greatly underrepresented by reference genomes.  

In the current version, 124,295 species-level taxonomic units (mOTUs) were constructed using sequences of 10 single-copy marker genes recovered from these genomes. 30,256 mOTUs are represented by an isolate genome, whereas 94,039 mOTUs are represented by MAGs only.


If you use the mOTUs profiler, please cite:

> **Reference genome-independent taxonomic profiling of microbiomes with mOTUs3**
> 
> Hans-Joachim Ruscheweyh* , Alessio Milanese*, Lucas Paoli, Nicolai Karcher, Quentin Clayssen,
> Marisa Isabell Metzger, Jakob Wirbel, Peer Bork, Daniel R. Mende, Georg Zeller# & Shinichi Sunagawa#
> 
> _Microbiome_ (2022)
> 
> doi: [10.1186/s40168-022-01410-z](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-022-01410-z)
 
If you use the mOTUs database, please cite:

>  **The mOTUs online database provides web-accessible genomic context to taxonomic profiling of microbial communities**
> 
> Marija Dmitrijeva* , Hans-Joachim Ruscheweyh* , Lilith Feer , Kang Li , Samuel Miravet-Verde , Anna Sintsova , Daniel R Mende , Georg Zeller , Shinichi Sunagawa#
> 
> _Nucleic Acids Research_ (2025)
> 
> doi: [https://doi.org/10.1093/nar/gkae1004](https://doi.org/10.1093/nar/gkae1004)

---

## 📦 Installation

The mOTUs profiler, written in Python 3, can be executed on a 64-bit Linux or MacOS system. However, there are external dependencies that need to be pre-installed. These dependencies can be manually installed or, more conveniently, using the conda package manager.


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


Install using conda:

```bash
$ git clone git@github.com:motu-tool/mOTUs4-dev.git
$ cd mOTUs4-dev
$ conda env create -f mOTUs4-dev-conda.yaml
```

This will install:

```bash
- python==3.12.0
- bwa==0.7.19
- fetchmgs==2.1.0
- pysam=0.23.3
- polars==1.32.2
- rapidfuzz==3.13.0
- biopython==1.85
```


---

## 🚀 Usage



After installation, you can test whether the tool was installed correctly by executing:


```bash
$ motus --help
```

**Note** Currently the command to execute mOTUs is `python motus/motus.py` which will be replaced with `motus` once the tool is installed via `pip`.

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


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

motus.py: error: the following arguments are required: command
```

### Commands

The `profile` function in mOTUs is the main function that executes `map_tax`, `calc_mgc`, and `calc_motu` in sequence. It takes short read metagenomic sequencing data as input and generates a taxonomic profile.

Helper functions include `download`, which provides users with programmatic access to the ~4 million genomes in the motus-db; `downloadDB`, which downloads the marker gene database of mOTUs; `merge`, which merges multiple taxonomic profiles; and `classify`, which assigns user-submitted genomes to existing mOTUs.

---

### Profile

```bash
$ motus profile
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


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
            Input file(s) for reads in forward orientation, fastq(.gz)-formatted

        -r, --reverse  FILE [FILE ...]
            Input file(s) for reads in reverse orientation, fastq(.gz)-formatted

        -s, --single  FILE [FILE ...]
            Input file(s) for unpaired reads, fastq(.gz)-formatted

        -n, --sample-name  STR
            Sample name (default: 'unnamed sample')

    Output options:
        -o, --output-file  FILE
            Output file name [required]

        -a, --relative-abundance
            Write a second output file with relative abundances (default: False)

    Algorithm options:
        -g, --marker-genes  INT
            Required number of marker genes for a mOTU to be called present: 1=higher recall, 6=higher precision, 10=maximum (default: 3)

        -l, --alignment-length  INT
            Minimum length of the alignment (bp) (default: 75)

        -t, --threads  INT
            Number of threads (default: 1)

        -y, --counting-mode  STR
            Which scale the abundances are reported in (default: INSERT_SCALED)
            Values: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]


motus.py: error: the following arguments are required: -o
```

### Map Tax

```bash
$ motus map_tax
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The map_tax command in mOTUs takes short read metagenomic sequencing data as input and
        maps reads to the mOTUs marker gene database.


    Usage:
       motus map_tax -f FILE [FILE ...] -r FILE [FILE ...] -s FILE [FILE ...] -o FILE [options]
       motus map_tax -f FILE [FILE ...] -r FILE [FILE ...] -o FILE [options]
       motus map_tax -s FILE [FILE ...] -o FILE [options]


    Input options:
        -f, --forward  FILE [FILE ...]
            Input file(s) for reads in forward orientation, fastq(.gz)-formatted

        -r, --reverse  FILE [FILE ...]
            Input file(s) for reads in reverse orientation, fastq(.gz)-formatted

        -s, --single  FILE [FILE ...]
            Input file(s) for unpaired reads, fastq(.gz)-formatted

    Output options:
        -o, --output-file  FILE
            Output file name [required]

    Algorithm options:
        -l, --alignment-length  INT
            Minimum length of the alignment (bp) (default: 75)

        -t, --threads  INT
            Number of threads (default: 1)

```

### Calc MGC

```bash
$ motus calc_mgc
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The calc_mgc command in mOTUs takes a file storing the alignments of sequencing reads
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

```

### Calc mOTU

```bash
$ motus calc_motu
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The calc_motu command in mOTUs takes a file containing marker gene cluster
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

        -a, --relative-abundance
            Write a second output file with relative abundances (default: False)

    Algorithm options:
        -g, --marker-genes  INT
            Required number of marker genes for a mOTU to be called present: 1=higher recall, 6=higher precision, 10=maximum (default: 3)

        -y, --counting-mode  STR
            Which scale the abundances are reported in (default: INSERT_SCALED)
            Values: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]
```

### Merge

```bash
$ motus merge
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The merge command in mOTUs takes multiple profiles produced after running the
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


motus.py: error: the following arguments are required: -i, -o
```

#### Required arguments


Input: `-i`: Specifies the mOTUs profile files to merge. These files must be generated from the same mOTUs version and with the same parameters. At least two profiles are required. The input can be provided as a text file with one line per profile or as a space-separated list containing multiple mOTUs profiles.

Output: `-o`: Specifies the path to the merged profile file.

---


### downloadMGDB

```bash
$ motus downloadMGDB
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The downloadMGDB command in mOTUs downloads the marker gene reference database used
        by the profile and map_tax commands.


    Usage:
        downloadMGDB [options]


    Options:
        -f, --force
            Force download even when database is already present

```




#### Required arguments


#### 🔧 Options

| Option               | Description                                                                                                                                                                                                                                                                |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-f`     | **Force**: Download the database even if it’s already present.


---

### classify

```bash
$ motus classify
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The classify command in mOTUs takes a list of genome sequence files as input and
        assigns these genomes to existing mOTUs in the database.


    Usage:
        motus classify -i FILE -o FILE [options]


    Input options:
        -i, --input-file
            Text file listing genome sequence files in .fasta(.gz) format to classify.
            One line per genome file [required]

    Output options:
        -o, --output-file
            Output file name. Each line contains a genome and its associated mOTU [required]

    Algorithm options:
        -t, --threads
            Number of threads (default: 1)


motus.py: error: the following arguments are required: -i, -o
```


#### Required arguments


Input: `-i`: A text file containing fasta-formatted (gzip allowed) genome files that will be associated with existing mOTUs.

Output:`-o`: The output file, containing one line per genome with its associated mOTU.
  
#### 🔧 Options

| Option               | Description                                                                                                                                                                                                                                                                |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-i`     | **Genomes**: A text file containing fasta-formatted (gzip allowed) genome files that will be associated with existing mOTUs.
| `-o`     | **Output file**: The output file, containing one line per genome with its associated mOTU.
| `-t`     | **Threads**: The number of threads to be used for the alignment step.



---



### prep_long

```bash
$ motus prep_long
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The prep_long command in mOTUs takes long-read sequencing data and converts it
        into the appropriate input format to be used by the profile and map_tax commands.


    Usage:
        motus prep_long -i FILE -o FILE [options]


    Input options:
        -i, --input-file  FILE
            Long-read sequencing file to convert, can be in fasta(.gz) or fastq(.gz) format [required]

    Output options:
        -o, --output-file  FILE
            Output file name. This converted file is ready to be used by motus profile [required]

    Algorithm options:
        -sl, --splitting-length  INT
            Target fragment length (in bp) for splitting long reads (default: 300)

        -ml, --minimum-length  INT
            Minimum read length after splitting. Shorter reads are discarded (default: 50)

           
      motus.py: error: the following arguments are required: -i, -o
```




#### Required arguments


Input: `-i`: The input file containing long reads. It can be in fastA(.gz) or fastQ(.gz) format.

Output: `-o`: The output file where the converted reads will be stored in fastA format.  

  
#### 🔧 Options

| Option               | Description                                                                                                                                                                                                                                                                |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-i`     | **Input file**: The input file containing long reads. It can be in fastA(.gz) or fastQ(.gz) format.
| `-o`     | **Output file**: The output file where the converted reads will be stored in fastA format.
| `-sl`     | **Split length**: The length of short reads. The default value is 300.
| `-ml`     | **Minimum length**: Reads shorter than this length will not be written to the output. The default value is 75.


---



### download

```bash
$ motus download
```

```bash
python motus.py download
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2
    
    
    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible 
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025). 
        doi: https://doi.org/10.1093/nar/gkae1004
    

    Summary:
        The download command in mOTUs downloads sequences of requested genomes from mOTUs-db.


    Usage:
        motus download -i FILE -o PATH [options]
        motus download -i STR [STR ...] -o PATH [options]


    Input options:
        -i, --input-genomes  FILE/STR
            Can be either a list of genome identifiers (1-n) separated by spaces or a text file
            listing the identifiers of genomes for download. One line per genome. The output of
            the motus genomes command can be used as input for this command.

    Output options:
        -o, --output-folder  PATH
            Path to output folder where the downloaded sequences will be saved.

        -r, --representatives
            Download only sequences from representative genomes.

           
motus.py: error: the following arguments are required: -o, -i
```

#### Required arguments


The `-i` parameter allows you to specify the names of the genomes to download. You can provide a single genome name, multiple genome names, or a file containing genome names, with one name per line. The output file of `motus find` is compatible with this parameter.


The `-o` parameter specifies the output folder where the genomes will be downloaded. If the folder doesn’t exist, it will be created. Files with the same names in this folder will be overwritten without warning.  
 

  
#### 🔧 Options

| Option               | Description                                                                                                                                                                                                                                                                |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-i`     | **Genomes**: This parameter allows you to specify the names of the genomes to download. You can provide a single genome name, multiple genome names, or a file containing genome names, with one name per line.
| `-o`     | **Output folder**: This parameter specifies the output folder where the genomes will be downloaded to.
| `-r`     | **representative only**: Download only representative genomes.


---



### find

```bash
$ motus genomes
```

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.3

    
    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible 
        genomic context to taxonomic profiling of microbial communities. Nucleic Acids Research (2025). 
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The genomes command in mOTUs queries mOTUs-db based on identifiers, functional,
        or taxonomic annotations and returns a list of genomes matching indicated query.


    Usage:    
        motus genomes -i FILE -o FILE [options]


    Input options:
        -i, --input-queries  FILE/STR
            Can be either a list of search queries (1-n) or a text file listing search queries
            with one line per query. Queries can be genome identifiers, PFAM, KEGG or EGGNOG ids
            or GTDB taxonomy names. If the query does not exactly match any database entry,
            alternative queries will be suggested.

    Output options:
        -o, --output-file  FILE
            Output file containing a list of genome identifiers matching search queries and their
            annotations as indicated by the -d parameter. This output file can be used as input
            for the motus download command.

        -d, --details  STR,[STR]
            List of annotations to report. Choose any combination of [KEGG, PFAM, EGGNOG, TAXONOMY],
            for example, -d KEGG,PFAM.
                            

motus.py: error: the following arguments are required: -i, -o
```
    

#### Required arguments


The `-i` parameter allows you to specify the queries used to search for genomes. You can provide a single query, multiple queries, or a file containing queries, with one name per line. Queries can be genome names, GTDB taxonomy or annotation identifiers such as KEGG, PFAM or EGGNOG. Queries that dont match the database exactly will be used to suggest alternatives using fuzzy search.

The `-o` parameter specifies the output file where genomes with or without annotations will be stored.  
 

  
#### 🔧 Options

| Option               | Description                                                                                                                                                                                                                                                                |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-i`     | **Search Queries**: allows you to specify the tokens used to search for genomes. You can provide a single token, multiple tokens, or a file containing tokens, with one name per line. 
| `-o`     | **Output file**: This parameter specifies the output file where genomes with or without annotations will be stored.
| `-r`     | **Report**: Decide on which annotations to report. Can be any combination of [KEGG, PFAM, EGGNOG, TAXONOMY], e.g. -r KEGG,PFAM



## ❓ Need Help?

Write a issue on GitHub
