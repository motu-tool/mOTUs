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
> _Nuclic Acids Research_ (2025)
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
    genomic context to taxonomic profiling of microbial communities. Nuclic Acids Research (2025).
    doi: https://doi.org/10.1093/nar/gkae1004


    motus <command> [options]

        -- Taxonomic profiling
              profile     Perform taxonomic profiling (map_tax + calc_mgc + calc_motu) in a single step

              map_tax     Map reads to the marker gene database
              calc_mgc    Calculate marker gene cluster (MGC) abundance
              calc_motu   Summarize MGC abundances into a mOTU profile


        -- Tool utilities
              downloadMGDB  Download the mOTUs marker gene database
              merge         Merge multiple taxonomic profiling results into one table
              classify      Classify user genomes into mOTUs
              prep_long     Prepare long reads to be profiled by mOTUs


       -- Genome accession
              genomes     Search the mOTUs-db by keyword (taxonomic, functional)
              download    Download sequence files from mOTUs-db

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
       genomic context to taxonomic profiling of microbial communities. Nuclic Acids Research (2025).
       doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
       The profile function in mOTUs is the main function that executes map_tax, calc_mgc, and calc_motu in sequence.
       It takes short read metagenomic sequencing data as input and generates a taxonomic profile.


    Usage:
       motus profile -f FILE [FILE ...] -r FILE [FILE ...] -s FILE [FILE ...] -o FILE [options]
       motus profile -f FILE [FILE ...] -r FILE [FILE ...] -o FILE [options]
       motus profile -s FILE [FILE ...] -o FILE [options]


    Input options:
       -f, --forward  FILE[ FILE]
           Input file(s) for reads in forward orientation, fastq(.gz)-formatted

       -r, --reverse  FILE[ FILE]
           Input file(s) for reads in reverse orientation, fastq(.gz)-formatted

       -s, --single  FILE[ FILE]
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

### Merge

```bash
$ motus merge
```

<details>
<summary>merge cli options</summary>

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


    References:
        Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
        taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
        doi: https://doi.org/10.1186/s40168-022-01410-z

        Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
        genomic context to taxonomic profiling of microbial communities. Nuclic Acids Research (2025).
        doi: https://doi.org/10.1093/nar/gkae1004


    Summary:
        The merge function in mOTUs takes multiple profiles produced by the profile function and
        combines them into a single table.


    Usage:
        motus merge -i FILE [FILE ...] -o FILE


    Input options:
        -i, --input-files  FILE [FILE ...]
            A list of mOTUs profile files or a text file containing the list of profiles to be merged
            with one line per mOTUs profile file [required].

    Output options:
        -o, --output-file  FILE
            Output file name [required].


motus.py: error: the following arguments are required: -i, -o
```
</details>



#### Required arguments


Input: `-i`: Specifies the mOTUs profile files to merge. These files must be generated from the same mOTUs version and with the same parameters. At least two profiles are required. The input can be provided as a text file with one line per profile or as a space-separated list containing multiple mOTUs profiles.

Output: `-o`: Specifies the path to the merged profile file.

---


### downloadDB

```bash
$ motus downloadDB
```

<details>
<summary>downloadDB cli options</summary>

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


    References:

    Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
    taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
    doi: https://doi.org/10.1186/s40168-022-01410-z

    Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
    genomic context to taxonomic profiling of microbial communities. Nuclic Acids Research (2025).
    doi: https://doi.org/10.1093/nar/gkae1004


     motus downloadDB [options]

         Options:

            -f, --force        Force download even when database is already present


```
</details>



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

<details>
<summary>classify cli options</summary>

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


    References:

    Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
    taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
    doi: https://doi.org/10.1186/s40168-022-01410-z

    Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
    genomic context to taxonomic profiling of microbial communities. Nuclic Acids Research (2025).
    doi: https://doi.org/10.1093/nar/gkae1004


     motus classify [options]

         Options:

            -i        Text file with fasta formatted (gzip allowed) genome
                      files which will be associated with existing mOTUs
            -o        Output file. One line per genome with associated mOTU.
            -t        Number of threads (default = 1)


motus.py: error: the following arguments are required: -i, -o
```
</details>



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

<details>
<summary>prep_long cli options</summary>

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2


    References:

    Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand
    taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022).
    doi: https://doi.org/10.1186/s40168-022-01410-z

    Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible
    genomic context to taxonomic profiling of microbial communities. Nuclic Acids Research (2025).
    doi: https://doi.org/10.1093/nar/gkae1004


     motus prep_long [options]

        Input options:
           -i  FILE   long read file to convert, can be fasta(.gz) or fastq(.gz)
        Output options:
           -o  FILE   converted file, ready to be used by motus profile
        Algorithm options:
           -sl INT    splitting length for the long reads. (default = 300)
           -ml INT    minimum read length, shorter are discarded. (default = 50)

           
      motus.py: error: the following arguments are required: -i, -o
```
</details>



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

<details>
<summary>download cli options</summary>

```bash
python motus.py download
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.2
    
    
    References:
    
    Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
    taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
    doi: https://doi.org/10.1186/s40168-022-01410-z

    Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible 
    genomic context to taxonomic profiling of microbial communities. Nuclic Acids Research (2025). 
    doi: https://doi.org/10.1093/nar/gkae1004
    

     motus download [options]

         Input options:
            -i  FILE/STR Can be either a list of genome names (1-n) or 
                            a text file with genomes to download. One line 
                            per genome name. The input file is
                            compatible with the output of motus find.

         Output options:
            -o  PATH     Output folder.

         Options:
            -r, --representatives           Download only representative genomes.

           
motus.py: error: the following arguments are required: -o, -i
```
</details>



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
$ motus find
```

<details>
<summary>find cli options</summary>

```bash
Program: motus - a tool for marker gene-based OTU (mOTU) profiling
    Version: 4.0.3

    
    References:
    
    Profiler: Ruscheweyh, Milanese et al. Cultivation-independent genomes greatly expand 
    taxonomic-profiling capabilities of mOTUs across various environments. Microbiome (2022). 
    doi: https://doi.org/10.1186/s40168-022-01410-z

    Database: Dmitrijeva, Ruscheweyh et al. The mOTUs online database provides web-accessible 
    genomic context to taxonomic profiling of microbial communities. Nuclic Acids Research (2025). 
    doi: https://doi.org/10.1093/nar/gkae1004
    

     motus genomes [options]

        Input options:
            -i  FILE/STR Can be either a list of search queries (1-n) or
                            a text file with queries. One line
                            per query name. Queries can be genome names,
                            PFAM, KEGG or EGGNOG ids or GTDB taxonomy
                            names. Will offer suggestions if queries dont
                            match database entries exactly.
        Output options:
            -o  FILE     Genome names with or without annotations that were
                            found to match search queries.

            -d, --details  STR,[STR] Annotation to report. Choose any combination of
                            [KEGG, PFAM, EGGNOG, TAXONOMY], e.g.
                            -r KEGG,PFAM
                            

           
motus.py: error: the following arguments are required: -i, -o
```
    

</details>



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
