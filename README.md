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

        -- Utilities
              download    Download genomes associated with mOTUs
              downloadDB  Download the mOTUs marker gene database
              merge       Merge multiple taxonomic profiling results into one table
              classify    Classify user genomes into mOTUs
              prep_long   Prepare long reads to be profiled by mOTUs


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

<details>
<summary>profile cli options</summary>

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


    motus profile [options]

    Input options:
       -f  FILE[ FILE]  input file(s) for reads in forward orientation, fastq(.gz)-formatted
       -r  FILE[ FILE]  input file(s) for reads in reverse orientation, fastq(.gz)-formatted
       -s  FILE[ FILE]  input file(s) for unpaired reads, fastq(.gz)-formatted
       -n  STR          sample name ['unnamed sample']

    Output options:
       -o  FILE         output file name [required]
       -c               Write second output file with relative abundances

    Algorithm options:
       -g  INT          number of marker genes cutoff: 1=higher recall, 6=higher precision, 10=maximum [3]
       -l  INT          min length of the alignment (bp) [75]
       -t  INT          number of threads [1]
       -y  STR          type of read counts [INSERT_SCALED]
                        Values: [INSERT_RAW, INSERT_NORM, INSERT_SCALED, BASE_RAW, BASE_NORM]

motus.py: error: the following arguments are required: -o
```
</details>

#### Required arguments


Input: `-f, -r, -s`: One or multiple fastq/fasta files, which can be gzipped. The order of input files matters if using paired-end data (`-f, -r`).

Output: `-o`: Path to the output file. This also serves as a prefix for intermediate files.






#### 🔧 Options

| Option               | Description                                                                                                                                                                                                                                                                |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-f,-r`     | **Input path - Paired**: One or more gzipped fasta/q files. The input files must have the same order in both -f and -r.
| `-s`     | **Input path - Single**: One or more gzipped fasta/q files. The order of the input files doesn’t matter for single-end files.
| `-o`     | **Output prefix**: Path to the output file. This prefix is also used for intermediate files.
| `-n`       | **Samplename**: Name of the sample. Required when merging samples. The default value is “unnamed sample”.
| `-c`       | **Relative Abundance**: Report relative abundance in addition to counts.  
| `-g`       | **Sensitivity**: The number of marker genes with abundance required to call a mOTU present. The default value is 3, with a minimum of 1 and a maximum of 10. A value of 1 results in high recall but low precision, while a value of 10 results in high precision but low recall.   
| `-l`       | **Length**: Filter alignments if their length is below this value. **Note**: Choose a value greater than or equal to the length of the reads. Default value is `75`.
| `-t`       | **Threads**: Number of threads to use for the alignment step. Default is 1.
| `-y`       | **Counting method**: mOTUs can count in different modes. For more details, see the Wiki. The default mode is INSERT\_SCALED. Other options include INSERT\_RAW, INSERT\_NORM, INSERT\_SCALED, BASE\_RAW, and BASE\_NORM.



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


    motus merge [options]

        Input options:
           -i  FILE[ FILE]  A list of mOTUs profile files or a text file with one line
                            per mOTUs profile files to be merged


        Output options:
           -o  FILE  output file name


motus.py: error: the following arguments are required: -i, -o
```
</details>



#### Required arguments


Input: `-i`: Specifies the mOTUs profile files to merge. These files must be generated from the same mOTUs version and with the same parameters. At least two profiles are required. The input can be provided as a text file with one line per profile or as a space-separated list containing multiple mOTUs profiles.

Output: `-o`: Specifies the path to the merged profile file.

---

### download

```bash
$ motus download
```

<details>
<summary>download cli options</summary>

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


     motus download [options]

         Output options:
            -s  FILE  Genome metadata file
            -o  PATH  Genome output folder. Only required when
                      -l is not set

         Options:
            -l        Skip genome download. Only create genome report file
            -r        Download only representative genomes
            -w   STR  Keyword: Can be mOTU, genome name or taxonomy.
                        Fuzzy search enabled for taxonomy


motus.py: error: the following arguments are required: -s, -w
```
</details>



#### Required arguments


Input: `-w`: This keyword is used to search for GTDB taxonomy, mOTUs, and genome names. If no hits are found, a fuzzy search is performed on taxonomy.

Output:

* `-s`: This option provides a summary of files that are or could be downloaded.
* `-o`: This option specifies the output folder where downloaded genomes will be stored.  
  
#### 🔧 Options

| Option               | Description                                                                                                                                                                                                                                                                |
| -------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-w`     | **Keyword**: This keyword is used to search for GTDB taxonomy, mOTUs, and genome names. If no hits are found, a fuzzy search is performed on taxonomy.
| `-o`     | **Output folder**: This option specifies the output folder where downloaded genomes will be stored.
| `-s`     | **Summary file**: This option provides a summary of files that are or could be downloaded.
| `-l`     | **Skip download**: Skip genome download. Only write links to the summary file.
| `-r`     | **Representative only**: Instead of downloading all genomes, only download the representative per mOTU

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

            -f        Force download even when database is already present


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







## ❓ Need Help?

Write a issue on GitHub
