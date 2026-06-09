

# download toy db
motus downloadMGDB --toy -db motus4.1-toy-db

# motus with one set of paired data - default params
motus profile -f input/ERR4507416/ERR4507416_1.motus.fastq.gz -r input/ERR4507416/ERR4507416_2.motus.fastq.gz -n ERR4507416 -t 1 -db motus4.1-toy-db/ -o output/ERR4507416-default





# motus with one set of paired data - different params


## allow /1 and /2 as read name suffixes


motus profile -f input/ERR4507418-suffix/ERR4507418_1.suffix.motus.fastq.gz -r input/ERR4507418-suffix/ERR4507418_2.suffix.motus.fastq.gz -n ERR4507418 -t 1 -db motus4.1-toy-db/ -o output/ERR4507418-withsuffix


## allow ignore randomly removed reads

motus profile -f input/ERR4507418-randomremovedreads/ERR4507418_1.rrr.motus.fastq.gz -r input/ERR4507418-randomremovedreads/ERR4507418_2.rrr.motus.fastq.gz -n ERR4507416 -t 1 -db motus4.1-toy-db/ -o output/ERR4507416-unbroken --skip-pair-check 


## -g

motus profile -f input/ERR4507416/ERR4507416_1.motus.fastq.gz -r input/ERR4507416/ERR4507416_2.motus.fastq.gz -n ERR4507416 -t 1 -db motus4.1-toy-db/ -o output/ERR4507416-g1 -g 1

motus profile -f input/ERR4507416/ERR4507416_1.motus.fastq.gz -r input/ERR4507416/ERR4507416_2.motus.fastq.gz -n ERR4507416 -t 1 -db motus4.1-toy-db/ -o output/ERR4507416-g10 -g 10

## -l

motus profile -f input/ERR4507416/ERR4507416_1.motus.fastq.gz -r input/ERR4507416/ERR4507416_2.motus.fastq.gz -n ERR4507416 -t 1 -db motus4.1-toy-db/ -o output/ERR4507416-l150 -l 150

motus profile -f input/ERR4507416/ERR4507416_1.motus.fastq.gz -r input/ERR4507416/ERR4507416_2.motus.fastq.gz -n ERR4507416 -t 1 -db motus4.1-toy-db/ -o output/ERR4507416-l40 -l 40



## -y
motus profile -f input/ERR4507416/ERR4507416_1.motus.fastq.gz -r input/ERR4507416/ERR4507416_2.motus.fastq.gz -n ERR4507416 -t 1 -db motus4.1-toy-db/ -o output/ERR4507416-yINSERT_RAW -y INSERT_RAW

motus profile -f input/ERR4507416/ERR4507416_1.motus.fastq.gz -r input/ERR4507416/ERR4507416_2.motus.fastq.gz -n ERR4507416 -t 1 -db motus4.1-toy-db/ -o output/ERR4507416-yINSERT_NORM -y INSERT_NORM

motus profile -f input/ERR4507416/ERR4507416_1.motus.fastq.gz -r input/ERR4507416/ERR4507416_2.motus.fastq.gz -n ERR4507416 -t 1 -db motus4.1-toy-db/ -o output/ERR4507416-yBASE_RAW -y BASE_RAW

motus profile -f input/ERR4507416/ERR4507416_1.motus.fastq.gz -r input/ERR4507416/ERR4507416_2.motus.fastq.gz -n ERR4507416 -t 1 -db motus4.1-toy-db/ -o output/ERR4507416-yBASE_NORM -y BASE_NORM



# single-end

motus profile -s input/ERR4507416/ERR4507416_1.motus.fastq.gz -n ERR4507416 -t 1 -db motus4.1-toy-db/ -o output/ERR4507416-default-single



# run second dataset with default parameters

motus profile -f input/ERR4507418/ERR4507418_1.motus.fastq.gz -r input/ERR4507418/ERR4507418_2.motus.fastq.gz -n ERR4507418 -t 1 -db motus4.1-toy-db/ -o output/ERR4507418-default


# merge the two default runs

motus merge -i output/ERR4507418-default output/ERR4507416-default -o output/merge-default -db motus4.1-toy-db/ 



# tests that should break

## merge two motus results that were created from different paramaters

motus merge -i output/ERR4507418-default output/ERR4507416-g10 -o output/merge-broken -db motus4.1-toy-db/

motus merge -i output/ERR4507418-default output/ERR4507416-yINSERT_NORM -o output/merge-broken -db motus4.1-toy-db/


motus profile -f input/ERR4507418-randomremovedreads/ERR4507418_1.rrr.motus.fastq.gz -r input/ERR4507418-randomremovedreads/ERR4507418_2.rrr.motus.fastq.gz -n ERR4507416 -t 1 -db motus4.1-toy-db/ -o output/ERR4507416-broken


