.PHONY: none testdata

comma := ,
empty :=
space := $(empty) $(empty)

ALLTARGETS = sample_data/ERR878216_both.bam \
	     sample_data/ERR878216_both.mgc \
	     sample_data/ERR878216_both.profile \
	     sample_data/ERR878216_both_counts.profile \
	     sample_data/ERR878216_both.biom \
	     sample_data/ERR878216_single_1.bam \
	     sample_data/ERR878216_single_1.mgc \
	     sample_data/ERR878216_single_1.profile \
	     sample_data/ERR878216_single_1.biom \
	     sample_data/ERR878216_single_2.bam \
	     sample_data/ERR878216_single_2.mgc \
	     sample_data/ERR878216_single_2.profile \
	     sample_data/ERR878216_single_2.biom \
	     sample_data/ERR878216_single_merge_1_2.biom \
	     sample_data/ERR878216_single_merge_1_2.profile \
 	     sample_data/ERR878216_both_snv.bam \
 	     sample_data/ERR878216_single_1_snv.bam

none:
	@echo "Run 'make testdata' to create/update all sample_data"
	@echo "Note you should only create/update sample_data when databases"
	@echo "change or bugs that affect the output are fixed."
	@echo "In any case visually inspect the result and use 'git diff' to see major changes"


testdata: $(ALLTARGETS)
	
sample_data/ERR878216_both.bam: sample_data/ERR878216_sample_1.fastq.gz sample_data/ERR878216_sample_2.fastq.gz
	@../motus profile -f $(word 1,$^) -r $(word 2,$^) -I $@ > /dev/null
sample_data/ERR878216_both.mgc: sample_data/ERR878216_sample_1.fastq.gz sample_data/ERR878216_sample_2.fastq.gz
	@../motus profile -f $(word 1,$^) -r $(word 2,$^) -M $@ > /dev/null
sample_data/ERR878216_both.profile: sample_data/ERR878216_sample_1.fastq.gz sample_data/ERR878216_sample_2.fastq.gz
	@../motus profile -f $(word 1,$^) -r $(word 2,$^) -o $@
sample_data/ERR878216_both_counts.profile: sample_data/ERR878216_sample_1.fastq.gz sample_data/ERR878216_sample_2.fastq.gz
	@../motus profile -f $(word 1,$^) -r $(word 2,$^) -c -o $@

sample_data/ERR878216_both.biom: sample_data/ERR878216_sample_1.fastq.gz sample_data/ERR878216_sample_2.fastq.gz
	@../motus profile -B -f $(word 1,$^) -r $(word 2,$^) -o $@
	
sample_data/ERR878216_single_1.bam: sample_data/ERR878216_sample_1.fastq.gz
	@../motus profile -s $^ -I $@ > /dev/null
sample_data/ERR878216_single_1.mgc: sample_data/ERR878216_sample_1.fastq.gz
	@../motus profile -s $^ -M $@ > /dev/null
sample_data/ERR878216_single_1.profile: sample_data/ERR878216_sample_1.fastq.gz
	@../motus profile -s $^ -o $@

sample_data/ERR878216_single_1.biom: sample_data/ERR878216_sample_1.fastq.gz
	@../motus profile -B -s $^ -o $@

sample_data/ERR878216_single_2.bam: sample_data/ERR878216_sample_2.fastq.gz
	@../motus profile -s $^ -I $@ > /dev/null
sample_data/ERR878216_single_2.mgc: sample_data/ERR878216_sample_2.fastq.gz
	@../motus profile -s $^ -M $@ > /dev/null
sample_data/ERR878216_single_2.profile: sample_data/ERR878216_sample_2.fastq.gz
	@../motus profile -s $^ -o $@

sample_data/ERR878216_single_2.biom: sample_data/ERR878216_sample_2.fastq.gz
	@../motus profile -B -s $^ -o $@

sample_data/ERR878216_single_merge_1_2.profile: sample_data/ERR878216_single_1.profile sample_data/ERR878216_single_2.profile
	../motus merge -i $(subst $(space),$(comma),$^) -o $@
sample_data/ERR878216_single_merge_1_2.biom: sample_data/ERR878216_single_1.profile sample_data/ERR878216_single_2.profile
	../motus merge -B -i $(subst $(space),$(comma),$^) -o $@

sample_data/ERR878216_both_snv.bam: sample_data/ERR878216_sample_1.fastq.gz sample_data/ERR878216_sample_2.fastq.gz
	@../motus map_snv -f $(word 1,$^) -r $(word 2,$^) -o $@
sample_data/ERR878216_single_1_snv.bam: sample_data/ERR878216_sample_1.fastq.gz
	@../motus map_snv -s $^ -o $@

