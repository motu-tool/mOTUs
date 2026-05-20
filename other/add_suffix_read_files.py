
import pysam
import sys
import gzip


r1_in = sys.argv[1]
r2_in = sys.argv[2]
r1_out = sys.argv[3]
r2_out = sys.argv[4]


def read_fastq(infile):
    d = []
    with pysam.FastxFile(infile) as fh:
        for entry in fh:
            name =entry.name
            sequence = entry.sequence
            quality = entry.quality
            d.append((name, sequence, quality))
    return d

def write_fastq(d, outfile, suffix):
    with gzip.open(outfile, 'wt') as handle:
        for (name, sequence, quality) in d:
            handle.write(f'@{name}{suffix}\n{sequence}\n+\n{quality}\n')


write_fastq(read_fastq(r1_in), r1_out, '/1')
write_fastq(read_fastq(r2_in), r2_out, '/2')