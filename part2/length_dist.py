#!/usr/bin/env python

import matplotlib.pyplot as plt
import argparse
import gzip
import numpy as np

# QAA Part 2 Question 7
# Plot the trimmed read length distributions for both R1 and R2 reads (on the same plot)

# accept filename as argument
parser = argparse.ArgumentParser(description="Get zipped FASTQ filename")
parser.add_argument("-r1", help="Read 1 zipped FASTQ file", required=True)
parser.add_argument("-r2", help="Read 2 zipped FASTQ file", required=True)
parser.add_argument("-o", help="Output file path and name", required=True)
args = parser.parse_args()
read1_path: str = args.r1
read2_path: str = args.r2
outname: str = args.o

# test files
# read1_path = "/projects/bgmp/bmeluch/bioinfo/Bi622/Demultiplex/Assignment-the-first/test.fastq"
# read2_path = "/projects/bgmp/bmeluch/bioinfo/Bi622/Demultiplex/Assignment-the-first/test.fastq"
# outname = "/projects/bgmp/bmeluch/bioinfo/Bi623/QAA/part2/test_plot.png"

# init line counter
lcount: int = 0
# init sequence length
seqlen: int = 0

# lengths dicts
read1_lengths: list = [0]*102
read2_lengths: list = [0]*102

# get sequence lengths
with gzip.open(read1_path, 'rt') as read1:
#with open(read1_path, 'rt') as read1:
    for line in read1:
        lcount += 1
        if lcount%4 == 2 or lcount == 2:
            seqlen = len(line.strip())
            read1_lengths[seqlen] += 1

with gzip.open(read2_path, 'rt') as read2:
#with open(read2_path, 'rt') as read2:
    for line in read2:
        lcount += 1
        if lcount%4 == 2 or lcount == 2:
            seqlen = len(line.strip())
            read2_lengths[seqlen] += 1

# plot histogram
plt.title("Read length distribution")
plt.xlabel("Read length")
plt.ylabel("Number of reads")
plt.yscale("log")
x = np.arange(102)
y1 = read1_lengths
y2 = read2_lengths
width = 0.40
plt.bar(x-width/2, y1, width)
plt.bar(x+width/2, y2, width)
plt.legend(["Read 1", "Read 2"])
plt.savefig(outname)
plt.show()