#!/usr/bin/env python

import argparse

# accept filename as argument
parser = argparse.ArgumentParser(description="Get SAM filename")
parser.add_argument("-f", help="SAM file to run this script on", required=True)
parser.add_argument("-o", help="output", required=True)
args = parser.parse_args()
fname: str = args.f
outname: str = args.o

# holder for alignment section fields, bit flag
align_fields: list = []
flag: int = 0

# list of unmapped reads
unmappedr: dict = {}
# list of mapped reads
mappedr: dict = {}

# read SAM file
with open(fname, 'r') as sam:
    for line in sam:
        # only act on alignment section, not header lines
        if line[0] != '@':
            align_fields = line.split('\t')
            flag = int(align_fields[1])

            # check if read is mapped
            if((flag & 4) != 4):
                mapped = True
            else:
                mapped = False

            # count if mapped or unmapped, but don't count the same read multiple times
            if mapped:
                if((flag & 256) != 256):
                    if align_fields[0] in mappedr:
                        mappedr[align_fields[0]] += 1
                    else:
                        mappedr[align_fields[0]] = 1
            else:
                if align_fields[0] in unmappedr:
                    unmappedr[align_fields[0]] += 1
                else:
                    unmappedr[align_fields[0]] = 1

# output results
with open(outname, 'w') as output:
    output.write("Mapped reads: "+str(len(mappedr))+'\n')
    output.write("Unmapped reads: "+str(len(unmappedr)))