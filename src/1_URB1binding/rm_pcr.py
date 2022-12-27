#!/bin/python
# -*- coding: utf-8 -*-
import sys
import pysam

def main():
    if len(sys.argv) != 3:
        sys.exit('rm_pcr.py barcode bam')
    barcode = {}
    with open(sys.argv[1], 'r') as barcode_f:
        for line in barcode_f:
            name = line.split()[0][1:]
            tag = line.split()[1]
            barcode[name] = tag
    prefix = sys.argv[2].split('.')[0]
    bam_f = pysam.Samfile(sys.argv[2], 'rb')
    out = pysam.Samfile(prefix + '_combine.bam', 'wb', header=bam_f.header)
    start, end = '', ''
    duplication = {}
    for read in bam_f.fetch():
        if read.pos != start and read.aend != end:
            remove_duplication(duplication, barcode, out)
            duplication = {}
            start = read.pos
            end = read.aend
        duplication[read.qname] = read
    remove_duplication(duplication, barcode, out)

def remove_duplication(duplication, barcode, out):
    tag_set = {}
    for name in duplication:
        tag_set[barcode[name]] = name
    for tag in tag_set:
        read = duplication[tag_set[tag]]
        out.write(read)

if __name__ == '__main__':
    main()