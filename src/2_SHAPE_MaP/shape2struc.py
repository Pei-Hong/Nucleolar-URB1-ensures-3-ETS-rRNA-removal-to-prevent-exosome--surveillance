#!/usr/bin/env python

import os
import subprocess
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--shape", type=str, help="input the mapped sam file")
parser.add_argument("-f", "--fa", type=str, help="input the fa file")
parser.add_argument("--circ",action='store_true', help="see the seq as circular RNA when predict RNA structure  ")
parser.add_argument("-o", "--out", type=str, help="output the paired fastq file")
args = parser.parse_args()


Name=args.fa.split("/")[-1].split(".")[0]
print Name
RNAfold_linear= "RNAfold -p -d2 --shape="+args.shape+" --shapeMethod=D <"+args.fa+">"+args.out+".out"
RNAfold_circ = "RNAfold -p -d2 --circ --shape="+args.shape+" --shapeMethod=D <"+args.fa+">"+args.out+".out"

#draw="relplot.pl "+Name+"_ss.ps "+Name+"_dp.ps > "+Name+"_rss.ps"

if args.circ is True:
	print RNAfold_circ
	os.popen(RNAfold_circ)
#	print draw
#	os.popen(draw)
else:
	print RNAfold_linear
	os.popen(RNAfold_linear)
#	print draw
#	os.popen(draw)

fileIn=open(args.out+".out")
dotOut=open(args.out+".dot",'w')
for line in fileIn:
    name=line
    seq=next(fileIn)
    dot=next(fileIn).split()[0]
    dotOut.write(name+seq+dot)
    break
fileIn.close()
dotOut.close()


dot2ct="dot2ct "+args.out+".dot"+" "+args.out+".ct"
print dot2ct
os.popen(dot2ct)
#draw="pvclient.py --ct "+Name+args.out+".ct"+" --shape "+args.shape+" --structures 1 --out "+Name+args.out
#print draw
#os.popen(draw)