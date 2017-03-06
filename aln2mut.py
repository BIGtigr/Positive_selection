#!/usr/bin/python

import sys
import os
import csv
from Bio import SeqIO
import ntpath



# USAGE: python aln2mut.py < fd result > < alignment file fasta format > < output >

def main():

    fd = dict()# KEY: gene name VALUE: positions changed by functional divergence
    
    with open(sys.argv[1],"r") as f_in:
        gene = sys.argv[1].split("/")[-1]
        gene = gene.split(".")[0]
        for line in f_in:
            if "p <" in line:
                line = line.replace("\n","")
                res = line.split("\t")[2]
                pos_fd = line.split("\t")[0] # add position values to array
                pos_fd = pos_fd.replace(" ","")
                res_ref = res.split("/")[1]
                res_ref = res_ref.replace(" ","")
                res_ref = res_ref.split("(")[0]
                res_sk = res.split("/")[0]
                res_sk = res_sk.replace(" ","")
                res_sk = res_sk.split("(")[0]
                res_final = res_sk+"-"+res_ref
                fd[pos_fd] = res_final


    fasta = SeqIO.parse(sys.argv[2],"fasta")
    for f in fasta:
        if f.id == "s288c":
            for k in fd.keys():
                gap = 0
                pos = int(k)
                for s in range(0,int(k)):
                    if f.seq[s] == "-":
                        gap += 1

                pos = int(pos) - gap
                print gene+" "+fd[k].split("-")[1]+str(pos)+fd[k].split("-")[0]

                    

    fasta.close()


if __name__ == "__main__" :
    main()
