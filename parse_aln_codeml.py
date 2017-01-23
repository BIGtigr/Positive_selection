#!/usr/bin/python



from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys
import re
import os


def get_header(seq_id):
    header = seq_id.split("/")[0]

    return header


def main():

    f_in = SeqIO.parse(sys.argv[1], "fasta")
    gene_name = sys.argv[1].split(".")[0]
    

    seqs = []
    headers = []
    

    for rec in f_in:
        h = get_header(rec.id)
        headers.append(h)
        seqs.append(str(rec.seq))
        #print h


    with open(gene_name+".phy", "a") as f_out:
        f_out.write("  "+str(len(headers))+"  "+str(len(seqs[0]))+"  "+"\n")
        f_out.write(headers[0]+"  "+seqs[0]+"\n")
        f_out.write(headers[1]+"  "+seqs[1]+"\n")
        f_out.write(headers[2]+"  "+seqs[2]+"\n")
        f_out.write(headers[3]+"  "+seqs[3]+"\n")
        f_out.write(headers[4]+"  "+seqs[4]+"\n")
        f_out.write(headers[5]+"  "+seqs[5]+"\n")
        f_out.write("C.glabrata"+"  "+seqs[6]+"\n")
        
    


    f_in.close()


if __name__ == "__main__" :
    main()
