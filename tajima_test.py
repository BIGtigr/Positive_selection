#!/usr/bin/python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from scipy.stats import chisqprob
import sys
import numpy as np
import os


'''
Programa para realizar el test de Tajima aplicado a dos especies que se comparen con respecto a una especie outgroup

INPUT: archivo secuencia del gen alineada en formato fasta. Debe estar ordenado: primero una especie, luego
la otra y al final del todo la secuencia outgroup
'''



def main():

    # Leer alineamiento en formato fasta y repartirlo en especies y outgroup
    # Cambiar valores si hay más o menos de 4 cepas por especie
    suffix = sys.argv[1] # sufijo de los archivos de alineamiento
    path = os.getcwd()
    aln_files = [f for f in os.listdir(path) if f.endswith(suffix)]

    for f in aln_files:
        gene = f.split(".")[0]
        alignment = AlignIO.read(f,"fasta")
        al1 = alignment[:4]
        al2 = alignment[4:8]
        out = alignment[-1]

        # obtener una secuencia consenso para cada especie para despues comparar con el outgroup
        summary_align1 = AlignInfo.SummaryInfo(al1)
        summary_align2 = AlignInfo.SummaryInfo(al2)
        cons1 = summary_align1.dumb_consensus()
        cons2 = summary_align2.dumb_consensus()
        m1 = 0
        m2 = 0
        for x in range(1,(len(out.seq))):
            nt1 = cons1[x]
            nt2 = cons2[x]
            nt3 = out.seq[x]
            if nt1 != "X" and nt2 != "X" and nt3 != "X":
            # caso 1 nt2 = nt3 pero distinto de nt1
                if nt2 == nt3 and nt2 != nt1:
                    m1 += 1
            # caso 2 nt1 = nt3 pero distinto de nt2
                if nt1 == nt3 and nt1 != nt2:
                    m2 += 1

    # calcular valor chi-cuadrado a partir de m1 y m2 y ver si es significativo
        val = ((m1-m2)**2)*1.0/(m1+m2)
        p_val = chisqprob(val, 1)
        print gene,p_val

    
if __name__ == "__main__":
    main()

