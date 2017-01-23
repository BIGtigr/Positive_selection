#!/usr/bin/python


# USAGE: python calc_p-value.py > results.csv (run script in the folder where codeml files are)

# Output: GENE NAME,p-value

# Read codeml outputs (alternative and null hypothesis) and get the p-value for every gene analized

'''
A significant result with the branch-site codon model means that
positive selection affected a subset of sites during a specific evolutionary time 

'''


import sys
import os
from os.path import basename
from scipy.stats import chisqprob


# Get the Ln value from output codeml file
def getLn(f_in):
    with open(f_in, "r") as file_in:
        if os.stat(f_in).st_size > 0: # if file is not empty
            for line in file_in:
                if line.startswith("lnL"):
                    lnl = line.split("  ")[2]
                    break

        else:
            lnl = "NA"
    return lnl

# Calculate p-value. Chi-square distribution with 1 degree of freedom
def calcPvalue(ln_1,ln_2):
    if ln_1 and ln_2 != "NA":
        val = 2*(float(ln_1)-(float(ln_2)))
        p_val = chisqprob(val, 1)
    else:
        p_val = "NA"

    return p_val



def main():

    path = os.getcwd()
    alt_files = [f for f in os.listdir(path) if f.endswith("_alt.txt")]
    null_files = [f for f in os.listdir(path) if f.endswith("_null.txt")]

    positivos = 0
    for f in alt_files:
        gene_name = f.split("_")[0]
        alt = f
        null = gene_name + "_null.txt"
        ln_alt = getLn(alt)
        ln_null = getLn(null)
        p_val = calcPvalue(ln_alt,ln_null)
        line = gene_name+","+str(p_val)
        print line
        if p_val < 0.05:
            positivos += 1

    print "TOTAL SIGNIFICANT RESULTS :", positivos                    

if __name__ == "__main__" :
    main()
