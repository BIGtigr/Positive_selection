#!/usr/bin/python


# USAGE: python calc_p-value. py > results.csv (call script from the folder where codeml files are)

# Output: GENE NAME,p-value

# Read codeml outputs (alternative and null hypothesis) and get the p-value for every gene analized

'''
A significant result with the branch-site codon model means that
positive selection affected a subset of sites during a specific evolutionary time 

'''


import sys
import os
import re
from os.path import basename
from scipy.stats import chisqprob



def getLn(f_in):
    with open(f_in, "r") as file_in:
        contents = file_in.readlines()
    targets = [s for s in contents if "lnL(ntime" in s]
    lnl = re.findall("-\d+\.\d+", targets[0])[0] # find negative number in line lnL
    
    return lnl
    

def calcPvalue(ln_1,ln_2):
    if ln_1 and ln_2:
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
        pv = calcPvalue(ln_alt,ln_null)
        line = gene_name+","+str(pv)
        print line
        if pv < 0.05:
            with open("significativos.csv", "a") as o_file:
                o_file.write(line+"\n")

    #print "TOTAL SIGNIFICANT RESULTS :", positivos                    

if __name__ == "__main__" :
    main()
