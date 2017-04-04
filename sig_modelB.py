#!/usr/bin/python


# USAGE: python sig_modelB.py  (call script from the folder where codeml files are)

# Output: GENE NAME omega value

# Read codeml outputs for MODEL B analysis (alternative and null hypothesis) and get the omega value for both Sk Sc branches

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
    targets = [s for s in contents if "lnL" in s]
    lnl = re.findall("-\d+\.\d+", targets[0])[0] # find negative number in line lnL
    
    return lnl


def getOmega(f_in):
    with open(f_in, "r") as file_in:
        contents = file_in.readlines()
    targets = [s for s in contents if "w (dN/dS) for branches:" in s]
    omega = []
    om_values = []
    omega = re.findall("\d+\.\d+", targets[0]) # find negative number in line lnL
    for x in omega:
        om_values.append(float(x))
    
    return om_values
    

def calcPvalue(ln_1,ln_2):
    if ln_1 and ln_2:
        val = 2*(float(ln_1)-(float(ln_2)))
        p_val = chisqprob(val, 2)
    else:
        p_val = "NA"

    return p_val



def main():

    path = os.getcwd()
    alt_files = [f for f in os.listdir(path) if f.endswith("_alt.txt")]
    null_files = [f for f in os.listdir(path) if f.endswith("_null.txt")]

    positivos = 0
    for f in alt_files:
        #print f
        gene_name = f.split("_")[0]
        alt = f
        null = gene_name + "_null.txt"
        ln_null = getLn(null)
        ln_alt = getLn(alt)
        pv = calcPvalue(ln_alt,ln_null)
        om_alt = getOmega(alt)
        if pv < 0.05:
            #elif om_alt[1] > 1 and om_alt[2] <= 1:
             #   print gene_name, om_alt[1] # Sk omega
            elif om_alt[2] > 1 and om_alt[1] <= 1:
                print gene_name, om_alt[2] # Sc omega
        

    #print "TOTAL SIGNIFICANT RESULTS :", positivos                    

if __name__ == "__main__" :
    main()
