#!/usr/bin/python

import sys
import os

def commonElements(list_a, list_b):
    common = list(set(list_a).intersection(list_b))
    
    return common


def main():
    # List both FD codeml files
    path_fd = os.path.dirname(sys.argv[1])
    fd_files = [f for f in os.listdir(path_fd) if f.endswith(".txt") and f.startswith("Y")]
    path_pv = os.path.dirname(sys.argv[2])
    pv_files = [z for z in os.listdir(path_pv) if z.endswith("_alt.txt") and z.startswith("Y")]


    fd = dict()# KEY: gene name VALUE: positions changed by functional divergence
    # read number of changes in FD files
    for f in fd_files:
        gene = f.split(".")[0]
        fd.setdefault(gene,[])
        pos_fd = [] 
        with open(os.path.join(path_fd,f),"r") as f_in:
            for line in f_in:
                if "p <" in line:
                    line = line.replace("\n","")
                    pos_fd.append(line.split("\t")[0]) # add position values to array
        
        for p in pos_fd:
            p = p.replace(" ","")
            fd[gene].append(p) # add position values to dictionary with gene name as key
        if len(fd) == 0:
            fd[gene].append("0") # add value zero if there is no FD change in that protein
            
    


    sel = dict() # KEY: gene name VALUE: positions under positive selection
    for h in pv_files:
        gene = h.split("_")[0]
        sel.setdefault(gene,[])
        pos_pv = []
        i = True
        # Read specific codeml output until find BEB results
        with open(os.path.join(path_pv,h),"r") as f_in: 
            for line in f_in:
                if line == "Positive sites for foreground lineages Prob(w>1):\n":
                    line = f_in.next()
                    if line != "\n":
                        while i:
                            pos_pv.append(line) # add specific lines to matrix
                            line = f_in.next()
                            if line == "\n":
                                i = False  # stop reading when there is a \n


        res = []
        for x in pos_pv:   # format position values, keeping only integers from line
            x = x.replace("\n","")
            x = x.replace("    ","")
            x = x.replace("   ","")
            res.append(x.split(" ")[0])

        for r in res:
            sel[gene].append(r) # add position values to dictionary with gene name as key


    for k in sel.keys():
        val_sel = []
        val_fd = []
        common = []
        if k in fd.keys():
            val_sel = sel.get(k)
            val_fd = fd.get(k)
            
        common = commonElements(val_sel,val_fd)
        common = [int(x) for x in common]
        common.sort()

        common = commonElements(val_sel,val_fd)
        if len(common) != 0:
            with open("overlapping_FD.csv","a") as f_out:
                f_out.write(k+",")
                for val in common:
                    f_out.write(str(val)+",")
                f_out.write("\n")

            



                    



if __name__ == "__main__" :
    main()

