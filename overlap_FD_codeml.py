#!/usr/bin/python

import sys
import os

def main():
    # List both FD codeml files
    path_fd = os.path.dirname(sys.argv[1])
    fd_files = [f for f in os.listdir(path_fd) if f.endswith(".txt") and f.startswith("Y")]
    path_pv = os.path.dirname(sys.argv[2])
    pv_files = [z for z in os.listdir(path_pv) if z.endswith("_alt.txt") and z.startswith("Y")]


    '''
    # read number of changes in FD files
    for f in fd_files:
        gene = f.split(".")[0]
        c = 0
        pos_fd = []
        with open(os.path.join(path_fd,f),"r") as f_in:
            for line in f_in:
                if "p <" in line:
                    line = line.replace("\n","")
                    pos_fd.append(line.split("\t")[0])

    '''

    sel = dict() # KEY: gene name VALUE: positions under positive selection
    for h in pv_files:
        pos_pv = []
        gene = h.split("_")[0]
        i = True
        # Read specific codeml output until find BEB results
        with open(os.path.join(path_pv,h),"r") as f_in: 
            for line in f_in:
                if line == "Positive sites for foreground lineages Prob(w>1):\n":
                    line = f_in.next()
                    if line != "\n":
                        while i:
                            pos_pv.append(line)
                            line = f_in.next()
                            if line == "\n":
                                i = False


        res = []
        sel.setdefault(gene,[])
        for x in pos_pv:
            x = x.replace("\n","")
            x = x.replace("    ","")
            x = x.replace("   ","")
            residuos.append(x.split(" ")[0])

        for r in res:
            sel[gene].append(r)


                    



if __name__ == "__main__" :
    main()
