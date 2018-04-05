#!/usr/bin/env python
import sys
import glob
import argparse
from collections import defaultdict


__author__ = "Thomas Hansen (tbh@mbg.au.dk)"
__lab__ = "ncRNAlab"
__version__ = '1.0.1'



circglob = {
    "acsf" : "**/circle_candidates_expr",
    "circrna_finder" : "**/s_filteredJunctions.bed",
    "circexplorer" : "*circ.txt",
    "circexplorer2" : "*.circ.txt",
    "ciri" : "*.ciri.out.txt",
    "ciri2" : "*.ciri.out.txt",
    "dcc" : "**/*",
    "find_circ" : "*.circ_candidates.bed",
    "find_circ_40" : "*.circ_candidates.bed",
    "knife" : "**/circReads/combinedReports/*_circJuncProbs.txt",
    "mapsplice" : "**/circular_RNAs.txt",
    "uroborus" : "**/circRNA_list.txt",
}


class circobject (object):

    def __init__(self, algorithm, csample = 2, cmerged = 2):
               
        self.algorithm = algorithm.lower()
        self.cutoff_sample = csample
        self.cutoff_merged = cmerged
        self.dReads = defaultdict (list)
        self.dPos = defaultdict (int)
        self.dStrand = defaultdict (str) #only dcc

    def get_files (self, files):


        if len(files) == 0:

            f = []
            for file in glob.glob(circglob[self.algorithm]):
                f.append (file)

            return f


        else:
            return files
            
    def add_circ (self, file, line):
        
        cols = line.strip().split("\t")

        if cols[0].strip() == '':
            return

        if cols[0][0] == "#":
            return
                
        if len (cols) < 4:
            return

        while len (cols) < 6:
            cols.append ("NA")


        chrom, start, end, bsj, strand = cols[0], cols[1], cols[2], cols[4], cols[5] # default bed 
        lin = 0

        ############# ACSF ####################
        if self.algorithm == "acsf":

            if cols[0] == "newid":
                return
            
            pos = cols[0].split ("_")
            
            if len (cols) >= 4:
                bsj = int (cols[2]) + int (cols[3])      
            else:
                bsj = int (cols[2])
                        
            chrom = "chr" + pos[0]
            si = len(pos) - 3
            ei = si+1
            start, end = str (min (int(pos[si]), int(pos[ei]))), str (max (int(pos[si]),int(pos[ei]))) 
            strand = pos[ei+1][0]
            
        ############# circrna_finder ####################
        elif self.algorithm == "circrna_finder":
            bsj, strand = cols[4], "NA"

        ############# circexplorer ####################
        elif self.algorithm == "circexplorer":
            
            bsj = cols[12]  # bsj[bsj.find ('/')+1:]
            #bsj = cols[3]
            #bsj = bsj[bsj.find ('/')+1:]

            if len (cols) > 13 and cols[13] == "Yes": #ciRNA
                return
        
        ############# circexplorer2 ####################
        elif self.algorithm == "circexplorer2":
            
            bsj = cols[3]
            bsj = bsj[bsj.find ('/')+1:]
                    
        ############# ciri/ciri2 ####################
        elif self.algorithm == "ciri" or self.algorithm == "ciri2":
            
            chrom, start, end, bsj, lin, strand = cols[1], cols[2], cols[3], cols[4], cols[6], cols[10]
                  
            if strand != "+" and strand != "-":
                strand = "NA"
            
            if lin == "":
                lin = '0'

        ############# ciri/ciri2 ####################
        elif self.algorithm == "dcc":

            if file.find ("Coordinates") != -1:

                junction_type = int(cols[4])            
                if junction_type < 1 or junction_type > 2:
                    return
                self.dStrand[(chrom, start, end)] = strand    
                bsj, strand = 0, "NA"    
                
            elif file.find ("RNACount") != -1: 
                   
                bsj = cols[3]
                    
    
        ############# find_circ (mapq 40) ####################
        elif self.algorithm == "find_circ_40":
            
            if len (cols) < 21 and (int (cols[8]) < 40 or int (cols[7]) < 40):  
                return

            elif int (cols[9]) < 40 or int (cols[8]) < 40:  
               return

            
        ############# knife ####################
        
        elif self.algorithm == "knife":
      
            pos_info = cols[0].split ('|')
            chrom, pos1, pos2 = pos_info[0], int (pos_info[1].split (':')[1]), int (pos_info[2].split (':')[1])
            start, end = str(min (pos1, pos2)), str(max (pos1, pos2))
            strand = pos_info[4]
            
            bsj = cols[5].strip ()

            if bsj == '0':
                return
            #p_value = max (float (cols[2]), float(cols[4]))
                  
        ############# mapsplice ####################
        elif self.algorithm == "mapsplice":
            
            chrom = chrom[:chrom.find ('~')]            
            start, end = str (min (int(cols[1]), int(cols[2]))), str (max (int(cols[1]),int(cols[2]))) 
            bsj = cols[4]
            strand = cols[5][0]
                        
        elif self.algorithm == "uroborus":
            bsj = cols[6]
            strand = cols[3]


            
        if not start.isdigit () or not end.isdigit ():
            return

        
        b,l = 0,0
        if (chrom, start, end, strand, file) in self.dReads:
            (b, l) = self.dReads[(chrom, start, end, strand, file)] # used in circrna_finder

        if int (bsj) + int(b) == 0 or int (bsj) + int(b)  < self.cutoff_sample:
            return

        self.dReads[(chrom, start, end, strand, file)] = (int(bsj)+int(b), int(lin)+int(l))
        self.dPos[(chrom, start, end, strand)] = 1
        

    


    def output (self, files, onlyBsj = True):

        output = ["#chrom", "start", "end", "name", "score", "strand"]

        if self.algorithm == "dcc":
            files = [f for f in files if f.find ("RNACount") != -1]

        for f in files:

            if onlyBsj:
                output += [f + "_BSJ"]
            else:
                output += [f + "_BSJ", f + "_LIN"]
                

        print "\t".join (output)

        circ_index = 1

        for (chrom, start, end, strand) in self.dPos:
                    
            output_data = []
                
            total_bsj = 0  
            
            for f in files:
            
                bsj,lin = 0,0
            
                if (chrom, start, end, strand, f) in self.dReads:
                    (bsj,lin) = self.dReads[(chrom, start, end, strand, f)]         

                if onlyBsj:
                    output_data += [str(bsj)]
                else:
                    output_data += [str(bsj), str(lin)]

                total_bsj = total_bsj + int(bsj)
            
            
            name = "_".join (["circ", str(circ_index)])
               
            if self.algorithm == "dcc":
                strand = self.dStrand[(chrom, start, end)]


            output = [chrom, str(start), str(end), name, str(total_bsj), strand] + output_data
        
            if total_bsj >= self.cutoff_merged:  
                print "\t".join (output)
                circ_index += 1

        return circ_index-1       




algorithms = [a for a in circglob]
algorithms.sort ()

parser = argparse.ArgumentParser(description='circM: Merge circRNA predictions')

parser.add_argument("-f", nargs="+", dest="files", default=[],type=str,help="Input files")
parser.add_argument("-cm", dest="cutoff_merged",type=int,default=2,help="Merged expression cutoff")
parser.add_argument("-cf", dest="cutoff_file",type=int,default=2,help="Expression cutoff for each sample")
parser.add_argument("-a", dest="algorithm",type=str,choices=algorithms,help="Prediction algorithm")

args = parser.parse_args()

circ = circobject (args.algorithm, args.cutoff_file, args.cutoff_merged)

files = circ.get_files (args.files)

files.sort ()


print >>sys.stderr, "Algorithm...:", args.algorithm
print >>sys.stderr, "Cutoff...:"
print >>sys.stderr, "\tcircRNA must have at least", args.cutoff_file, "reads in one sample to be included"
print >>sys.stderr, "\tcircRNA must have at least", args.cutoff_merged, "reads in total to be included"
print >>sys.stderr, "Found...:", len(files), "files for analysis"


circ.output   

for file in files:
   
    print >>sys.stderr, "Analyzing...:", file  
    f = open(file, 'r')
   
    for line in f:

        circ.add_circ (file, line)

        
ncirc = circ.output(files)

print >>sys.stderr, "Found", ncirc, "circRNAs"

       
    