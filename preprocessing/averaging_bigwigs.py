#!/usr/bin/env python

'''This program will take the output of multiBamSummary --outRawCounts.
It will use all columns after the third one to calculate a median value,
which is returned. 
Possible improvements to this script include:
- accept user-supplied header instead of an example bigWig, e.g. obtained via
awk '{OFS="\t";print $1,$2,$2+1}' GRCm38.primary_assembly.genome.fa.fai | \
    egrep ^chr  | bedtools sort -i | cut -f1,2 > mm38.fa.fai
    This would ensure the same sorting as for the tab file.
- filter tab lines based on the chromosomes, i.e., discard every line
that does not start with a chromosome known from the header

'''

import sys
import argparse
import getopt
from statistics import median
from statistics import mean
import re
import pyBigWig

## define the arguments ================================================

def get_args():
    
    parser = argparse.ArgumentParser(description='Merge multiple bigWig files into one.')
    
    parser.add_argument('--infile', '-in', type = str, required=True, help="The outRawCounts output of multiBamSummary, i.e. a tab-separated file where the first three columns correspond to chr, star, end. Should be sorted, can be gzipped.")
    parser.add_argument('--summaryType', '-st', type = str, default = 'median', help="What type of operation should be used for calculating the merged value. Choices are: mean, median, sum, min, max. Default: median.")
    parser.add_argument('--outfile', '-out', type = str, required=True, help = 'Prefix for the output file')
    parser.add_argument('--outfileType', '-ot', type = str, default='bedGraph', help = 'Indicate whether a bigWig or a bedGraph file should be returned. Default: bedGraph')
    parser.add_argument('--oriBigWig', type = str, required=False, help = 'If you choose to output a bigWig file, indicate here the name of one of the bigWig files that was the basis for the multiBigWigSummary result. This will be used to extract the chromosome information that will be stored in the resulting bigWig.')
        
    args=parser.parse_args()
    
    return args


## fct for creating bigWig file=========================================
def initBigWig(exampleBW, outname):
    """!
    Initiates a bigWig file
    
    @param exampleBW String: name of the bigWig file from which to use the header info
    @param outname String: prefix for the file name of the bigWig file to be created

    @return Opens a bigWig file and adds a header so that it's ready to
    receive more information.
    """
    
    # read the header info from an example bigWig
    bw_check= pyBigWig.open(exampleBW)
    # extract the chromosome information
    chrom_info = bw_check.chroms()
    # wrangle the chrom info into a list of tuples
    chrom_info = [tuple(i) for i in chrom_info.items()]
    
    #open a new bigWig file
    bw_out = pyBigWig.open(outname + '.bw', "w")   
    bw_out.addHeader(chrom_info)

    return bw_out
    

## main function =======================================================
def main():
    
    args = get_args()
    
    if ".gz" in args.infile:
        in_tab = gzip.open(args.infile, "r")
    else:
        in_tab = open(args.infile, "r")
   
   
    if args.outfileType == 'bigWig':
        out = initBigWig(args.oriBigWig, args.outfile)
    else:
        out = open(args.outfile + '.bg', 'w')
    
    for x in in_tab.readlines():
        x = x.strip().split("\t") # get the individual columns
    
        if not bool(re.match(r'#',x[0])): # ignore commented lines
            coords = x[0:3]
            x = x[3:] # drop first three columns
            xf = [float(s) for s in x] # turn strings into floats
            
            if args.summaryType == 'mean':
                x_out = mean(xf)
            elif args.summaryType == 'median':
                x_out = median(xf)
            elif args.summaryType == 'sum':
                x_out = sum(xf)
            elif args.summaryType == 'min':
                x_out = min(xf)
            elif args.summaryType == 'max':
                x_out = max(xf)
                
            if args.outfileType == 'bedGraph':
                out.write( "{c}\t{s}\t{e}\t{v:9.6f}\n".format(c=coords[0], s=coords[1], e=coords[2], v=x_out) )
            elif args.outfileType == 'bigWig':
                print(coords)
                print(x_out)
                out.addEntries( [coords[0]], [int(coords[1])], ends = [int(coords[2])], values = [x_out])
            
    
    out.close()

if __name__ == '__main__':
    main()
