#imports 
import os 
from collections import defaultdict #Slightly different dictionary, saves some time writting code
#location of files

#making a dictionary mapping file to contigs, this way you get all contigs you need from one file access.
#This does use more memory however, so if you are collecting hundreds of million of identifiers, it might be a problem
#sets also make data access instant
file_contigs = defaultdict(set) 
#'with' automatically will close the file for you
with open('contigs_with_top_diamond_viral.txt','r') as f:
    for line in f:
        contig = '>' + line.strip()
        #it is always best to do all the work you can over the same loop
        contig_file = contig.strip('>').split('_')[0] + '.txt'
        file_contigs['CM_all.txt'].add(contig)

#open output file 
o = open('Virus_top_hit_PolyA.fasta','w')

#this way, you only open a file once and get all the contigs you need!
#If this is still slow, we can use the multiprocessing module if need be
for f in file_contigs:
	#Open file based on start of contig name
    g = open(f,'r')#
    #search for line with contig name as first column
    for line in g:
        if line.strip().split('\t')[0] in file_contigs[f]:
			#write to output
            o.write(line.strip().split('\t')[0])
            o.write('\n')
            o.write(line.strip().split('\t')[1].replace('"',''))
            o.write('\n')
    g.close()

o.close()