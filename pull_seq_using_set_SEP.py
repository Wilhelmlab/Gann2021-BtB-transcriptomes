#imports 
import os 
from collections import defaultdict #Slightly different dictionary, saves some time writting code
#location of files
path = 'C:/Users/egann/Desktop/viral_group_accesion'
files = os.listdir(path)

#making a dictionary mapping file to contigs, this way you get all contigs you need from one file access.
#This does use more memory however, so if you are collecting hundreds of million of identifiers, it might be a problem
#sets also make data access instant

for file in files: 
    file_contigs = defaultdict(set) 
    #'with' automatically will close the file for you
    with open(os.path.join(path,file),'r') as f:
        for line in f:
            contig = line.strip()
            #it is always best to do all the work you can over the same loop
            contig_file = contig.strip('>').split('_')[0] + '.txt'
            file_contigs['CM_Viruses.txt'].add(contig)

    #open output file 
    out_name = file.strip('.txt') + '_bt.txt'
    o = open(out_name,'w')

    #this way, you only open a file once and get all the contigs you need!
    #If this is still slow, we can use the multiprocessing module if need be
    for f in file_contigs:
    	#Open file based on start of contig name
        g = open(f,'r')#
        #search for line with contig name as first column
        for line in g:
            if line.strip().split('\t')[1].split('.')[0] in file_contigs[f]:
    			#write to output
                o.write(line)
        g.close()

    o.close()