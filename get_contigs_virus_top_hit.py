#imports 
import os 
from collections import defaultdict

#open virus accessions 
g = open('Viral_Protein_Accessions.txt','r')

virus_accessions = set()

for line in g:
	virus_accessions.add(line.strip())

g.close()

#open out contig list file 
o = open('contigs_with_top_diamond_viral.txt','w')

#path to diamond blast tables 
db_path = 'C:/Users/egann/Desktop/working/db_tables'

#get all files in the path 
db_files = os.listdir(db_path)

#make smaller temp tables by pulling out all fasta headers 
for file in db_files:
	out = []
	#small table 
	accessions_seen = set() #holds lines already seen 
	smaller_table = []
	with open(os.path.join(db_path,file),'r') as f: 
		#if seen accession do not write line  
		for line in f: 
			if line.strip().split('\t')[0] not in accessions_seen:
				accessions_seen.add(line.strip().split('\t')[0])
				smaller_table.append(line)

	for line in smaller_table:
		if line.strip().split('\t')[1].split('.')[0] in virus_accessions:
			out.append(line.strip().split('\t')[0])


	for accession in out:
		o.write(accession)
		o.write('\n')

	print(file)