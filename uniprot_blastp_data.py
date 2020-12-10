#Make a table from the blastp against the uniprot/swiss prot 
#database with the protein file, the protein from uniprot/swiss prot
#and the name of the protein, only if the percent identity > 30% 
#and the evalue < 1e-50 


#imports 
import os 
import csv 

#in directory 
directory = 'C:/Users/egann/Desktop/Aureococcus_strains/Comparing Proteomes/uniprot_blastp'

#open the blast table and just pull top hit
#if the percent identity line[2] > 30 and the evalue line[10] < 1e-50
top_hits_table = []
seen = set()

with open(os.path.join(directory,'Aa_strains_uniprot.txt'),'r') as f:
	for line in f:
		if line.strip().split('\t') not in seen:
			seen.add(line.strip().split('\t'))
			data = line.strip().split('\t')
			if int(data[2]) >= 30:
				if int(data[10]) <= 1e-50:
					top_hits_table.append(data)

for x in top_hits_table:
	print(x)
