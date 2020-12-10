#imports 
import os 
import csv

#make a smaller blast table of only the top hits and the accession 
accessions_seen = set()
top_hits_table = []
viral_proteins = set()

with open('Virus_contigs_blastx_against_Viruses.txt','r') as f: 
	for line in f: 
		if line.strip().split('\t')[0] not in accessions_seen:
			top_hits_table.append([line.strip().split('\t')[0],line.strip().split('\t')[1].split('.')[0]])
			accessions_seen.add(line.strip().split('\t')[0])
			viral_proteins.add(line.strip().split('\t')[1].split('.')[0])

#open family to group file 
fam_to_group = []

with open('family_to_group.txt','r') as f: 
	for line in f: 
		fam_to_group.append(line.strip().split('\t'))


#open viral proteins downloaded from NCBI
viral_proteins_table = []

with open('viral_proteins.csv','r') as f:
	for line in f:
		if line.strip().split(',')[0] in viral_proteins:
			viral_proteins_table.append(line.strip().split(','))

#write an out file 


with open('virus_contigs_NCBI_data.txt','w') as o:
	out_table = []

	for line_th in top_hits_table:
		out_line = []
		out_line.append(line_th[0])
		out_line.append(line_th[1])
		for line_vp in viral_proteins_table:
			if line_th[1] == line_vp[0]:
				out_line.append(line_vp[3])
				out_line.append(line_vp[5])
				for line_ftg in fam_to_group:
					if line_vp[5] == line_ftg[0]:
						out_line.append(line_ftg[1])						

		out_table.append(out_line)

	writer = csv.writer(o,delimiter='\t')
	writer.writerows(out_table)

