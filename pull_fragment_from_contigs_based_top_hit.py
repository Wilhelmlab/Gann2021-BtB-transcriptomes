#Date 7.24.2020
#Purpose: Pulls fragments of contigs in fasta file from alignmed 
#portion of top BLASTx hit

#imports 
import os 
from collections import defaultdict

#open out contig list file 
o = open('contigs_Viral_Ribozero_fragments_all.fasta','w')


#make smaller temp tables by pulling out all fasta headers
#using the query alignment start column[6] and query alignment
#end column[7] to pull fragment 
headers_with = list()

#small table 
accessions_seen = set() #holds lines already seen 
smaller_table = []
with open('Blastx_Viruses_Ribozero_against_Viruses.txt','r') as f: 
	#if seen accession do not write line  
	for line in f: 
		if line.strip().split('\t')[0] not in accessions_seen:
			accessions_seen.add(line.strip().split('\t')[0])
			smaller_table.append(line.strip().split('\t'))

#pull out fragment 
for line in smaller_table:
	if int(line[6]) < int(line[7]):
			#make sure fragment is larger than 150 bp 
			if int(line[7])-int(line[6]) >= 150:
				frag_name = '>' + line[0] + '-' + line[6] + '-' + line[7]
				headers_with.append(frag_name)

	if int(line[7]) < int(line[6]): 
			#make sure fragment is larger than 150 bp 
			if int(line[6])-int(line[7]) >= 150:
				frag_name = '>' + line[0] + '-' + line[7] + '-' + line[6]
				headers_with.append(frag_name)

print(len(smaller_table))
print(len(headers_with))

#open the fasta file with all contigs add to a dictionary
fasta_full = dict()
key =""

with open('Virus_top_hit_Ribo_reduced.fasta','r') as f:
		for line in f:
				if line.startswith('>'):
					key = line.strip()
					fasta_full[key] = ""
				else:
					fasta_full[key] += line.strip()

#search fasta full dict and pull out only the fragments
#write to an outfile
for fasta in fasta_full:
	for term in headers_with:
		if term.split('-')[0] == fasta:
			o.write(term)
			o.write('\n')
			o.write(fasta_full[fasta][int(term.split('-')[1]):int(term.split('-')[2])])
			o.write('\n')
quit()