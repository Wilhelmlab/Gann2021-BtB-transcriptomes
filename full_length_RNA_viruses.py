###Find full length RNA viruses

#imports
import os
import csv

#directory 
directory = 'C:/Users/egann/Desktop/BtB Transcriptome/RNA Viruses Contigs/Determining_full_length'

#open BLASTx table that the blast was all Riboviria contigs 
#queried against all Riboviria proteins from RefSeq in NCBI 
#add to a list of lists
blastx_table = []

with open(os.path.join(directory,'riboviria_contigs_riboviria_proteins.txt'),'r') as f: 
	for line in f: 
		blastx_table.append(line.strip().split('\t'))


#open the full contigs list and get the names
full_contigs = dict()
key = ""

with open(os.path.join(directory,'riboviria_contigs.fasta'),'r') as f: 
	for line in f:
		if line.startswith('>'):
			key = line.strip().replace('>',"")
			full_contigs[key] = ""
		else:
			full_contigs[key] += line.strip()


#get the list of contig names 
contig_names = set()

for line in blastx_table:
	if line[0] not in contig_names:
		contig_names.add(line[0])



#for each contig subset the blast x table 
out_dict = dict()
contigs_with_multiple = []

for contig in contig_names:

	final_fragments = []
	#subsetted blastx table 
	contig_table = []
	locations_of_alignment = []

	for line in blastx_table:
		if line[0] == contig:
			contig_table.append(line)
			location_pair = [int(line[6]),int(line[7])]
			location_pair.sort()
			location_pair.append(location_pair[1] - location_pair[0])
			if location_pair not in locations_of_alignment:
				locations_of_alignment.append(location_pair)
	
	#sort locations_of_alignment by size 
	locations_of_alignment.sort(key=lambda x: x[2],reverse=True)
	
	#find all alignment locations not found within the largest fragment
	large_fragment = locations_of_alignment[0]

	not_in = []

	for alignment_pair in locations_of_alignment:
		if alignment_pair[0] < large_fragment[0]:
			if alignment_pair[1] < large_fragment[0]:
				not_in.append(alignment_pair)
		elif alignment_pair[0] > large_fragment[1]:
			if alignment_pair[1] > large_fragment[1]:
				not_in.append(alignment_pair)

	final_fragments.append(large_fragment)

	if len(not_in) != 0:

		large_fragment = not_in[0]

		not_in_take_2 = []

		for alignment_pair in not_in:
			if alignment_pair[0] < large_fragment[0]:
				if alignment_pair[1] < large_fragment[0]:
					not_in_take_2.append(alignment_pair)
			elif alignment_pair[0] > large_fragment[1]:
				if alignment_pair[1] > large_fragment[1]:
					not_in_take_2.append(alignment_pair)

		final_fragments.append(large_fragment)
		
		if len(not_in_take_2) != 0:
			
			large_fragment = not_in_take_2[0]

			not_in_take_3 = []

			for alignment_pair in not_in_take_2:
				if alignment_pair[0] < large_fragment[0]:
					if alignment_pair[1] < large_fragment[0]:
						not_in_take_3.append(alignment_pair)
				elif alignment_pair[0] > large_fragment[1]:
					if alignment_pair[1] > large_fragment[1]:
						not_in_take_3.append(alignment_pair)

			final_fragments.append(large_fragment)

	#if there are multiple fragments within the contig
	#pull those to be then run through the PFAM

	#In total there are 214, with multiple PFAMs, pull those 
	if len(final_fragments) != 1:
		for full_contig in full_contigs:
			if full_contig == contig:
				sequence = full_contigs[full_contig]
		contigs_with_multiple.append([contig,len(final_fragments),len(sequence)])
		for fragment in final_fragments:
			frag_name = '>' + contig + '_' + str(fragment[0]) + '_' + str(fragment[1])
			out_dict[frag_name] = sequence[fragment[0]:fragment[1]]

#write out dict to an outfile 
with open(os.path.join(directory,'Riboviria_contig_fragments_have_multiple.fasta'),'w') as o:
	for contig in out_dict:
		o.write(contig)
		o.write('\n')
		o.write(out_dict[contig])
		o.write('\n')

#write contigs with multiple to an out as well 
with open(os.path.join(directory,'contigs_with_multiple_hits.txt'),'w') as o:
	writer = csv.writer(o,delimiter='\t')
	writer.writerows(contigs_with_multiple)


#get top hit from BLASTx table of Riboviria fragments 
seen = set() 
subject_top_hits = []
query_subject = []

with open(os.path.join(directory,'riboviria_fragments_blastx_table_plus_frame.txt'),'r') as f:
	for line in f: 
		if line.strip().split('\t')[0] not in seen:
			seen.add(line.strip().split('\t')[0])
			subject_top_hits.append(line.strip().split('\t')[1])
			query_subject.append([line.strip().split('\t')[0],line.strip().split('\t')[1].replace('|','').replace("ref",'')])

subject_top_hits = list(set(subject_top_hits))

#pull out the top hits to then go and search for 

#open the Riboviria refseq file 
riboviria_proteins = dict()
key = ""

with open(os.path.join(directory,'Riboviria_proteins.fasta'),'r') as f:
	for line in f:
		if line.startswith('>'): 
			key = line.strip().split(' ')[0]
			riboviria_proteins[key] = ""
		else: 
			riboviria_proteins[key] += line.strip()

with open(os.path.join(directory,'Riboviria_references_to_get_PFAM.fasta'),'w') as o:
	for hit in subject_top_hits:
		count = 0
		for protein in riboviria_proteins:
			if protein.strip('>') == hit.replace('|','').replace("ref",''):
				o.write(protein)
				o.write('\n')
				o.write(riboviria_proteins[protein])
				o.write('\n')


#Open the PFAM domain search table and get the PFAM accession and description 
PFAM_domain_table = []

with open(os.path.join(directory,'PFAM_domain_subject_hits.txt'),'r') as f:
	for line in f:
		PFAM_domain_table.append(line.strip().split('\t'))

#pull out the information from the PFAM table and append it to query_subject to write to an out file 
out_table = []

for query in query_subject:
	out_line = query
	count = 0
	for line in PFAM_domain_table:
		if query[1] == line[0]:
			out_line.append(line[4])
			out_line.append(line[8])
	out_table.append(out_line)

#write to an outfile
with open(os.path.join(directory,'riboviria_contig_fragments_PFAMs.txt'),'w') as o:
	writer = csv.writer(o,delimiter='\t')
	writer.writerows(out_table)


#get final list of contigs with multiple hits, that include a structural protein, an RDRP, or something else 

#first open the list of contigs with multiple hits 
contigs_with_multiple = []

with open(os.path.join(directory,'contigs_with_multiple_hits.txt'),'r') as f:
	for line in f:
		if line != '\n':
			contigs_with_multiple.append(line.strip().split('\t'))

#open the contigs_fragments_with PFAM file
fragments_with_PFAMS = []

with open(os.path.join(directory,'contig_fragments_with_PFAM.txt'),'r') as f:
	for line in f:
		if line != '\n':
			fragments_with_PFAMS.append(line.strip().split('\t'))


#groups 
Groups = ['Structural','RDRP','Peptidase','Helicase','no_PFAM']

#make a final out table 
out_table = []

for contig in contigs_with_multiple:
	#out line
	out_line = contig

	#get subset of fragment tables
	begins_with = contig[0] + '_'
	contig_subset = []
	for line in fragments_with_PFAMS:
		if line[0].startswith(begins_with):
			contig_subset.append(line)
	
	for group in Groups: 
		out_group = []
		for line in contig_subset:
			if line[4] == group:
				to_add = line[0] + ' (' + line[2] + ')'
				out_group.append(to_add)
		out_line.append(str(out_group).strip('[').strip(']').replace("'",""))


	for data in full_contigs:
		if contig[0] == data:
			out_line.append(full_contigs[data])

	out_table.append(out_line)

Groups.append('seqeunce')
Groups.insert(0,'length')
Groups.insert(0,'alignment_fragment_count')
Groups.insert(0,'contig')

out_table.insert(0,Groups)

#write to an out file 
with open('temp.txt','w') as o:
	writer = csv.writer(o,delimiter='\t')
	writer.writerows(out_table)

with open(os.path.join(directory,'contigs_with_multiple_proteins_by_type.txt'),'w') as o:
	with open('temp.txt','r') as f:
		for line in f:
			if line != '\n':
				o.write(line)

os.remove('temp.txt')