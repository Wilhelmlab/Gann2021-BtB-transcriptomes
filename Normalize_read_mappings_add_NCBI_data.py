#imports
import os 
import csv 

#open library size file and write to a list
lib_size = []

with open('Library_size.txt','r') as f:
	for line in f:
		line_data = line.strip().split('\t')
		no_commas = line_data[1].replace(',','')
		lib_size.append([line_data[0],float(no_commas)])

#get list of all files from specific directory 
path = 'C:/Users/egann/Desktop/BtB Transcriptome-7.29.2020/working/Read Mappings cd-hit clusters from that transcriptome'

transcriptome_files = os.listdir(path)


#out file path directory 
out_path = 'C:/Users/egann/Desktop/BtB Transcriptome-7.29.2020/working/out'

#for each read mapping file
#first get length of the library from lib_size
#then divide each contig's # of reads mapped by the lib_size
#and the contig length 

for file in transcriptome_files:
	#get the length of the library from lib_size
	trimmed_reads_size = 0
	for lib in lib_size:
		if file.split('_')[0] == lib[0]:
			trimmed_reads_size = lib[1]

	#open each read mappings file make readlines list
	unedited_readmappings = []

	with open(os.path.join(path,file),'r') as f:
		for line in f:
			unedited_readmappings.append(line.strip().split('\t'))

	#remove headers
	del unedited_readmappings[0]

	#make a searchable list that just has the contig name 
	#in the first position without the fraction and _mapping
	#and the normalized reads mapped 
	new_list = []

	for line in unedited_readmappings:
		contig = line[0].split('-')[0]
		consensus_length = float(line[1])
		normalizing_factor = consensus_length*trimmed_reads_size
		normalized_reads = float(line[2])/normalizing_factor
		new_list.append([contig,normalized_reads])

	#get information from the virus_contigs_NCBI_data file 
	#to get the organism, family, and group 
	out_list = []
	for data in new_list:

		with open('virus_contigs_NCBI_data.txt','r') as f:

			for line in f:
				if line.strip().split('\t')[0] == data[0]:
					out_list.append([data[0],data[1],line.strip().split('\t')[2],line.strip().split('\t')[3],line.strip().split('\t')[4]])


	#write to an out file

	out_name = file.split('_')[0] + '_normalized_reads_NCBI_data.txt'

	with open(os.path.join(out_path,out_name),'w') as o:
		writer = csv.writer(o,delimiter='\t')
		writer.writerows(out_list)