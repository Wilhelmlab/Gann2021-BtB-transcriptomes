#imports
import os
import csv 

#directory 
directory = 'C:/Users/egann/Desktop/BtB Transcriptome/Add Backs/output'


#get files in the directory
files = os.listdir(directory)

with open(os.path.join(directory,files[0]),'r') as f:
	#pull the names of all the mappings 
	CDS_names = []
	for line in f:
		CDS_names.append(line.strip().split('\t')[0])

del(CDS_names[0])

#get the lengths of each library 
library_lengths = []

with open('help.txt','r') as f:
	for line in f:
		library_lengths.append(line.strip().split('\t'))

out = []
#for each file: open it, get the MCP 
for file in files:
	out_line = [file.split('_')[0]]

	lib_length = 0
	file_table = []
	#open file and add get the table 
	with open(os.path.join(directory,file),'r') as f:
		for line in f:
			file_table.append(line.strip().split('\t'))
	#delete the first row 
	del(file_table[0])
	#get the length of the library 
	for library in library_lengths:
		if file.startswith(library[0]):
			lib_length = float(library[1])
	#get the major capsid protein
	#and normalize it by lib size and length of CDS
	for line in file_table:
		if '_YP_009052173.1_' in line[0]:
			MCP_length = float(line[1])
			MCP_read_counts = float(line[2])
			MCP_normalized = MCP_read_counts/(MCP_length*lib_length)

	for name in CDS_names:
		for line in file_table:
			if name == line[0]:
				line_length = float(line[1])
				line_read_counts = float(line[2])
				line_normalized = line_read_counts/(line_length*lib_length)
				to_add = line_normalized/MCP_normalized
				out_line.append(to_add)

	print(file.split('_')[0])
	out.append(out_line)

CDS_names.insert(0,'names')
out.insert(0,CDS_names)

with open('temp.txt','w') as o:
	writer = csv.writer(o,delimiter='\t')
	writer.writerows(zip(*out))

with open(os.path.join(directory,'normalized_MCP_all.txt'),'w') as o:
	with open('temp.txt','r') as f:
		for line in f:
			if line != '\n':
				o.write(line)