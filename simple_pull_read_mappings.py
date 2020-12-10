#combine read mapping files 

#imports 
import os 
import csv

#get the path to the location you want 
path = 'C:/Users/egann/Desktop/BtB Transcriptome/Read Mappings/Read Mappings to FEMS 2016 Amplicons (0.9,0.97) ignored'

#get the names of all files in that location 
files = os.listdir(path)

#get the names of all the CDS mapped by opening the first file 

with open(os.path.join(path,files[0]),'r') as f:
	#make rows from in file 
	order = []

	for line in f:
		if len(line) > 5:
			order.append(line.strip().split('\t')[0])

	
#pull all data from all files in the correct order 

out = []

for file in files: 

	out_by_file = [file.split('_')[0]]

	with open(os.path.join(path,file),'r') as f:
		
		for line in f: 
				out_by_file.append(line.strip().split('\t')[2])
		

	out.append(out_by_file)


	print(file)

order.insert(0,"Transcriptome")
out.insert(0,order)

with open('FEMS_MCP_read_mappings_all.txt','w') as o:
	writer = csv.writer(o,delimiter='\t')
	writer.writerows(zip(*out))