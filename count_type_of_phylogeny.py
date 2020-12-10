families = []
groups = []
organisms = []

with open('virus_contigs_NCBI_data.txt','r') as f: 
	for line in f:
		families.append(line.strip().split('\t')[3])
		groups.append(line.strip().split('\t')[4])
		organisms.append(line.strip().split('\t')[2])

families = list(set(families))
groups = list(set(groups))
organisms = list(set(organisms))


out_counts_family = []
out_counts_groups = []
out_counts_organisms = []

for data in families:
	data_count = 0
	with open('virus_contigs_NCBI_data.txt','r') as f: 
		for line in f:
			if line.strip().split('\t')[3] == data:
				data_count = data_count + 1
	out_counts_family.append([data,data_count])

for data in groups:
	data_count = 0
	with open('virus_contigs_NCBI_data.txt','r') as f: 
		for line in f:
			if line.strip().split('\t')[4] == data:
				data_count = data_count + 1
	out_counts_groups.append([data,data_count])

for data in organisms:
	data_count = 0
	with open('virus_contigs_NCBI_data.txt','r') as f: 
		for line in f:
			if line.strip().split('\t')[2] == data:
				data_count = data_count + 1
	out_counts_organisms.append([data,data_count])


import csv 

with open('out_counts_family.txt','w') as o:
	writer = csv.writer(o,delimiter='\t')
	writer.writerows(out_counts_family)

with open('out_counts_groups.txt','w') as o:
	writer = csv.writer(o,delimiter='\t')
	writer.writerows(out_counts_groups)

with open('out_counts_organisms.txt','w') as o:
	writer = csv.writer(o,delimiter='\t')
	writer.writerows(out_counts_organisms)