#imports
import os

#open file to split and make into dictionary
fasta_dict = dict()
key = ""

with open('contigs_Viral_Ribozero_fragments_all.fasta','r') as f:
	for line in f:
		if line.startswith('>'):
			key = line.strip()
			fasta_dict[key] = ""
		else:
			fasta_dict[key] += line.strip()

#get the different names of the contigs
libraries = []

for header in fasta_dict:
	libraries.append(header.strip('>').split('_')[0])

libraries = list(set(libraries))

out_path = 'C:/Users/egann/Desktop/working/out'

for library in libraries:
	out_name = library + '_viral_contig_CDS.fasta'

	with open(os.path.join(out_path,out_name),'w') as o:
		#separate based on start of headers
		for header in fasta_dict:
			if header.strip('>').split('_')[0] == library:
				o.write(header)
				o.write('\n')
				o.write(fasta_dict[header])
				o.write('\n')