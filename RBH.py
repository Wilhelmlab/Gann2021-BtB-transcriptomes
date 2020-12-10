#open files 
Poly_a_blast_table = []
Ribozero_blast_table = []

with open('PolyA_Ribozero.txt','r') as f:
	already_seen = set()
	for line in f:
		if line.strip().split('\t')[0] not in already_seen:
			Poly_a_blast_table.append([line.strip().split('\t')[0],line.strip().split('\t')[1]])
			already_seen.add(line.strip().split('\t')[0])

with open('Ribozero_PolyA.txt','r') as f:
	already_seen = set()
	for line in f:
		if line.strip().split('\t')[0] not in already_seen:
			Ribozero_blast_table.append([line.strip().split('\t')[0],line.strip().split('\t')[1]])
			already_seen.add(line.strip().split('\t')[0])

with open('Ribozero_PolyA_clustered_RBH.txt','w') as o:

	for pair_pa in Poly_a_blast_table:
		for pair_rz in Ribozero_blast_table:
			if pair_pa[0] == pair_rz[1]: 
				if pair_pa[1] == pair_rz[0]:
					o.write(pair_pa[0])
					o.write('\n')
					o.write(pair_pa[1])
					o.write('\n')