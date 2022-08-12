#Script written by ACL to grab NTD seqs from R2G files based on a set of post-Guidance files
#and categorizes them in a couple of different ways. 09/16/21


import os, sys
from Bio import SeqIO


R2G_handle = 'NTDfiles300' #absolute or relative path to folder of NTD R2G files
clades_of_interest = ['Ba', 'Za', 'Op', 'EE', 'Ex', 'Am', 'Sr', 'Pl'] #2-5 dig. code
postguidance_handle = '408_meiosis_postguidance95gapTrimmed' #absolute or relative path to folder of postguidance files

#grab all of the R2G seqs we may want
seq_ids = { record.description : str(record.seq) for file in os.listdir(R2G_handle) for record in SeqIO.parse(R2G_handle + '/' + file, 'fasta') if file.endswith('.fasta') and file.startswith(clades_of_interest[0])}
for i in range(1,len(clades_of_interest)) :
	seq_ids.update({ record.description : str(record.seq) for file in os.listdir(R2G_handle) for record in SeqIO.parse(R2G_handle + '/' + file, 'fasta') if file.endswith('.fasta') and file.startswith(clades_of_interest[i])})

#Make output directory for the seqs filtered and categorized by postguidance file
if(not os.path.isdir('Filtered_NTDs_by_PGFile')):
	os.mkdir('Filtered_NTDs_by_PGFile')

#A file to hold ALL filtered NTD seqs
all_decontam = open('all_filtered_NTD_seqs.fasta', 'w')

#Write out all of the seqs of interest from each postguidance file
pg_seqs = []
for file in os.listdir(postguidance_handle):
	if(file.endswith('.fas')):
		with open('Filtered_NTDs_by_PGFile/' + file + '_filtered_NTD.fasta', 'w') as o:
			for rec in SeqIO.parse(postguidance_handle + '/' + file, 'fasta'):
				for clade in clades_of_interest:
					if(rec.description.startswith(clade)):
						try:
							o.write('>' + rec.description + '\n' + seq_ids[rec.description] + '\n\n')
							all_decontam.write('>' + rec.description + '\n' + seq_ids[rec.description] + '\n\n')
							pg_seqs.append(rec.description)
						except KeyError:
							print('\nUh oh... ' + rec.description + ' could not be found in the R2G files! Skipping this sequence.\n')

all_decontam.close()

#Make output directory for the seqs filtered and categorized by taxon
if(not os.path.isdir('Filtered_NTDs_by_Taxon')):
	os.mkdir('Filtered_NTDs_by_Taxon')

#Write out all of the seqs of interest from each postguidance file, but now separating out by taxon
for taxon in list(dict.fromkeys([rec[:10] for rec in pg_seqs])):
	print(taxon)
	with open('Filtered_NTDs_by_Taxon/' + taxon + '_filtered_NTD.fasta', 'w') as o:
		for rec in pg_seqs:
			if(rec[:10] == taxon):
				o.write('>' + rec + '\n' + seq_ids[rec] + '\n\n')