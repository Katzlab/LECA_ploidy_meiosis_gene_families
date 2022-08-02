import os
import sys


output_handle = '../tip_counts.csv'


def get_args():
	use_minor_clades = False
	directory = False
	default_trees_handle = ''
	try:
		if(sys.argv[1] == '--input_file'):
			default_trees_handle = sys.argv[2]
		elif(sys.argv[1] == '--input_dir'):
			directory = True
			default_trees_handle = sys.argv[2]
		else:
			bad_script_call()
	except IndexError:
		bad_script_call()
	
	try:	
		if(sys.argv[3] == '--key'):
			key = sys.argv[4]
			taxon_list = False
		elif(sys.argv[3] == '--taxon_list'):
			key = sys.argv[4]
			taxon_list = True
		else:
			taxon_list = False
			key = 'NULL'
	except IndexError:
		taxon_list = False
		key = 'NULL'
		
	return [default_trees_handle, directory, key, taxon_list]
	
	
def get_all_taxa(file):

	taxa = []

	for line in open(file, 'r'):	
		if(line[0] == '('):
			line = line.split(':')
			for i, entry in enumerate(line):
				if(len(entry) > 25):
					taxa.append(entry.split('(')[-1].split(',')[-1][:10])
					
	return taxa


def get_tip_counts(file, dir_bool, spec_counts_per_og, codes):

	og = file.split('_bestTree.')[1][:10]
		
	if(og not in spec_counts_per_og):
		spec_counts_per_og.update({ og : codes.copy() })
		
	for line in open(file, 'r'):	
		if(line[0] == '('):
			line = line.split(':')
			for i, entry in enumerate(line):
				if(len(entry) > 25):
					entry = entry.split('(')[-1].split(',')[-1][:10]
					
					#If the following condition is not satisfied, either the taxon is not in TAXASELECTION or it does not match the key
					if(entry in spec_counts_per_og[og]):				
						spec_counts_per_og[og][entry] += 1
						
	return spec_counts_per_og
	

def main():

	dth, directory, key, taxon_list = get_args()
	
	codes = {}
	
	if(taxon_list == False):
		if(directory == False):
			taxa = get_all_taxa(dth)
		else:
			taxa = []
			for f, file in enumerate(os.listdir(dth)):
				if('ds_store' not in file.lower()):
					taxa.extend(get_all_taxa(dth + '/' + file))
		
		for taxon in taxa:
			if(key in taxon and key != 'NULL'):
				codes.update({ taxon : 0 })
			elif(key == 'NULL'):
				codes.update({ taxon : 0 })
	else:
		for line in open(key, 'r'):
			codes.update({ line.split('\n')[0].split('\t')[0].split(' ')[0] : 0 })
			
	spec_counts_per_og = {}
		
	if(directory == False):
		spec_counts_per_og = get_tip_counts(dth, directory, spec_counts_per_og, codes)
	else:
		for f, file in enumerate(os.listdir(dth)):
			if('ds_s' not in file.lower()):
				fname = dth + '/' + file
				spec_counts_per_og = get_tip_counts(fname, directory, spec_counts_per_og, codes)
								
	with open(output_handle, 'w') as o:
		o.write(',')
		for og in spec_counts_per_og:
			o.write(og + ',')
		o.write('\n')
				
		for taxon in codes:
			o.write(taxon + ',')
			for og in spec_counts_per_og:
				o.write(str(spec_counts_per_og[og][taxon]) + ',')
			o.write('\n')
				
		
main()



