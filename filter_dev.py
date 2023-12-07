from Bio import SeqIO
from itertools import compress
import re
import os
import argparse

'''
Script usage
------------------------------------------------------------------------------------------
Takes a directory with fasta files isoforms in which to be filtered, 
yeilds the output directory with filtered files.

Example of usage:
cmd>>> python3 filter.py -i ./input_data_dir/ -o ./output_data_dir/
'''

# Defining command line arguments
# Directories names must be provided with slash (/) at the end
parser = argparse.ArgumentParser(description="Protein isoforms filtration")
parser.add_argument("-i", "--indir", help="Directory with unfiltered fasta", required=True)
parser.add_argument("-o", "--outdir", help="Directory for output fasta", required=True)
args = parser.parse_args()

data_dir = args.indir
output_dir = args.outdir

# Create output directory if it doesn't exist
if not os.path.isdir(output_dir):
	os.mkdir(output_dir)

# Get all files from input data directory
files = os.listdir(data_dir)

for file in files:
	print(f'Processing file: {file}') # notify about processing start
	records = list(SeqIO.parse(open(data_dir + file), 'fasta'))
	regex = re.compile("\[gene=[\w]*\]")  # \w: [a-zA-Z0-9_]
	
	gene_names = [] # a list to track if we already met specific gene name (length equals to the number of unique gene names)
	names = [] # a list to build filter list further (length equals to the number of records)

	for record in records:
		gene = regex.search(record.description)
		
		if gene is not None: 
			gene = gene.group()
			gene_name = gene[gene.find('=')+1: -1] # extract name of gene
			
			# If the gene name has not been encountered yet
			if gene_name not in gene_names: 
				gene_names.append(gene_name)
				names.append(gene_name)
				
			# If gene name had already been encountered:
			else: 
				length1 = len(records[names.index(gene_name)].seq) # check the length of the previous isoform
				
				# If the length of current isoform if bigger than it of the previous one:
				if len(record.seq) > length1:
					names[names.index(gene_name)] = 'drop this' # delete gene name from the first isoform position
					names.append(gene_name) # add gene name for current isoform position
				else:
					names.append('drop this')
		
		else: # in case record does not have a gene name
			names.append('no name')

	filter_list = [x != 'drop this' for x in names] # create a boolean list, where False is for record to be dropped and True is for record to be left
	filtered_records = list(compress(records, filter_list)) # drop certain records

	new_file = output_dir + file[:-4] + '.filtered.fa' # create output filename
	with open(new_file, 'w') as nf:
		SeqIO.write(filtered_records, nf, 'fasta')
	print(f'Processing done: {file}\n...\n') # notify about processing completed
	
