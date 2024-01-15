from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import argparse

'''
Run OrthoFinder and parse the results:
>> python3 prekoala.py -w [project directory] -p [OrthoFinder path] -i [OrthoFinder input data directory]

Parse existing OrthoFinder results:
>> python3 prekoala.py -w [project directory] -i [OrthoFinder input data directory] -r [results directory]

All paths should be given in respect to working project directory.
'''

# By default:
of_path = False
threads = 4

# Parse arguments:
parser = argparse.ArgumentParser(description="OrthoFinder Output Parsing")
parser.add_argument("-w", "--workdir", help="Working directory for the project", required=False)
parser.add_argument("-p", "--orth", help="OrthoFinder program path", required=False)
parser.add_argument("-i", "--input", help="Inpur directory for OrthoFinder", required=False)
parser.add_argument("-t", "--threads", help="Number of threads", required=False)
parser.add_argument("-r", "--results", help="Directory for OrthoFinder results being parsed", required=False)
args = parser.parse_args()

work_dir = args.workdir     # Working directory for the project
of_path = args.orth         # OrthoFinder program path
in_dir = args.input         # Inpur directory for OrthoFinder
threads = args.threads      # Number of threads
results_dir = args.results  # Directory for OrthoFinder results being parsed (within in_dir)


### Set working directory:
if work_dir:
    os.chdir(work_dir)


### Run OrthoFinder:
if of_path:
    os.system(of_path + ' -f ' + in_dir + ' -t ' + threads)


### Import orthologs summary from OrthoFinder output -- select pathway: 
# Path to parse the only one OrthoFinder output:
if len(os.listdir(os.path.join('.', in_dir, 'OrthoFinder'))) == 1: 
    path_to_ortho_finder_data = os.path.join('.', 
                                             in_dir,
                                             'OrthoFinder', 
                                             os.listdir(os.path.join('.','OrthoFinder'))[0], 
                                             'Orthogroups', 
                                             'Orthogroups.tsv')
# Path to parse exactly given OrthoFinder output:
elif results_dir:
    path_to_ortho_finder_data = os.path.join('.', 
                                             in_dir,
                                             'OrthoFinder', 
                                             results_dir, 
                                             'Orthogroups', 
                                             'Orthogroups.tsv')
# Ask which OrthoFinder output to parse:
else:
    available_results_dirs = os.listdir(os.path.join('.', in_dir, 'OrthoFinder'))
    os.system(f'echo "Specify results directory. Available options: {available_results_dirs}"')


# Create output directory for the script:
if not os.path.exists(script_output := os.path.join('.', 'prekoala')):
    os.mkdir(os.path.join('.', 'prekoala'))



### Read OrthoFinder orthologs CSV into dataframe and process it:
orthogroups = pd.read_csv(path_to_ortho_finder_data, sep='\t')

# Create a dataframe with number of orthologs in one species in given orthogroup:
counts = orthogroups.copy()

# If species has no orthologs in the orthogroup, 
# OrthoFinder returns NaN. Replace NaN with 0:
counts = counts.replace(np.nan, 0)

# If a species had paralogs, the value in the dataframe
# was a string separated into identifiers by comma.
# Convert string with identifiers to list:
counts = counts.map(lambda x: len(x.split(', ')) if type(x) == str else x)

# Make a subset of initial dataset based on counts of orthologs
# (they should be equal to 1 since we need 1:1 orthologs)
idents = orthogroups.iloc[:, 1:][counts.iloc[:, 1:] == 1]
idents = pd.concat([orthogroups['Orthogroup'], idents], axis=1)

# Since orthologs represented in all organisms are to be analyzed,
# we drop orthogroups with NaN values:
idents_all_sp = idents.dropna()

# Remove '.filtered' from headers (when we delete isoforms,
# script filter.py adds this to new filename)
idents_all_sp.columns = map(lambda x: x.replace('.filtered', ''),
                            idents_all_sp.columns)



# To run BlastKOALA we need to provide a FASTA file with sequences
# to be functionally annotated. The function of genes within orthogroup
# we consider the same. So to get a list of KO-numbers (we will need
# them further during annotation process) for all orthogroups we can
# extract one sequence from each orthogroup by means of extraction of all
# sequences from an arbitrary species listed in idents_all_sp dataframe.
# In this case arbitrary species is Corynebacterium_amycolatum.

idents_1_sp = idents_all_sp.iloc[:, 1].tolist() # list of identificators

path_to_filtered_genome = idents_all_sp.iloc[:, 1].name # species name

# Parse all records from a species selected:
all_recs_1_sp = list(SeqIO.parse(os.path.join('.',
                                              in_dir,
                                              path_to_filtered_genome + '.filtered.fa'), 
                                              'fasta'))

# Make a list of required records:
recs_1_sp = []
for rec in all_recs_1_sp:
    if rec.id in idents_1_sp:
        recs_1_sp.append(rec)

# Write a list of sequences of one of the species to FASTA:
SeqIO.write(recs_1_sp, os.path.join(script_output, 'fasta_1_sp.fa'), 'fasta')
# Save a dataframe with 1:1 orthologs in all species to CSV:
idents_all_sp.to_csv(os.path.join(script_output, 'orthologs.csv'))