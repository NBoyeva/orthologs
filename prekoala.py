from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import argparse

of_path = False
threads = 4

parser = argparse.ArgumentParser(description="OrthoFinder Output Parsing")
parser.add_argument("-w", "--workdir", help="Working directory for the project", required=False)
parser.add_argument("-p", "--orth", help="OrthoFinder path", required=False)
parser.add_argument("-i", "--input", help="Inpur directory for OrthoFinder", required=False)
parser.add_argument("-t", "--threads", help="Number of threads", required=False)
parser.add_argument("-r", "--results", help="Directory for OrthoFinder results being parsed", required=False)
args = parser.parse_args()

work_dir = args.workdir
of_path = args.orth
in_dir = args.input
threads = args.threads
results_dir = args.results

if work_dir:
    os.chdir(work_dir)
#os.chdir('/home/timurk/bio/orthomaster/data/Boyeva_all_files/abs')
# os.chdir("/home/timurk/bsu/laboratory_analysis/proteoms/")

# Running OrthoFinder
if of_path:
    os.system(of_path + ' -f ' + in_dir + ' -t ' + threads)

# import orthologs summary from OrthoFinder output
if len(os.listdir(os.path.join('.','OrthoFinder'))) == 1: 
    path_to_ortho_finder_data = os.path.join('.', 
                                             'OrthoFinder', 
                                             os.listdir(os.path.join('.','OrthoFinder'))[0], 
                                             'Orthogroups', 
                                             'Orthogroups.tsv')
elif results_dir:
    path_to_ortho_finder_data = os.path.join('.', 
                                             'OrthoFinder', 
                                             results_dir, 
                                             'Orthogroups', 
                                             'Orthogroups.tsv')
else:
    available_results_dirs = os.listdir(os.path.join('.','OrthoFinder'))
    os.system(f'echo "Specify results directory. Available options: {available_results_dirs}"')


#path_to_ortho_finder_data = './OrthoFinder/Results_Dec14/Orthogroups/Orthogroups.tsv'
# path_to_ortho_finder_data = '/home/timurk/bsu/laboratory_analysis/proteoms/oresults/Orthogroups/Orthogroups.tsv'

if not os.path.exists(script_output := os.path.join('.', 'orthologs')):
    os.mkdir(os.path.join('.', 'orthologs'))

orthogroups = pd.read_csv(path_to_ortho_finder_data, sep='\t')

# a dataframe with number of orthologs in one species in given orthogroup
counts = orthogroups.copy()

# if species has no orthologs in the orthogroup, OrthoFinder returns NaN
# here we replace NaN by 0
counts = counts.replace(np.nan, 0)

# If a species had paralogs, the value in the dataframe
# was a string separated into identifiers by comma.
# So here we make a list from it
counts = counts.map(lambda x: len(x.split(', ')) if type(x) == str else x)

# make a subset of initial dataset based on counts of orthologs
# (they should be equal to 1 since we need 1:1 orthologs)
idents = orthogroups.iloc[:, 1:][counts.iloc[:, 1:] == 1]
idents = pd.concat([orthogroups['Orthogroup'], idents], axis=1)

# since we need orthologs represented in all organisms,
# we drop orthogroups with NaN values
idents_all_sp = idents.dropna()

# removing '.filtered' from headers (when we delete isoforms,
# script filter.py adds this to new filename)
idents_all_sp.columns = map(lambda x: x.replace('.filtered', ''),
                            idents_all_sp.columns)

# To run BlastKOALA we need to provide a fasta file with sequences
# to be functionally annotated. The function of genes within orthogroup
# we consider the same. So to get a list of KO-numbers (we will need
# them further during annotation process) for all orthogroups we can
# extract one sequence from each orthogroup by means of extraction of all
# sequences from an arbitrary species listed in idents_all_sp dataframe.
# In this case arbitrary species is Corynebacterium_amycolatum.

idents_1_sp = idents_all_sp.iloc[:, 1].tolist()

path_to_filtered_genome = idents_all_sp.iloc[:, 1].name

all_recs_1_sp = list(SeqIO.parse(os.path.join('.',
                                              'filtered',
                                              'path_to_filtered_genome' + '.filtered.fa'), 
                                              'fasta'))

recs_1_sp = []

for rec in all_recs_1_sp:
    if rec.id in idents_1_sp:
        recs_1_sp.append(rec)

# Writing list of sequences of one of the specious to fasta
SeqIO.write(recs_1_sp, os.path.join(script_output, 'fasta_1_sp.fa'), 'fasta')
idents_all_sp.to_csv(os.path.join(script_output, 'orthologs.csv'))