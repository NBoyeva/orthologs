from ortho_functions import (_get_kegg, get_kegg_info, extract_maps, get_fasta_by_group, pairwise_paml, get_matrix,
                             get_dN_dS, plot_heatmap, plot_heatmap_2)
import pandas as pd
import os
import numpy as np

# os.chdir("/home/timurk/bsu/laboratory_analysis/proteoms/")
os.chdir('/home/timurk/bio/orthomaster/data/Boyeva_all_files/abs')

script_output = os.path.join('.', 'orthologs')
idents_all_sp = pd.read_csv(os.path.join(script_output, 'orthologs.csv'))
path_to_ko = './orthologs/user_ko.txt'

# import BlastKOALA output
ko = pd.read_csv(path_to_ko, sep='\t', header=None)
# add orthogroups to the dataframe
ko = pd.concat([ko, pd.Series(list(idents_all_sp['Orthogroup']), name='Orthogroup')], axis=1)
# drop orthogroups without KO-numbers
ko = ko.dropna()

# Takes time to run, so we saved the result to csv.
# _get_kegg_v = np.vectorize(_get_kegg)
# all_kegg = get_kegg_info(list(ko[1]), _get_kegg_v)
# all_kegg.to_csv(os.path.join(script_output, 'all_kegg.csv'))

all_kegg = pd.read_csv(os.path.join('.', 'orthologs', 'all_kegg.csv'))
all_kegg = all_kegg.dropna(subset='PATHWAY')  # drop not annotated records
all_kegg['PATHWAY'] = all_kegg['PATHWAY'].apply(extract_maps)

# write the main pathway name to new column
all_kegg['TYPE'] = all_kegg['PATHWAY'].map(lambda x: x[0])

# Group the KEGG output dataframe by the main metabolic pathway
pathway_groups = all_kegg.groupby('TYPE')
# list of pathways full names
# e.g.
pathway_name_full_list = []

# we will analyse metabolic groups with more than thres orthogroups
thres = 10

for i in pathway_groups:
    if len(i[1]) >= thres:
        pathway_name_full_list.append(i[0])

# list of pathways ids (e.g.)
pathway_id_list = list(map(lambda x: x[:8], pathway_name_full_list))

# list of pathways short names
# e.g.
pathway_name_list = list(map(lambda x: x[9:], pathway_name_full_list))

# pathways with indexes 1 and 2 are removed from these lists and further pipeline since they
# caused some errors in PAML and I had no time to figure out why
pathway_name_full_list = [pathway_name_full_list[0]] + pathway_name_full_list[3:]
pathway_id_list = [pathway_id_list[0]] + pathway_id_list[3:]
pathway_name_list = [pathway_name_list[0]] + pathway_name_list[3:]

print('INFO: Pathways, that were found in orthogroups with orthologs 1:1')
for i in pathway_name_list:
    print(i)

# I cleared cell output as it takes way too much space

if not os.path.exists(os.path.join('.', 'pathways')):
    os.mkdir(os.path.join('.', 'pathways'))

for i in range(len(pathway_id_list)):

    pathway_name_full = pathway_name_full_list[i]
    pathway = pathway_id_list[i]
    pathway_name = pathway_name_list[i]

    pathway_df = pathway_groups.get_group(pathway_name_full)

    get_fasta_by_group(pathway_df, pathway, ko, idents_all_sp)

    ogs = []
    for kon in pathway_df['Unnamed: 0']:  # iterating over ko-numbers of certain metabolic pathway
        og = list(ko[ko[1] == kon]['Orthogroup'])[0]  # inferring orthogroup number from ko-number
        ogs.append(og)

    for og in ogs:
        pairwise_paml(og, pathway,
                      path_to_ogs='./OrthoFinder/Results_Dec14/Orthogroup_Sequences/',
                      # '/home/timurk/bsu/laboratory_analysis/proteoms/oresults/Orthogroup_Sequences/',
                      idents_all_sp=idents_all_sp)

for pathway_full, pathway in zip(pathway_name_full_list, pathway_id_list):
        plot_heatmap_2(pathway_full, pathway, idents_all_sp, pathway_groups, ko, all_kegg)
