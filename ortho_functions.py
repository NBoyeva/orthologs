from Bio.AlignIO.PhylipIO import SequentialPhylipWriter
from Bio.KEGG import REST as kegg
from Bio import SeqIO
from Bio.Phylo.PAML import yn00
from Bio import AlignIO
import re
import pandas as pd
import os
from Bio.Phylo.PAML._paml import PamlError
import matplotlib.pyplot as plt
import seaborn as sns
import Bio
from Bio.Phylo.TreeConstruction import dis
import shutil


def find_substring(s):
    pattern = r'WP_.*?\.1'
    matches = re.findall(pattern, s)
    return matches


def get_sp_by_prot_ident(og, ident, idents_all_sp):
    """
    Input:
    og -- orthogroup
    ident -- protein sequence identifiers

    Output:
    Returns a species to which the sequence with identifiers given belongs.

    Example:
    >> get_sp_by_prot_ident('OG0000310', 'lcl|NZ_CP066023.1_prot_WP_197915016.1_491')
    >>'Corynebacterium_amycolatum'
    """

    # retrieving the list of protein identifiers from the orthogroup
    og_idents = list(idents_all_sp[idents_all_sp['Orthogroup'] == og].iloc[0])[1:]

    # retrieving the list of species from the orthogroup
    species = list(idents_all_sp[idents_all_sp['Orthogroup'] == og])[1:]

    # creating a dictionary with identifiers:species pairs
    ident_dict = dict(zip(og_idents, species))

    return ident_dict[ident]


def _get_kegg(kegg_id):
    kegg_output = kegg.kegg_get(kegg_id).read()
    results = {}
    for line in kegg_output.split('\n'):
        splits = line.split()
        if not line.startswith(' '):
            if len(splits) > 0:
                key = splits[0]
                value = ' '.join(splits[1:])
                results[key] = value
        else:
            results[key] += ' '.join(splits)
    return pd.DataFrame(results, index=[kegg_id])


def extract_maps(text):
    return [match.group(0) for match in re.finditer(r'(map\d+.*?)(?=map|$)', text)]


def get_kegg_info(kegg_ids, _get_kegg_v):
    if isinstance(kegg_ids, str):
        kegg_ids = [kegg_ids]
    return pd.concat(_get_kegg_v(kegg_ids), sort=False)


def get_desc_by_og(og, ko, all_kegg):
    """
    Input:
    og -- orthogroup id

    Output:
    Returns the description of orthogroup sequences function.

    Example:
    >> get_desc_by_og('OG0000310')
    >> 'phosphorrosslerite carboxylase (GTP) [EC:4.1.1.32]'
    """

    kon = ko[ko['Orthogroup'] == og][1].tolist()[0]
    desc = all_kegg[all_kegg['Unnamed: 0'] == kon]['NAME'].tolist()[0]
    return desc


def get_fasta_by_group(df, pathway, ko, idents_all_sp):
    """
    Input:
    df -- dataframe consisting info about pathway
    pathway -- short pathway name (pathway id)

    Output:
    Writes CDS sequences of orthologs within pathway groups for each orthogroup separately.
    """

    print('Starting creating fasta-files for orthogroups in pathway: ', pathway)

    output_folder = os.path.join('.', 'pathways')
    pathway_folder = os.path.join(output_folder, pathway)
    if not os.path.exists(pathway_folder):
        os.mkdir(pathway_folder)

    dna_seqs_path = os.path.join('.', 'filtered_dna')  ##### Replace on filtered dna

    # iterating over ko-numbers of certain metabolic pathway
    for kon in df['Unnamed: 0']:

        dna_records = []

        # inferring orthogroup number from ko-number
        og = list(ko[ko[1] == kon]['Orthogroup'])[0]

        # retrieving the list of protein identifiers from the orthogroup
        og_idents = list(idents_all_sp[idents_all_sp['Orthogroup'] == og].iloc[0])[2:]  ### Why 2????

        # the list of species containing these proteins respectively
        species = list(idents_all_sp[idents_all_sp['Orthogroup'] == og])[2:]  ### Why 2????

        ### writing file with DNA sequences ###

        for n, ident in enumerate(og_idents):  # iterating over protein identifiers within the orthogroup

            dna_ident = ident.replace('prot', 'cds')  # converting protein identifiers to CDS one

            sp_dna_filename = species[n] + '.fa'  # retrieving filename for the species that contains the protein
            sp_dna_path = os.path.join(dna_seqs_path, sp_dna_filename)  # creating path for this file
            sp_dna = list(SeqIO.parse(sp_dna_path, 'fasta'))  # parsing this file

            for record in sp_dna:  # iterating over DNA sequences within species
                if dna_ident in record.id:
                    dna_records.append(record)  # adding CDSs with specified identifiers

        if not os.path.exists(os.path.join(pathway_folder, og)):
            os.makedirs(os.path.join(pathway_folder, og))

        with open(os.path.join(pathway_folder, og, og + '.fa'), 'w') as f:
            SeqIO.write(dna_records, f, 'fasta')
            print('Done for ', og)

    print('Finished creating fasta-files for orthogroups in pathway: ', pathway)


def pairwise_paml(og, pathway, path_to_ogs, idents_all_sp, ext_sp='Mycobacterium_leprae'):
    """
    Input:
    og -- orthogroup, within which we will create pairwise fasta
    pathway -- the metabolic pathway within which we are working
    ext_sp -- species name, with which all the orthologs from other species will be aligned

    Workflow:
    Creates pairwise protein and DNA fasta within orthogroup.
    Makes protein alignment.
    Converts it first to codon-based DNA alignment in fasta format, further to phyla format.
    Runs PAML on codon-based DNA alignment.

    Output:
    Returns PAML program output (including dN/dS)
    """

    print(f'Working with {og} in pathway: {pathway}')

    # path to protein fasta file with records of given orthogroup
    og_prot_path = os.path.join(path_to_ogs, og + '.fa')

    # path to DNA fasta file with records of given orthogroup
    og_dna_path = os.path.join('.', 'pathways', pathway, og, og + '.fa')

    # making directories for pairwise protein and DNA fasta files
    # (pairwise file contains a protein from a certain species from given orthogroup
    # and protein from external species in this orthogroup)
    pairwise_prot_path = os.path.join('.', 'pathways', pathway, og, 'pairwise_prot')
    pairwise_dna_path = os.path.join('.', 'pathways', pathway, og, 'pairwise_dna')

    for folder in [pairwise_prot_path, pairwise_dna_path]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    ### Working with protein sequences ###

    # parse protein sequences from the whole orthogroup
    og_prot = list(SeqIO.parse(og_prot_path, 'fasta'))

    ### Retrieving external species protein record ###

    for record in og_prot:  # iterating over protein sequences from the whole orthogroup
        record_sp = get_sp_by_prot_ident(og, record.id, idents_all_sp)
        if record_sp == ext_sp:
            record_ext = record  # save the sequence from the external species to record_ext
            break

    ### Writing pairwise protein fasta files ###

    for record in og_prot:

        records_to_write = [record_ext]

        record_sp = get_sp_by_prot_ident(og, record.id, idents_all_sp)
        if record_sp != ext_sp:
            records_to_write.append(record)
            file_path = os.path.join(pairwise_prot_path, record_sp + '.fa')
            SeqIO.write(records_to_write, file_path, 'fasta')

    ### Working with DNA sequences ###

    og_dna = list(SeqIO.parse(og_dna_path, 'fasta'))
    ext_sp_dna_id = record_ext.id.replace('prot', 'cds')

    ### Retrieving external species protein record ###
    record_ext_dna = ''
    for record in og_dna:
        if record.id == ext_sp_dna_id:
            record_ext_dna = record
            break

    ### Writing pairwise DNA fasta files ###

    for record in og_dna:
        records_to_write = [record_ext_dna]
        prot_like_ident = record.id.replace('cds', 'prot')
        record_sp = get_sp_by_prot_ident(og, prot_like_ident, idents_all_sp)
        if record_sp != ext_sp:
            records_to_write.append(record)
            file_path = os.path.join(pairwise_dna_path, record_sp + '.fa')
            SeqIO.write(records_to_write, file_path, 'fasta')

            ### Creating all necessary directories ###

    pairwise_alignments_path = os.path.join('.', 'pathways', pathway, og, 'pairwise_alignments')
    pal2nal_output_path = os.path.join('.', 'pathways', pathway, og, 'pal2nal')
    phylip_alignments_path = os.path.join('.', 'pathways', pathway, og, 'phylip_alignments')
    paml_output_path = os.path.join('.', 'pathways', pathway, og, 'paml_output')

    for folder in [pairwise_alignments_path,
                   pal2nal_output_path,
                   phylip_alignments_path,
                   paml_output_path]:
        if not os.path.exists(folder):
            os.makedirs(folder)

    ### Creating pairwise alignments ###

    for prot in os.listdir(pairwise_prot_path):  # prot -- protein fasta file

        prot_basename = prot.split('.')[0]  # prot_basename -- the name of species
        prot_path = os.path.join(pairwise_prot_path, prot)
        alignment_path = os.path.join(pairwise_alignments_path, prot_basename + '.clw')
        os.system('clustalw align -infile=' + prot_path + ' -outfile=' + alignment_path + ' -QUIET')

        ### Accessing DNA sequences of orthogroup to get full identifiers ###

        og_dna = list(SeqIO.parse(og_dna_path, 'fasta'))
        full_ids = []
        for record in og_dna:
            full_ids.append(record.id)

        ### Rewrite alignment in fasta with full identifiers as in DNA sequences ###

        alignment = AlignIO.read(alignment_path, "clustal")
        for seq in alignment:
            cut_id = seq.id
            cut_id = cut_id.replace('prot', 'cds')
            for full_id in full_ids:
                if cut_id in full_id:
                    seq.id = full_id

        new_alignment_path = os.path.join(pairwise_alignments_path, prot_basename + '.fa')
        AlignIO.write(alignment, new_alignment_path, 'fasta')

        ### Running pal2nal ###

        dna_path = os.path.join(pairwise_dna_path, prot_basename + '.fa')

        pal2nal_path = os.path.join(pal2nal_output_path, prot_basename + '.fa')

        os.system(
            os.path.join('.', 'pal2nal.v14', 'pal2nal.pl') + ' ' + new_alignment_path + ' '
            + dna_path + ' -codontable 11 -output fasta > ' + pal2nal_path)

        ### Converting alignment to phylip format ###'
        aln = AlignIO.read(pal2nal_path, 'fasta')

        for seq in aln:
            seq.id = seq.id[:15]

        phylip_path = os.path.join(phylip_alignments_path, prot_basename + '.phy')

        with open(phylip_path, 'w') as output_handle:
            SequentialPhylipWriter(output_handle).write_alignment(aln, id_width=30)

        ### Running PAML ###

        paml_output_sp_path = os.path.join(paml_output_path, prot_basename)

        if not os.path.exists(paml_output_sp_path):
            os.makedirs(paml_output_sp_path)

        o_file = os.path.join(paml_output_sp_path, 'paml.out')
        yn = yn00.Yn00(alignment=phylip_path,
                       working_dir=paml_output_sp_path,
                       out_file=o_file)

        try:
            yn.run(verbose=True)
        except PamlError:
            pass

        print('Finished PAML for ', prot_basename)


def get_dN_dS(pathway, og, species):
    """
    Input:
    pathway -- pathway -- pathway id (in form of map_____ (5 numbers))
    og -- orthogroup
    species -- species of ortholog for which dN/dS is calculated

    Output:
    For a specific species in orthogroup within pathway returns dN/dS value.

    """

    dN_path = os.path.join('.', 'pathways', pathway, og, 'paml_output', species, '2YN.dN')
    dS_path = os.path.join('.', 'pathways', pathway, og, 'paml_output', species, '2YN.dS')

    try:
        dN = pd.read_csv(dN_path).iloc[1, 0].split(' ')[-1]
    except FileNotFoundError:
        print(f"Can't find {dN_path}")
        return -10

    try:
        dS = pd.read_csv(dS_path).iloc[1, 0].split(' ')[-1]
    except FileNotFoundError:
        print(f"Can't find {dS_path}")
        return -10

    dN_dS = float(dN) / float(dS)
    return dN_dS


def get_matrix(pathway_full, pathway, idents_all_sp, pathway_groups, ko, all_kegg):
    """
    Input:
    pathway_full -- full name of the pathway as it is mentioned in KEGG
    pathway -- pathway id (in form of map_____ (5 numbers))

    Workflow:
    Calculates a matrix od dN/dS values for all species in orthogroups within pathway.

    Output:
    Returns matrix, list of species and annotations for orthogroups.
    """

    species_list = list(map(lambda x: x.split('.')[0], idents_all_sp.columns.tolist()[2:]))
    # print(all_species)
    species_ext = idents_all_sp.columns[-1]
    species_list.remove(species_ext)

    pathway_df = pathway_groups.get_group(pathway_full)

    ogs = []
    names = []
    for kon in pathway_df['Unnamed: 0']:  # iterating over ko-numbers of certain metabolic pathway
        og = list(ko[ko[1] == kon]['Orthogroup'])[0]  # inferring orthogroup number from ko-number
        ogs.append(og)
        names.append(get_desc_by_og(og, ko, all_kegg))

    matrix_list = []

    for og in ogs:
        og_row = []
        for species in species_list:
            og_row.append(get_dN_dS(pathway, og, species))
        matrix_list.append(og_row)

    matrix = pd.DataFrame(matrix_list)
    matrix.columns = species_list
    matrix.index = names

    return matrix, species_list, names


def plot_heatmap(pathway_full, pathway, idents_all_sp, pathway_groups, ko, all_kegg):
    """
    Input:
    pathway_full -- full name of the metabolic pathway as it is mentioned in KEGG
    pathway -- pathway id (in form of map_____ (5 numbers))

    Output:
    Saves and shows a heatmap showing dN/dS in orthogroups of a given pathway.
    """

    matrix, species_list, names = get_matrix(pathway_full, pathway, idents_all_sp, pathway_groups, ko, all_kegg)

    sns.set()

    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300

    fig, ax = plt.subplots()
    fig.set_size_inches(10, 5)
    cbar_ax = fig.add_axes([.02, .1, .04, .8])

    ax = sns.heatmap(matrix, ax=ax, cbar_ax=cbar_ax)

    ax.xaxis.tick_top()  # x-axis on top
    ax.xaxis.set_label_position('top')
    ax.yaxis.tick_right()  # y-axis on right
    ax.yaxis.set_label_position('right')

    xlab = ['sp' + str(i + 1) for i in range(matrix.shape[1])]
    print(len(names), matrix.shape)
    ax.set_xticklabels(xlab)
    ax.set_yticklabels(names, rotation=0)

    ax.set_title(pathway_full[9:])

    plt.savefig(f'pathways/{pathway}.png', bbox_inches='tight')


def plot_heatmap_2(pathway_full, pathway, idents_all_sp, pathway_groups, ko, all_kegg):
    """
    Input:
    pathway_full -- full name of the metabolic pathway as it is mentioned in KEGG
    pathway -- pathway id (in form of map_____ (5 numbers))

    Output:
    Saves and shows a heatmap showing dN/dS in orthogroups of a given pathway.
    """

    matrix, species_list, names = get_matrix(pathway_full, pathway, idents_all_sp, pathway_groups, ko, all_kegg)

    sns.set()

    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300

    if len(names) > 27:
        # Split the matrix and names into three
        third = len(names) // 3
        matrices = [matrix[:third], matrix[third:2 * third], matrix[2 * third:]]
        names_list = [names[:third], names[third:2 * third], names[2 * third:]]

        fig, axs = plt.subplots(3, 1, figsize=(10, 15), sharex=True)
        cbar_ax = fig.add_axes([.02, .1, .04, .8])

        for i in range(3):
            ax = sns.heatmap(matrices[i], ax=axs[i], cbar_ax=cbar_ax)
            ax.yaxis.tick_right()  # y-axis on right
            ax.yaxis.set_label_position('right')
            ax.set_yticklabels(names_list[i], rotation=0)

    else:
        fig, ax = plt.subplots(figsize=(10, 5))
        cbar_ax = fig.add_axes([.02, .1, .04, .8])

        ax = sns.heatmap(matrix, ax=ax, cbar_ax=cbar_ax)

        ax.xaxis.tick_top()  # x-axis on top
        ax.xaxis.set_label_position('top')
        ax.yaxis.tick_right()  # y-axis on right
        ax.yaxis.set_label_position('right')

        ax.set_yticklabels(names, rotation=0)

    xlab = ['sp' + str(i + 1) for i in range(matrix.shape[1])]
    ax.set_xticklabels(xlab)

    ax.set_title(pathway_full[9:])

    plt.savefig(f'pathways/{pathway}.png', bbox_inches='tight')
