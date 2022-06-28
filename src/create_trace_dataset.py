import pandas as pd
import csv
import numpy as np
import os
import networkx as nx
import node2vec

#########################################################

WALK_LENGTH = 80
NUM_WALKS = 10

def symbol_to_ensg():
    with open("../data/hgnc_complete_set.txt", encoding="utf8") \
            as hgnc_complete_set:
        f = hgnc_complete_set
        dict = {}
        next(f)  # skip the first line (headlines)
        for line in f:
            ensg_name = ''
            for word in line.split():
                if word.startswith('ENSG'):
                    ensg_name = word
            gene_symbol = line.split("\t")[1]
            dict.update({gene_symbol: ensg_name})
    return dict


def expression_v8_cpm():
    print('compute expression features')
    exp_df = pd.read_table("../data/medians_output_v8_cpm.csv", delimiter=',', index_col=0)

    protein_coding_list = protein_coding_genes()
    # exp_df = exp_df.filter(items=protein_coding_list, axis=0)
    # exp_df = exp_df.drop_duplicates(keep=False)

    exp_df = exp_df[exp_df.index.isin(protein_coding_list)]

    # exp_df.columns = [str(col) + '_expression' for col in exp_df.columns]

    return exp_df


def df_pref_v8():
    print('compute preferential expression features')
    pref_df = \
        pd.read_table("../data/preferential_output_my_calc_V8_not_plus_one_all_18_9.csv", delimiter=',', index_col=0)

    protein_coding_list = protein_coding_genes()
    # pref_df = pref_df.filter(items=protein_coding_list, axis=0)
    # pref_df = pref_df.drop_duplicates(keep=False)
    pref_df = pref_df[pref_df.index.isin(protein_coding_list)]

    # pref_df.columns = [str(col) + '_preferential_expression' for col in pref_df.columns]

    return pref_df


def protein_coding_genes():
    df = pd.read_table('../data/biomart_protein_coding.txt',
                       delimiter=',')
    pc_list = df['Gene stable ID'].tolist()

    return pc_list


def create_df_paralogs_highest_identity():
    df1 = pd.read_table("../data/median-highest identity paralog.csv", delimiter=',')
    df1 = df1.set_index("ensgID")
    del df1['gene Description']
    del df1['paralog name']
    del df1['paralog Description']

    df2 = pd.read_table( "../data/median-rest.csv", delimiter=',')
    df2 = df2.set_index("ensgID")
    del df2['gene Description']
    del df2['paralog name']
    del df2['paralog Description']

    df = pd.concat([df1, df2])

    return df


def create_df_paralogs_all():
    df1 = pd.read_table( "../data/median-all paralogs.csv", delimiter=',')
    df1 = df1.set_index("ensgID")
    del df1['gene Description']

    return df1


def df_moran_pathways():
    df1 = pd.read_table( "../data/df_features_moran_pathways.csv", delimiter=',', index_col=0)

    return df1


def filter_expression(expression_df, th, fc_data):

    specificity = ""
    if fc_data:
        specificity = "_fc"

    df = expression_df.copy()
    local_treshold = th

    df[(df <= local_treshold)] = np.nan

    # remove genes that don't express in any tissue
    df['num_tissues'] = df.count(axis=1)
    df = df[0 != df['num_tissues']]
    del df['num_tissues']

    return df


def create_df_goldset_multi_label(merge_with_paralogs_data=False):

    df = pd.read_csv("../data/goldset_cor1.txt", delimiter='\t')

    gene_symbol_ensg_dict = symbol_to_ensg()

    df.rename(columns={'ensg_id': 'ensgID'}, inplace=True)
    df['causal_gene'] = ""
    disease_tissues = df['disease_tissues'].unique()
    new_disease_tissues = []
    for tissues in disease_tissues:
        tissue_list = tissues.split("|")
        new_disease_tissues = new_disease_tissues + tissue_list
    disease_tissues = list(set(new_disease_tissues))
    genes = df['ensgID'].unique()

    if merge_with_paralogs_data:
        df_para = pd.read_csv("../data/barshir_2017.csv",  delimiter=',')
        df_para['ensgID'] = ""

        for index_para in df_para.index:
            gene_para = df_para.at[index_para, 'causal_gene']
            if gene_para in gene_symbol_ensg_dict:
                df_para.at[index_para, 'ensgID'] = gene_symbol_ensg_dict[gene_para]

        disease_tissues_para = df_para['disease_tissue'].unique()
        genes_para = df_para['ensgID'].unique()

        disease_tissues = np.append(disease_tissues, disease_tissues_para)
        genes = np.append(genes, genes_para)

        disease_tissues = np.unique(disease_tissues)
        genes = np.unique(genes)

    new_df = pd.DataFrame(index=genes, columns=disease_tissues)
    new_df[:] = 'FALSE'

    for index in df.index:
        gene = df.at[index, 'ensgID']
        tissues = df.at[index, 'disease_tissues'].split('|')
        for tissue in tissues:
            new_df.at[gene, tissue] = 'TRUE'

    if merge_with_paralogs_data:
        for index_para in df_para.index:
            gene_para = df_para.at[index_para, 'ensgID']
            tissue_para = df_para.at[index_para, 'disease_tissue']
            new_df.at[gene_para, tissue_para] = 'TRUE'

    new_df.columns = [str(col) + '_causal' for col in new_df.columns]
    new_df = new_df.applymap(lambda x: 1 if x == 'TRUE' else x)
    new_df = new_df.applymap(lambda x: 0 if x == 'FALSE' else x)

    return new_df


def df_new_brain_labels():
    new_labels = pd.read_csv("../data/new_brain_labels.csv", index_col=0)

    return new_labels


def diff_net_interactome_dict(file_name):

    with open("../data/diffnet/" + file_name) as interactome:

        reader = csv.reader(interactome, delimiter='\t')
        dict = {}

        for line in reader:
            gene_a = line[0]
            gene_b = line[1]
            score = line[5]
            if gene_a in dict:
                dict[gene_a].append(score)
            else:
                dict.update({gene_a: [score]})
            if gene_b in dict:
                dict[gene_b].append(score)
            else:
                dict.update({gene_b: [score]})

    return dict


def my_protein_net_interactome_dict():

    with open("../data/GlobalInteractome.tsv") as interactome:

        reader = csv.reader(interactome, delimiter='\t')
        dict = {}

        for line in reader:
            gene_a = line[0]
            gene_b = line[1]
            if gene_a in dict:
                dict[gene_a].append(gene_b)
            else:
                dict.update({gene_a: [gene_b]})
            if gene_b in dict:
                dict[gene_b].append(gene_a)
            else:
                dict.update({gene_b: [gene_a]})

    for gene in dict:
        dict[gene] = set(dict[gene])
        dict[gene] = list(dict[gene])

    return dict


def my_protein_net_interactome_set():

    with open("../data/GlobalInteractome.tsv") as interactome:

        reader = csv.reader(interactome, delimiter='\t')
        interactome_set = set()

        for line in reader:
            gene_a = line[0]
            gene_b = line[1]

            ls = sorted([gene_a, gene_b])
            t = tuple(ls)
            interactome_set.add(t)

    return interactome_set


def create_interactome_dict_expression(df_filtered_expression):

    print('computing interactome dict')

    old_dict = my_protein_net_interactome_dict()
    new_dict = {}
    for gene in old_dict.keys():
        if gene in df_filtered_expression.index.values:
            new_dict.update({gene: []})
            for interactor in old_dict[gene]:
                if interactor in df_filtered_expression.index.values:
                    new_dict[gene].append(interactor)

    for gene in new_dict:
        new_dict[gene] = set(new_dict[gene])
        new_dict[gene] = list(new_dict[gene])

    return new_dict


def create_interactome_set_expression(df_filtered_expression):
    print('computing interactome set')

    int_set = my_protein_net_interactome_set()
    new_set = set()
    for t in int_set:
        if (t[0] in df_filtered_expression.index.values) and (t[1] in df_filtered_expression.index.values):
            new_set.add(t)

    return new_set


def create_tissue_interactome_dict_expression(results_dir, tissue, df_expression, interactome_dict_expression,
                                              threshold):

    tissue_interactome_dict = {}
    for gene in interactome_dict_expression.keys():
        if df_expression.loc[gene, tissue] > threshold:
            tissue_interactome_dict.update({gene: []})
            for interactor in interactome_dict_expression[gene]:
                if df_expression.loc[interactor, tissue] > threshold:
                    tissue_interactome_dict[gene].append(interactor)

    with open('%s/%s_interactome.csv' % (results_dir, tissue), 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in tissue_interactome_dict.items():
            writer.writerow([key, (",".join(value))])

    return tissue_interactome_dict


def create_tissue_interactome_dict_fc(tissue, df_expression, df_expression_fc, interactome_dict_expression,
                                      threshold, fc_threshold):
    tissue_interactome_dict = {}
    for gene in interactome_dict_expression.keys():
        if gene not in df_expression_fc.index.values:
            continue
        if df_expression.loc[gene, tissue] > threshold:
            tissue_interactome_dict.update({gene: []})
            for interactor in interactome_dict_expression[gene]:
                if interactor not in df_expression_fc.index.values:
                    continue
                if df_expression_fc.loc[interactor, tissue] > fc_threshold:
                    tissue_interactome_dict[gene].append(interactor)

    return tissue_interactome_dict


def create_tissue_specific_interactions_set(tissues_interactions_dict, tissue_enhanced_treshold):
    ts_interactions_set = set()
    for pair in tissues_interactions_dict:
        if len(tissues_interactions_dict[pair]) <= tissue_enhanced_treshold:
            ts_interactions_set.add(pair)

    return ts_interactions_set


def create_tissue_specific_interactions_dict(tissue, df_expression, interactome_dict_expression,
                                             tissue_specific_interactions_set, threshold):
    tissue_TS_interactome_dict = {}
    for gene in interactome_dict_expression.keys():
        if df_expression.loc[gene, tissue] > threshold:
            tissue_TS_interactome_dict.update({gene: []})
            for interactor in interactome_dict_expression[gene]:
                if df_expression.loc[interactor, tissue] > threshold:
                    ls = sorted([gene, interactor])
                    t = tuple(ls)
                    if t in tissue_specific_interactions_set:
                        tissue_TS_interactome_dict[gene].append(interactor)

    return tissue_TS_interactome_dict


def load_LCV():

    df = pd.read_table(
        "../data/gene_noise_500_protein_coding_0.csv",
        delimiter=',', index_col=0)

    return df


def load_developmental_Median_CPM():

    df = pd.read_table(
        "../data/Tissues_Median_CPM_Development4ML.csv",
        delimiter=',', index_col=0)

    return df


def load_developmental_CV_CPM():

    df = pd.read_table(
        "../data/Tissues_CV_CPM_Development4ML.csv",
        delimiter=',', index_col=0)

    return df


def create_df_num_interactors_in_tissues(results_dir, tissues, df_expression, interactome_dict_expression, threshold):

    df = pd.DataFrame(columns=df_expression.columns, index=df_expression.index)

    for tissue in tissues:
        tissue_interactome_dict = create_tissue_interactome_dict_expression(results_dir, tissue, df_expression,
                                                                            interactome_dict_expression, threshold)

        for gene in tissue_interactome_dict.keys():
            df.at[gene, tissue] = len(tissue_interactome_dict[gene])

    return df


def create_df_num_elevated_interactors_in_tissues(tissues, df_expression, df_expression_fc,
                                                  interactome_dict_expression, threshold, fc_threshold):
    df = pd.DataFrame(columns=df_expression.columns, index=df_expression.index)

    for tissue in tissues:
        tissue_fc_interactome_dict = \
            create_tissue_interactome_dict_fc(tissue, df_expression, df_expression_fc, interactome_dict_expression,
                                              threshold, fc_threshold)
        for gene in tissue_fc_interactome_dict.keys():
            df.loc[gene, tissue] = len(tissue_fc_interactome_dict[gene])

    return df


def create_df_num_specific_interactions_in_tissues(tissues, df_expression, interactome_dict_expression,
                                                   tissue_specific_interactions_set, threshold):
    df = pd.DataFrame(columns=df_expression.columns, index=df_expression.index)

    for tissue in tissues:
        tissue_specific_interactions_dict = \
            create_tissue_specific_interactions_dict(tissue, df_expression, interactome_dict_expression,
                                                     tissue_specific_interactions_set, threshold)
        for gene in tissue_specific_interactions_dict.keys():
            df.at[gene, tissue] = len(tissue_specific_interactions_dict[gene])

    return df


def create_tissues_interactions_dict(results_dir, tissues, df_expression, interactome_dict_expression,
                                     interactome_set_expression, threshold):

    dict = {}

    for tissue in tissues:
        tissue_interactome_dict = \
            create_tissue_interactome_dict_expression(results_dir, tissue, df_expression, interactome_dict_expression,
                                                      threshold)
        for pair in interactome_set_expression:
            if (pair[0] in tissue_interactome_dict) and (pair[1] in tissue_interactome_dict):
                if pair not in dict:
                    dict.update({pair: [tissue]})
                else:
                    dict[pair].append(tissue)

    return dict


def create_num_interactions_dif_median(tissues, old_df):

    df = old_df.copy()

    for gene in list(old_df.index.values):
        for tissue in tissues:
            lst = [float(i) for i in old_df.loc[gene]]
            df.at[gene, tissue] = \
                float(old_df.at[gene, tissue]) - np.nanmedian(lst)

    return df


def create_num_interactions_dif_mean(tissues, old_df):

    df = old_df.copy()

    for gene in list(old_df.index.values):
        for tissue in tissues:
            lst = [float(i) for i in old_df.loc[gene]]

            if all(v is np.nan for v in lst):
                df.at[gene, tissue] = np.nan
            else:
                df.at[gene, tissue] = \
                    float(old_df.at[gene, tissue]) - np.nanmean(lst)

    return df


def diff_net():
    df_diff = pd.DataFrame()
    for filename in os.listdir("../data/diffnet"):
        if filename.endswith('.tsv'):
            tissue_name = filename[:-len('.tsv')]
            print(tissue_name)

            diff_net_dict = diff_net_interactome_dict(filename)

            tissue_diff_max_dict = {}
            tissue_diff_min_dict = {}
            tissue_diff_med_dict = {}
            tissue_diff_mean_dict = {}

            for gene in diff_net_dict:
                tissue_diff_max_dict.update({gene: max(diff_net_dict[gene])})
                tissue_diff_min_dict.update({gene: min(diff_net_dict[gene])})
                tissue_diff_med_dict.update({gene: np.median([float(i) for i in diff_net_dict[gene]])})
                tissue_diff_mean_dict.update({gene: np.mean([float(i) for i in diff_net_dict[gene]])})

            tissue_df_max = pd.DataFrame.from_dict(
                tissue_diff_max_dict, orient='index', columns=[tissue_name.replace(' ', '_') + '_diff_net_max'])

            tissue_df_min = pd.DataFrame.from_dict(
                tissue_diff_min_dict, orient='index', columns=[tissue_name.replace(' ', '_') + '_diff_net_min'])

            tissue_df_med = pd.DataFrame.from_dict(
                tissue_diff_med_dict, orient='index', columns=[tissue_name.replace(' ', '_') + '_diff_net_med'])

            tissue_df_mean = pd.DataFrame.from_dict(
                tissue_diff_mean_dict, orient='index', columns=[tissue_name.replace(' ', '_') + '_diff_net_mean'])

            df_diff = pd.concat([df_diff, tissue_df_max, tissue_df_min, tissue_df_med, tissue_df_mean], axis=1, sort=False)

    return df_diff


def egenes():

    df_egenes = pd.DataFrame()
    for filename in os.listdir("../data/GTEx_Analysis_v7_eQTL"):
        if filename.endswith('.egenes.txt'):
            tissue_name = filename.split(sep='.')[0]

            feature_name = tissue_name + "_egene"

            with open("../data/GTEx_Analysis_v7_eQTL/" + filename) \
                    as tissue_egenes:

                dict = {}
                next(tissue_egenes)  # skip the first line (headlines)
                for gene in tissue_egenes:
                    ensg_name = gene.split(sep='\t')[0].split(sep='.')[0]
                    # print(ensg_name)

                    qval = gene.split(sep='\t')[28]

                    dict.update({ensg_name: float(qval)})

                tissue_df = pd.DataFrame.from_dict(
                    dict, orient='index', columns=[feature_name])

                df_egenes = pd.concat([df_egenes, tissue_df], axis=1, sort=False)

    return df_egenes


def goldset_dataset(df, features_name):
    df_for_goldset = df.copy()
    df_for_goldset.columns = [str(col) + '_%s' % features_name for col in df_for_goldset.columns]

    return df_for_goldset


def create_net_embedding(file_name=''):

    if file_name != '':
        df = pd.read_csv(file_name, skiprows=1, header=None,
                         index_col=0, sep=" ")

    else:

        df = pd.read_csv("../data/EMBEDDING.txt", skiprows=1, header=None,
                         index_col=0, sep=" ")

    num_cols = len(list(df.columns))
    rng = range(1, num_cols+1)
    col_names = ['Net_embedding_' + str(i) for i in rng]
    df.columns = col_names

    return df


def load_tissue_ppi_network(results_dir, tissue_name):
    df = pd.read_csv('%s/%s_interactome.csv' % (results_dir, tissue_name), header=None, index_col=0, sep=",")

    G = nx.Graph()

    edges = []

    for gene in df.index.values:

        if isinstance(df.loc[gene, 1], float):
            continue
        for interactor in df.loc[gene, 1].split(','):
            first = gene
            second = interactor
            edges.append((first, second))

    G.add_edges_from(edges)

    return G


def tissue_net_embedding(results_dir, tissue_name):
    df = pd.read_csv('%s/%s_EMBEDDING.txt' % (results_dir, tissue_name), skiprows=1, header=None,
                     index_col=0, sep=" ")

    num_cols = len(list(df.columns))
    rng = range(1, num_cols + 1)
    col_names = ['%s_Net_embedding_' % tissue_name + str(i) for i in rng]
    df.columns = col_names

    return df


def load_ppi_network():

    G = nx.Graph()

    f = open("../data/GlobalInteractome.tsv", "r")
    lines = f.readlines()
    edges = []
    for interaction in lines:
        first = interaction.split('\t')[0]
        second = interaction.split('\t')[1]
        edges.append((first, second))
    f.close()

    G.add_edges_from(edges)

    return G

#################################################################################################
# main


def run_create_trace_dataset(job_name='', threshold=7, fc_threshold=2, net_embedding=True):

    results_dir = '../results/%s' % job_name
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    df_expression = expression_v8_cpm()
    df_expression_fc = df_pref_v8()

    df_filtered_expression = filter_expression(df_expression, threshold, fc_data=False)
    df_filtered_expression_fc = filter_expression(df_expression_fc, fc_threshold, fc_data=True)

    df_goldset_multi_label = create_df_goldset_multi_label(merge_with_paralogs_data=True)
    df_new_brain_curation = df_new_brain_labels()
    df_goldset_multi_label = pd.concat([df_goldset_multi_label, df_new_brain_curation], axis=1, sort=True)

    interactome_dict_expression = create_interactome_dict_expression(df_filtered_expression)
    interactome_set_expression = create_interactome_set_expression(df_filtered_expression)

    tissues = list(df_expression.columns.values)
    tissue_enhanced_treshold = int(round(0.2 * len(tissues)))

    print('create_df_num_interactors_in_tissues')
    df_num_interactors_in_tissues = \
        create_df_num_interactors_in_tissues(results_dir, tissues, df_expression, interactome_dict_expression,
                                             threshold)
    df_num_elevated_interactors_in_tissues = \
        create_df_num_elevated_interactors_in_tissues(tissues, df_expression, df_expression_fc,
                                                      interactome_dict_expression, threshold, fc_threshold)

    print('create_tissues_interactions_dict')
    tissues_interactions_dict = \
        create_tissues_interactions_dict(results_dir, tissues, df_expression, interactome_dict_expression,
                                         interactome_set_expression, threshold)
    tissue_specific_interactions_set = \
        create_tissue_specific_interactions_set(tissues_interactions_dict, tissue_enhanced_treshold)

    df_num_specific_interactions_in_tissues = \
        create_df_num_specific_interactions_in_tissues(tissues, df_expression, interactome_dict_expression,
                                                       tissue_specific_interactions_set, threshold)

    df_num_interactors_in_tissues_dif_median = \
        create_num_interactions_dif_median(tissues, df_num_interactors_in_tissues)
    df_num_interactors_in_tissues_dif_mean = create_num_interactions_dif_mean(tissues, df_num_interactors_in_tissues)

    df_num_elevated_interactors_in_tissues_dif_median = \
        create_num_interactions_dif_median(tissues, df_num_elevated_interactors_in_tissues)

    df_num_elevated_interactors_in_tissues_dif_mean = \
        create_num_interactions_dif_mean(tissues, df_num_elevated_interactors_in_tissues)

    df_num_specific_interactions_in_tissues_dif_median = \
        create_num_interactions_dif_median(tissues, df_num_specific_interactions_in_tissues)

    df_num_specific_interactions_in_tissues_dif_mean = \
        create_num_interactions_dif_mean(tissues, df_num_specific_interactions_in_tissues)

    print('df_paralogs_highest_identity')
    df_paralogs_highest_identity = create_df_paralogs_highest_identity()
    df_paralogs_all = create_df_paralogs_all()

    print('df_moran_pathways')
    df_pathways = df_moran_pathways()

    print('diff_net')
    df_differential_net = diff_net()

    print('egenes')
    df_egenes = egenes()

    print('load_LCV')
    df_LCV = load_LCV()

    print('developmental')
    df_developmental1 = load_developmental_CV_CPM()
    df_developmental2 = load_developmental_Median_CPM()

    if net_embedding:
        print('load_ppi_network')
        ppi_network = load_ppi_network()

        print('Node2Vec')

        # node2vec_try = node2vec.Node2Vec(ppi_network, dimensions=64, walk_length=10, num_walks=20)

        node2vec_try = node2vec.Node2Vec(ppi_network, dimensions=64, walk_length=WALK_LENGTH, num_walks=NUM_WALKS)

        print('node2vec_try.fit')
        model = node2vec_try.fit()

        print('save_word2vec_format')
        model.wv.save_word2vec_format('%s/EMBEDDING.txt' % results_dir)

        print('net_embedding')
        df_net_embedding = create_net_embedding('%s/EMBEDDING.txt' % results_dir)

        for tissue in tissues:
            print(tissue)
            tissue_ppi_network = load_tissue_ppi_network(results_dir, tissue)

            node2vec_try_tissue = node2vec.Node2Vec(tissue_ppi_network, dimensions=64, walk_length=80, num_walks=10)

            model = node2vec_try_tissue.fit()

            model.wv.save_word2vec_format('%s/%s_EMBEDDING.txt' % (results_dir, tissue))

            df_net_embedding_tissue = tissue_net_embedding(results_dir, tissue)

            df_net_embedding = pd.concat([df_net_embedding, df_net_embedding_tissue], axis=1, sort=True)



    df_expression_for_goldset = goldset_dataset(df_expression, 'expression')
    df_expression_fc_for_goldset = goldset_dataset(df_expression_fc, 'preferential_expression')
    df_num_interactors_in_tissues_for_goldset = goldset_dataset(df_num_interactors_in_tissues, 'num_interactors')
    df_num_interactors_in_tissues_dif_median_for_goldset = \
        goldset_dataset(df_num_interactors_in_tissues_dif_median, 'num_interactors_dif_med')
    df_num_interactors_in_tissues_dif_mean_for_goldset = \
        goldset_dataset(df_num_interactors_in_tissues_dif_mean, 'num_interactors_dif_mean')
    df_num_preferential_interactors_in_tissues_for_goldset = \
        goldset_dataset(df_num_elevated_interactors_in_tissues, 'num_elevated_interactors')
    df_num_preferential_interactors_in_tissues_dif_median_for_goldset = \
        goldset_dataset(df_num_elevated_interactors_in_tissues_dif_median, 'num_elevated_interactors_dif_median')
    df_num_preferential_interactors_in_tissues_dif_mean_for_goldset = \
        goldset_dataset(df_num_elevated_interactors_in_tissues_dif_mean, 'num_elevated_interactors_dif_mean')
    df_num_specific_interactions_in_tissues_for_goldset = \
        goldset_dataset(df_num_specific_interactions_in_tissues, 'num_specific_interactions')
    df_num_specific_interactions_in_tissues_dif_median_for_goldset = \
        goldset_dataset(df_num_specific_interactions_in_tissues_dif_median, 'num_specific_interactions_dif_median')
    df_num_specific_interactions_in_tissues_dif_mean_for_goldset = \
        goldset_dataset(df_num_specific_interactions_in_tissues_dif_mean, 'num_specific_interactions_dif_mean')
    df_differential_net_for_goldset = \
        goldset_dataset(df_differential_net, '')
    df_egenes_for_goldset = \
        goldset_dataset(df_egenes, '')
    df_paralogs_highest_identity_for_goldset = \
        goldset_dataset(df_paralogs_highest_identity, 'paralogs_ratio_highest_identity')
    df_paralogs_all_for_goldset = goldset_dataset(df_paralogs_all, 'paralogs_ratio_all')
    df_pathways_for_goldset = goldset_dataset(df_pathways, 'pathways')
    df_LCV_for_goldset = goldset_dataset(df_LCV, 'LCV')
    df_developmental1_for_goldset = goldset_dataset(df_developmental1, '')
    df_developmental2_for_goldset = goldset_dataset(df_developmental2, '')


    omics_datasets_list = [df_expression_for_goldset,
                           df_expression_fc_for_goldset,
                           df_num_interactors_in_tissues_for_goldset,
                           df_num_interactors_in_tissues_dif_median_for_goldset,
                           df_num_interactors_in_tissues_dif_mean_for_goldset,
                           df_num_preferential_interactors_in_tissues_for_goldset,
                           df_num_preferential_interactors_in_tissues_dif_median_for_goldset,
                           df_num_preferential_interactors_in_tissues_dif_mean_for_goldset,
                           df_num_specific_interactions_in_tissues_for_goldset,
                           df_num_specific_interactions_in_tissues_dif_median_for_goldset,
                           df_num_specific_interactions_in_tissues_dif_mean_for_goldset,
                           df_differential_net_for_goldset,
                           df_egenes_for_goldset,
                           df_paralogs_highest_identity_for_goldset,
                           df_paralogs_all_for_goldset,
                           df_pathways_for_goldset,
                           df_LCV_for_goldset,
                           df_developmental1_for_goldset,
                           df_developmental2_for_goldset]

    df_complete_dataset = pd.concat(omics_datasets_list, axis=1, sort=True)

    if net_embedding:
        df_complete_dataset = pd.concat([df_complete_dataset, df_net_embedding], axis=1, sort=True)

    print('save df_complete_dataset')
    df_complete_dataset = pd.concat([df_complete_dataset, df_goldset_multi_label], axis=1, sort=True)
    df_complete_dataset.to_csv('%s/df_complete_dataset.csv' % results_dir)

    return
