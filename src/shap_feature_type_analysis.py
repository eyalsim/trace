import pandas as pd


def infer_feature_type(name):

    name = name.lower()

    new_name = ''

    if 'diff_proc_act' in name:
        new_name = 'Tissue Differential Process Activity'
    elif 'diff_proc.act' in name:
        new_name = 'Tissue Differential Process Activity'
    elif 'diff_PPI' in name:
        new_name = 'Tissue Differential PPIs'
    elif 'net_embed' in name:
        new_name = 'Tissue Network Embedding'
    elif 'pref_expression' in name:
        new_name = 'Tissue Preferential Expression'
    elif 'num_PPIs' in name:
        new_name = 'Tissue PPIs'
    elif 'num_pref_PPIs' in name:
        new_name = 'Tissue PPIs'
    elif 'specific_PPIs' in name:
        new_name = 'Tissue PPIs'
    elif 'dev_expression_var' in name:
        new_name = 'Tissue Expression Variability'
    elif 'dev_expression' in name:
        new_name = 'Tissue Expression (development) '
    elif 'expression_var' in name:
        new_name = 'Tissue Expression Variability'
    elif 'eQTL' in name:
        new_name = 'Tissue eQTL'
    elif 'paralog' in name:
        new_name = 'Tissue Paralogs Relationships'
    elif 'expression' in name:
        new_name = 'Tissue Expression'

    return new_name


def infer_origin_tissue(name):

    name = name.lower()

    new_name = ''

    if 'heart' in name:
        new_name = 'Heart'
    elif 'brain' in name:
        new_name = 'Brain'
    elif 'kidney' in name:
        new_name = 'Kidney'
    elif 'esophagus' in name:
        new_name = 'Esophagus'
    elif 'spleen' in name:
        new_name = 'Spleen'
    elif 'pituitary' in name:
        new_name = 'pituitary'
    elif 'adrenal gland' in name:
        new_name = 'Adrenal Gland'
    elif 'skin' in name:
        new_name = 'Skin'
    elif 'uterus' in name:
        new_name = 'Uterus'
    elif 'artery' in name:
        new_name = 'Artery'
    elif 'muscle' in name:
        new_name = 'Muscle'
    elif 'stomach' in name:
        new_name = 'Stomach'
    elif 'small intestine' in name:
        new_name = 'Small Intestine'
    elif 'pancreas' in name:
        new_name = 'Pancreas'
    elif 'thyroid' in name:
        new_name = 'Thyroid'
    elif 'vagina' in name:
        new_name = 'Vagina'
    elif 'ovary' in name:
        new_name = 'Ovary'
    elif 'adipose' in name:
        new_name = 'Adipose'
    elif 'liver' in name:
        new_name = 'Liver'
    elif 'cervix' in name:
        new_name = 'Cervix'
    elif ('transformed lymphocytes' in name) or ('transformed_lymphocytes' in name):
        new_name = 'Transformed Lymphocytes'
    elif 'nerve' in name:
        new_name = 'Nerve'
    elif 'testis' in name:
        new_name = 'Testis'
    elif 'blood' in name:
        new_name = 'Blood'
    elif 'minor salivary gland' in name:
        new_name = 'Minor Salivary Gland'
    elif 'colon' in name:
        new_name = 'Colon'
    elif 'lung' in name:
        new_name = 'Lung'
    elif 'fallopian tube' in name:
        new_name = 'Fallopian Tube'
    elif 'fibroblasts' in name:
        new_name = 'Fibroblasts'
    elif 'breast' in name:
        new_name = 'Breast'
    elif 'bladder' in name:
        new_name = 'Bladder'
    elif 'prostate' in name:
        new_name = 'Prostate'

    return new_name


#####################################################################################################################


def run_feature_type_analysis(job_name=''):

    results_dir = '../results/%s' % job_name

    tissues = ['whole_blood', 'whole_brain', 'heart', 'liver',
               'muscle_skeletal', 'nerve_tibial', 'skin', 'testis']

    feature_types = ['Tissue Expression',
                     'Tissue Expression (development)',
                     'Tissue Preferential Expression',

                     'Tissue Expression Variability',

                     'Tissue Network Embedding',
                     'Tissue PPIs',
                     'Tissue Differential PPIs',

                     'Tissue Differential Process Activity',

                     'Tissue eQTL',
                     'Tissue Paralogs Relationships'
                     ]

    feature_tissues = ['Blood', 'Brain', 'Heart',
                       'Liver',  'Muscle', 'Nerve',
                       'Skin', 'Testis']

    df_for_feature_plot = pd.DataFrame(index=feature_types, columns=tissues)
    df_for_tissue_origin_plot = pd.DataFrame(index=feature_tissues, columns=tissues)

    for tissue in tissues:
        print(tissue)

        df = pd.read_csv("%s/shap_importance_df_%s_causal__infer_tissues.csv" % (results_dir, tissue), index_col=0)

        for row_name in df.index.values:

            feature_name = df.loc[row_name, 'column_name']

            feature_type = infer_feature_type(feature_name)

            origin_tissue = infer_origin_tissue(feature_name)

            df.loc[row_name, 'feature_type'] = feature_type
            df.loc[row_name, 'origin_tissue'] = origin_tissue

        df.to_csv("%s/shap_importance_df_%s_causal_general_feature_analysis.csv" % (results_dir, tissue))

        shap_sum = df.loc[:, 'shap_importance'].sum(axis=0)

        for feature in feature_types:
            feature_sum = df.loc[df['feature_type'] == feature, 'shap_importance'].sum()
            feature_normalized_sum = feature_sum / shap_sum
            print(feature, feature_normalized_sum)

            df_for_feature_plot.loc[feature, tissue] = feature_normalized_sum

        for feature_tissue in feature_tissues:
            feature_sum = df.loc[df['origin_tissue'] == feature_tissue, 'shap_importance'].sum()
            feature_normalized_sum = feature_sum / shap_sum
            print(feature_tissue, feature_normalized_sum)

            df_for_tissue_origin_plot.loc[feature_tissue, tissue] = feature_normalized_sum

    df_for_feature_plot.rename(columns={'heart': 'Heart',
                                        'muscle_skeletal': 'Muscle',
                                        'skin': 'Skin',
                                        'liver': 'Liver',
                                        'testis': 'Testis',
                                        'whole_blood': 'Blood',
                                        'nerve_tibial': 'Nerve',
                                        'whole_brain': 'Brain',
                                        'cortex': 'Cortex',
                                        'cerebellum': 'Cerebellum'},
                               inplace=True
                               )

    df_for_feature_plot.to_csv("%s/df_for_feature_plot.csv" % results_dir)


    df_for_tissue_origin_plot.rename(columns={'heart': 'Heart',
                                              'muscle_skeletal': 'Muscle',
                                              'skin': 'Skin',
                                              'liver': 'Liver',
                                              'testis': 'Testis',
                                              'whole_blood': 'Blood',
                                              'nerve_tibial': 'Nerve',
                                              'whole_brain': 'Brain',
                                              'cortex': 'Cortex',
                                              'cerebellum': 'Cerebellum'},
                                     inplace=True
                                     )

    df_for_tissue_origin_plot.to_csv("%s/df_for_tissue_origin_plot.csv" % results_dir)

    return

