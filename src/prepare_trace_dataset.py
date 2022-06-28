import os

import numpy as np
import pandas as pd
from fancyimpute import IterativeImputer
from sklearn import preprocessing

#########################################################


job_name = 'expression_and_ppi'
results_dir = '../results/%s' % job_name
f = "%s/df_complete_dataset.csv" % results_dir

if not os.path.exists(results_dir):
    os.makedirs(results_dir)

remove_features = False
keep_features = False
features_to_keep = ['_expression', '_interactors', '_interactions']  # choose here one at a time

only_disease_causing_genes = False

merge_brain_labels = True
merge_heart_labels = True
merge_all_causal = False
unclassified_is_false = True

network_na_is_zero = True

iterative_impute = True
na_is_zero = not iterative_impute
impute_lcv = not iterative_impute
impute_paralogs = not iterative_impute

boolean_is_binary = True
add_pca = False
power_transform_expression_data = True
scale_data = True

prepare_data = True

#########################################################


def load_dataset(file_path='', remove_non_expressed=True):
    df = pd.read_csv(file_path, index_col=0)

    num_expression_tissues = count_num_expression_tissues(df)
    print('num_expression_tissues:', num_expression_tissues)

    num_lables = sum('_causal' in s for s in list(df.columns.values))
    print('num_lables :', num_lables)

    num_features = len(list(df.columns.values)) - num_lables
    print('num_features: ', num_features)

    if remove_non_expressed:
        df = remove_non_expressed_genes(df, num_expression_tissues)
        print('num of expressed genes in dataset:', df.shape)

    print('num genes in dataset:', df.shape)

    return df, num_features


def load_adapted_dataset(results_dir):
    df = pd.read_csv("%s/df_complete_dataset_ready_adapted_no_missing_values.csv" % results_dir, index_col=0)

    num_lables = sum('_causal' in s for s in list(df.columns.values))
    print('num_lables :', num_lables)

    num_features = len(list(df.columns.values)) - num_lables
    print('num_features: ', num_features)

    return df, num_features


def remove_non_expressed_genes(input_df, num_expression_tissues):
    new_df = input_df
    new_df = new_df.dropna(subset=list(input_df.columns.values)[:num_expression_tissues], how='all')

    return new_df


def remove_non_causal_genes(input_df, num_features):
    new_df = input_df.dropna(subset=[list(input_df.columns.values)[num_features:]], inplace=False, how='all')

    return new_df


def merge_heart_causal(input_df):
    new_df = input_df

    new_df.loc[(new_df['heart_atrial_appendage_causal'] == 1) | (new_df['heart_left_ventricle_causal'] == 1),
               'heart_causal'] = 1

    return new_df


def merge_brain_causal(input_df):
    new_df = input_df

    new_df.loc[(new_df['brain-not_specific_causal'] == 1) | (new_df['brain-0_causal'] == 1) |
               (new_df['brain-1_causal'] == 1) | (new_df['brain-2_causal'] == 1),
               'brain_causal'] = 1

    return new_df


def remove_features_containing(input_df, feature, total_num_features):

    feature_to_remove = feature
    cols = list(input_df.columns.values)
    cols_to_remove = [s for s in cols if ((feature_to_remove in s) and ('causal' not in s))]
    new_df = input_df.drop(cols_to_remove, axis=1)
    return_num_features = total_num_features - len(cols_to_remove)

    return new_df, return_num_features


def keep_features_containing(input_df, feature_strings_to_keep, total_num_features):

    cols = list(input_df.columns.values)
    feature_strings_to_keep.append('causal')
    print('feature_strings_to_keep:', feature_strings_to_keep)

    cols_to_remove = set()

    for col in cols:
        print(col)
        if all(feature_to_keep not in col for feature_to_keep in feature_strings_to_keep):
            cols_to_remove.add(col)

    cols_to_remove_list = list(cols_to_remove)
    print(cols_to_remove_list)
    new_df = input_df.drop(cols_to_remove_list, axis=1)

    return_num_features = total_num_features - len(cols_to_remove_list)

    return new_df, return_num_features


def power_transform_expression(input_df):
    cols = list(input_df.columns.values)
    expression_cols = [s for s in cols if 'expression' in s]
    print(expression_cols)

    transformer = preprocessing.PowerTransformer(method='yeo-johnson')
    input_df.loc[:, expression_cols] = transformer.fit_transform(input_df.loc[:, expression_cols])

    new_df = input_df

    return new_df


def turn_unclassified_to_false(input_df, num_features):
    new_df = input_df
    new_df[list(input_df.columns.values)[num_features:]] = \
        input_df[list(input_df.columns.values)[num_features:]].fillna(False)

    return new_df


def lcv_na_to_median(input_df):
    cols = [s for s in input_df.columns.values if 'LCV' in s]
    new_df = input_df
    new_df[cols] = new_df[cols].fillna(new_df[cols].median())

    return new_df


def paralogs_na_to_median(input_df):
    cols = [s for s in input_df.columns.values if 'paralogs' in s]
    new_df = input_df
    new_df[cols] = new_df[cols].fillna(new_df[cols].median())

    return new_df


def na_to_zero(input_df, cols):
    new_df = input_df
    new_df[cols] = input_df[cols].fillna(0)

    return new_df


def boolean_to_binary(input_df):
    new_df = input_df

    new_df = new_df.applymap(lambda x: 1 if x == True else x)
    new_df = new_df.applymap(lambda x: 0 if x == False else x)

    return new_df


def count_num_expression_tissues(input_df):
    cols = list(input_df.columns.values)
    count = 0
    for ele in cols:
        if 'expression' in ele and 'preferential' not in ele:
            count = count + 1
    return count


def run_prepare_dataset(job_name='', features_to_keep=None):

    results_dir = '../results/%s' % job_name
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    print('preparing data')

    if prepare_data:

        df_complete, num_features = load_dataset(f)

        if only_disease_causing_genes:
            df_complete = remove_non_causal_genes(df_complete, num_features)

        if merge_heart_labels:
            df_complete = merge_heart_causal(df_complete)

        if merge_brain_labels:
            df_complete = merge_brain_causal(df_complete)

        if remove_features:
            feature_to_remove = 'preferential'  # choose here one at a time
            df_complete, num_features = remove_features_containing(df_complete, feature_to_remove, num_features)

        if keep_features:
            df_complete, num_features = keep_features_containing(df_complete, features_to_keep, num_features)

            print(df_complete)
            # num_features = len([s for s in list(df_complete.columns.values) if (feature_to_keep in s)])

            print('new num_features:', num_features)

        if unclassified_is_false:
            df_complete = turn_unclassified_to_false(df_complete, num_features)

        if boolean_is_binary:
            df_complete = boolean_to_binary(df_complete)

        if impute_lcv:
            df_complete = lcv_na_to_median(df_complete)

        if impute_paralogs:
            df_complete = paralogs_na_to_median(df_complete)
            # df_complete = paralogs_ratio_to_1(df_complete)

        if na_is_zero:
            df_complete = na_to_zero(df_complete, list(df_complete.columns.values)[:num_features])

        if network_na_is_zero:

            cols = list(df_complete.columns.values)
            cols_to_impute = [s for s in cols if 'num' in s]

            df_complete = na_to_zero(df_complete, cols_to_impute)

        if power_transform_expression_data:
            df_complete = power_transform_expression(df_complete)

        if scale_data:
            scaler = preprocessing.MinMaxScaler(feature_range=(-1, 1))
            df_complete.iloc[:, :num_features] = scaler.fit_transform(df_complete.iloc[:, :num_features])

        if iterative_impute:
            df_complete.iloc[:, :num_features] = IterativeImputer(max_iter=10,
                                                                  tol=0.01,
                                                                  verbose=2,
                                                                  initial_strategy="median",
                                                                  n_nearest_features=100,
                                                                  # sample_posterior=True,
                                                                  random_state=0).fit_transform(df_complete.iloc[:,
                                                                                                :num_features])

        # test dataset for NaNs and non-finits
        mat = df_complete.iloc[:, :num_features]
        print(np.any(np.isnan(mat)))
        print(np.all(np.isfinite(mat)))

        df_complete.to_csv('%s/df_complete_dataset_ready_adapted_no_missing_values.csv' % results_dir)

    else:
        df_complete, num_features = load_adapted_dataset(results_dir)

        if keep_features:

            df_complete, num_features = keep_features_containing(df_complete, features_to_keep, num_features)

            print(df_complete)

            print('new num_features:', num_features)

            df_complete.to_csv('%s/df_complete_dataset_ready_adapted_no_missing_values.csv' % results_dir)

    return
