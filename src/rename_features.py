from prepare_trace_dataset import load_dataset


def rename_feature_name(df, old_features_name, new_feature_name):

    new_df = df.copy()
    new_df.columns = df.columns.astype(str).str.replace(old_features_name, new_feature_name)

    return new_df


def rename_dataset_features(job_name=''):

    results_dir = '../results/%s' % job_name
    f = '%s/df_complete_dataset.csv' % results_dir,

    df_complete, num_features = load_dataset(f, remove_non_expressed=False)

    df_complete = rename_feature_name(df_complete, 'Net_embedding_', 'net_embed_')

    df_complete = rename_feature_name(df_complete, '.4wpc_Development', '_dev_expression_4wpc')
    df_complete = rename_feature_name(df_complete, '.5wpc_Development', '_dev_expression_5wpc')
    df_complete = rename_feature_name(df_complete, '.6wpc_Development', '_dev_expression_6wpc')
    df_complete = rename_feature_name(df_complete, '.7wpc_Development', '_dev_expression_7wpc')
    df_complete = rename_feature_name(df_complete, '.8wpc_Development', '_dev_expression_8wpc')
    df_complete = rename_feature_name(df_complete, '.9wpc_Development', '_dev_expression_9wpc')
    df_complete = rename_feature_name(df_complete, '.10wpc_Development', '_dev_expression_10wpc')
    df_complete = rename_feature_name(df_complete, '.11wpc_Development', '_dev_expression_11wpc')
    df_complete = rename_feature_name(df_complete, '.12wpc_Development', '_dev_expression_12wpc')
    df_complete = rename_feature_name(df_complete, '.13wpc_Development', '_dev_expression_13wpc')
    df_complete = rename_feature_name(df_complete, '.14wpc_Development', '_dev_expression_14wpc')
    df_complete = rename_feature_name(df_complete, '.15wpc_Development', '_dev_expression_15wpc')
    df_complete = rename_feature_name(df_complete, '.16wpc_Development', '_dev_expression_16wpc')
    df_complete = rename_feature_name(df_complete, '.17wpc_Development', '_dev_expression_17wpc')
    df_complete = rename_feature_name(df_complete, '.18wpc_Development', '_dev_expression_18wpc')
    df_complete = rename_feature_name(df_complete, '.19wpc_Development', '_dev_expression_19wpc')
    df_complete = rename_feature_name(df_complete, '.20wpc_Development', '_dev_expression_20wpc')
    df_complete = rename_feature_name(df_complete, '.newborn_Development', '_dev_expression_newborn')
    df_complete = rename_feature_name(df_complete, '.infant_Development', '_dev_expression_infant')
    df_complete = rename_feature_name(df_complete, '.toddler_Development', '_dev_expression_toddler')
    df_complete = rename_feature_name(df_complete, '.school_Development', '_dev_expression_school')
    df_complete = rename_feature_name(df_complete, '.teenager_Development', '_dev_expression_teenager')
    df_complete = rename_feature_name(df_complete, '.youngAdult_Development', '_dev_expression_youngAdult')
    df_complete = rename_feature_name(df_complete, '.youngMidAge_Development', '_dev_expression_youngMidAge')
    df_complete = rename_feature_name(df_complete, '.olderMidAge_Development', '_dev_expression_olderMidAge')
    df_complete = rename_feature_name(df_complete, '.senior_Development', '_dev_expression_senior')

    df_complete = rename_feature_name(df_complete, 'Development_CV', 'dev_expression_var')

    df_complete = rename_feature_name(df_complete, '_LCV', '_expression_var')

    df_complete = rename_feature_name(df_complete, '_mean_tipa_pathways', '_diff_proc_act_mean')
    df_complete = rename_feature_name(df_complete, '_median_tipa_pathways', '_diff_proc_act_med')
    df_complete = rename_feature_name(df_complete, '_min_tipa_pathways', '_diff_proc_act_min')
    df_complete = rename_feature_name(df_complete, '_max_tipa_pathways', '_diff_proc_act_max')
    df_complete = rename_feature_name(df_complete, '_paths_num_pathways', '_num_processes')

    df_complete = rename_feature_name(df_complete, '_paralogs_ratio_highest_identity', '_paralog_ratio_single')

    df_complete = rename_feature_name(df_complete, '_egene', '_eQTL')

    df_complete = rename_feature_name(df_complete, '_diff_net_mean', '_diff_PPI_mean')
    df_complete = rename_feature_name(df_complete, '_diff_net_med', '_diff_PPI_med')
    df_complete = rename_feature_name(df_complete, '_diff_net_min', '_diff_PPI_min')
    df_complete = rename_feature_name(df_complete, '_diff_net_max', '_diff_PPI_max')

    df_complete = rename_feature_name(df_complete, '_num_specific_interactions_dif_mean', '_num_specific_PPIs_dif_mean')
    df_complete = rename_feature_name(df_complete, '_num_specific_interactions_dif_median', '_num_specific_PPIs_dif_med')
    df_complete = rename_feature_name(df_complete, '_num_specific_interactions', '_num_specific_PPIs')

    df_complete = rename_feature_name(df_complete, '_num_specific_interactions', '_num_specific_PPIs')
    df_complete = rename_feature_name(df_complete, '_num_specific_interactions_dif_median', '_num_specific_PPIs_dif_med')
    df_complete = rename_feature_name(df_complete, '_num_specific_interactions_dif_mean', '_num_specific_PPIs_dif_mean')

    df_complete = rename_feature_name(df_complete, '_num_elevated_interactors_dif_mean', '_num_pref_PPIs_dif_mean')
    df_complete = rename_feature_name(df_complete, '_num_elevated_interactors_dif_median', '_num_pref_PPIs_dif_med')
    df_complete = rename_feature_name(df_complete, '_num_elevated_interactors', '_num_pref_PPIs')

    df_complete = rename_feature_name(df_complete, '_num_interactors_dif_mean', '_num_PPIs_dif_mean')
    df_complete = rename_feature_name(df_complete, '_num_interactors_dif_med', '_num_PPIs_dif_med')
    df_complete = rename_feature_name(df_complete, '_num_interactors', '_num_PPIs')

    df_complete = rename_feature_name(df_complete, 'preferential_expression', 'pref_expression')

    # df_complete.to_csv(f.split('.csv')[0] + '_cor_features_names.csv')
    df_complete.to_csv(f)

    # load_dataset(f.split('.csv')[0] + '_cor_features_names.csv', remove_non_expressed=False)

    return
