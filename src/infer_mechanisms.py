import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap
from sklearn.feature_selection import SelectFromModel
from sklearn.svm import LinearSVC
from xgboost.sklearn import XGBClassifier

#####################################################################################################################

TISSUES_LIST_GENE_MODE = ['heart_causal',
                          'muscle_skeletal_causal'
                          ]

DISEASE_GENES_DICT_GENE_MODE = {
    "heart_causal": ['ENSG00000151067'],
    "muscle_skeletal_causal": ['ENSG00000198947']
}

TISSUES_LIST_TISSUE_MODE = ['heart_causal',
                            'muscle_skeletal_causal',
                            'skin_causal',
                            'whole_blood_causal',
                            'liver_causal',
                            'nerve_tibial_causal',
                            'testis_causal',
                            'whole_brain_causal'
                            ]

DISEASE_GENES_DICT_TISSUE_MODE = {}


def select_features_shap(X_train, Y, X_test, SEED):
    print("Feature selection")

    lsvc = LinearSVC(C=2, penalty="l1", dual=False, class_weight={0: 0.01, 1: 1}, tol=0.01,
                     max_iter=10000, random_state=SEED).fit(X_train, Y)

    model = SelectFromModel(lsvc, prefit=True, max_features=50)

    X_train_new = model.transform(X_train)
    if X_test is not None:
        X_test_new = model.transform(X_test.values.reshape(1, -1))
    else:
        X_test_new = None

    feature_idx = model.get_support()

    return X_train_new, X_test_new, feature_idx


def infer_models(clf=None, X_train_new=None, X_test_new=None, y_train=None, feature_names=None,
                 features_df_new=None, tissue=None, gene='', run_mode='', results_dir=''):
    print('training')
    model = clf.fit(X_train_new, y_train)

    print('shap')
    # compute SHAP values
    explainer = shap.TreeExplainer(model)

    train_shap_values = explainer.shap_values(X_train_new)
    expected_value = explainer.expected_value
    if X_test_new is not None:
        test_shap_values = explainer.shap_values(X_test_new)

    shap.summary_plot(train_shap_values, X_train_new, feature_names=feature_names,
                      max_display=20,
                      show=False)

    plt.savefig('%s/shap_summary_dotplot_%s_%s_%s.png' % (results_dir, tissue, gene, run_mode), dpi=150,
                bbox_inches='tight',
                pad_inches=0
                )
    plt.close()

    shap.summary_plot(train_shap_values, X_train_new, feature_names=feature_names,
                      max_display=20,
                      plot_type='bar',
                      color='grey',
                      show=False)

    plt.savefig('%s/shap_summary_barplot_%s_%s_%s.png' % (results_dir, tissue, gene, run_mode), dpi=150,
                bbox_inches='tight',
                pad_inches=0
                )
    plt.close()

    shap_sum = np.abs(train_shap_values).mean(axis=0)
    # print(shap_sum)

    d = {'column_name': features_df_new.columns.tolist(), 'shap_importance': shap_sum.tolist()}
    importance_df = pd.DataFrame(d)
    # print(importance_df)
    importance_df = importance_df.sort_values('shap_importance', ascending=False)
    importance_df.to_csv('%s/shap_importance_df_%s_%s_%s.csv' % (results_dir, tissue, gene, run_mode))

    #####################################
    # dependence plot

    shap.dependence_plot("rank(0)", train_shap_values, X_train_new,
                         feature_names=feature_names,
                         show=False)

    plt.savefig('%s/shap_dependence_top_feature_%s_%s_%s.png' % (results_dir, tissue, gene, run_mode), dpi=150,
                bbox_inches='tight',
                pad_inches=0)
    plt.close()

    shap.dependence_plot("rank(0)", train_shap_values, X_train_new, interaction_index="rank(1)",
                         feature_names=feature_names,
                         show=False)

    plt.savefig('%s/shap_dependence_top_and_second_feature_%s_%s_%s.png' %
                (results_dir, tissue, gene, run_mode), dpi=150,
                bbox_inches='tight',
                pad_inches=0)
    plt.close()

    shap.dependence_plot("rank(0)", train_shap_values, X_train_new, interaction_index=None,
                         feature_names=feature_names,
                         show=False)

    plt.savefig('%s/shap_dependence_top_feature_no_interaction_%s_%s_%s.png' %
                (results_dir, tissue, gene, run_mode), dpi=150,
                bbox_inches='tight',
                pad_inches=0)
    plt.close()

    ######################################

    if X_test_new is not None:

        shap.summary_plot(test_shap_values, X_test_new, feature_names=feature_names,
                          max_display=20,
                          show=False,
                          # plot_size=[12, 6]
                          )
        plt.savefig('%s/shap_test_summary_%s_%s_%s.png' % (results_dir, tissue, gene, run_mode),
                    dpi=150, bbox_inches='tight',
                    pad_inches=0)
        plt.close()

        shap.decision_plot(expected_value, test_shap_values, X_test_new,  # link="logit",
                           feature_names=feature_names,
                           show=False
                           )
        plt.savefig('%s/shap_test_decision_plot_%s_%s_%s.png' % (results_dir, tissue, gene, run_mode),
                    dpi=150, bbox_inches='tight',
                    pad_inches=0)
        plt.close()


def run_infer_mechanisms(job_name='', run_mode='', tissues=[], disease_genes={}):

    if run_mode is 'infer_genes':
        tissues = TISSUES_LIST_GENE_MODE
        disease_genes = DISEASE_GENES_DICT_GENE_MODE

    elif run_mode is 'infer_tissues':
        tissues = TISSUES_LIST_TISSUE_MODE
        disease_genes = DISEASE_GENES_DICT_TISSUE_MODE

    results_dir = '../results/%s' % job_name

    # TODO: WARNING - hard-coded to pre-computed dataset
    dataset = pd.read_csv("../article_results/dataset/df_complete_dataset_ready_adapted_no_missing_values.csv",
                          index_col=0)

    num_lables = sum('_causal' in s for s in list(dataset.columns.values))

    num_features = len(list(dataset.columns.values)) - num_lables

    # split data into X and y
    X = dataset.iloc[:, :num_features]
    Y = dataset.iloc[:, num_features:]
    features_df = X.copy()

    clf = XGBClassifier(
        objective='binary:logistic',
        learning_rate=0.1,
        gamma=0,
        max_depth=9,
        min_child_weight=1,
        colsample_bytree=1,
        subsample=1,
        # n_estimators=best_num_trees,
        n_estimators=150,
        scale_pos_weight=0.01,
        random_state=123,
        reg_lambda=0
    )

    for tissue in tissues:
        print(tissue)

        if Y[tissue].sum() < 50:
            continue

        if run_mode is 'infer_genes':

            for gene in disease_genes[tissue]:
                X_test = features_df.loc[gene, :]
                X_train = features_df[~features_df.index.isin([gene])]
                y_train = Y[~Y.index.isin([gene])].loc[:, tissue]

                X_train_new, X_test_new, cols = select_features_shap(X_train=X_train, Y=y_train,
                                                                     X_test=X_test, SEED=123)

                feature_names = list(features_df.columns.values[cols])
                features_df_new = features_df.iloc[:, cols]

                infer_models(clf=clf,
                             X_train_new=X_train_new,
                             X_test_new=X_test_new,
                             y_train=y_train,
                             feature_names=feature_names,
                             features_df_new=features_df_new,
                             tissue=tissue,
                             gene=gene,
                             run_mode=run_mode,
                             results_dir=results_dir)

        elif run_mode is 'infer_tissues':

            X_train = X
            y_train = Y[tissue]
            X_test = None

            X_train_new, X_test_new, cols = select_features_shap(X_train=X_train, Y=y_train,
                                                                 X_test=X_test, SEED=123)

            feature_names = list(features_df.columns.values[cols])
            features_df_new = features_df.iloc[:, cols]

            infer_models(clf=clf,
                         X_train_new=X_train_new,
                         X_test_new=X_test_new,
                         y_train=y_train,
                         feature_names=feature_names,
                         features_df_new=features_df_new,
                         tissue=tissue,
                         gene='',
                         run_mode=run_mode,
                         results_dir=results_dir)

    return


