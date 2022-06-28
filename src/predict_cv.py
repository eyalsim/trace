import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy import interp
from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import *
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from sklearn.neural_network import MLPClassifier
from sklearn.svm import LinearSVC
from sklearn.utils.class_weight import compute_sample_weight
from xgboost.sklearn import XGBClassifier

#########################################################


SPECIFIED_TISSUES = ['whole_brain_causal',
                     'heart_causal', 'muscle_skeletal_causal', 'skin_causal',
                     'liver_causal', 'nerve_tibial_causal', 'testis_causal', 'whole_blood_causal',
                     'cortex_causal', 'cerebellum_causal', 'brain_other_causal']


def select_features(X_train, Y, X_test, SEED):

    print("Feature selection")

    lsvc = LinearSVC(C=0.1, penalty="l1", dual=False, class_weight={0: 0.01, 1: 1}, tol=0.001,
                     max_iter=10000, random_state=SEED).fit(X_train, Y)

    model = SelectFromModel(lsvc, prefit=True)

    X_train_new = model.transform(X_train)
    X_test_new = model.transform(X_test)

    print('X_train_new.shape:', X_train_new.shape)
    print('X_test_new.shape:', X_test_new.shape)

    return X_train_new, X_test_new


def load_adapted_dataset(dataset_dir):
    if dataset_dir == '':
        df = pd.read_csv("../article_results/dataset/df_complete_dataset_ready_adapted_no_missing_values.csv",
                         index_col=0)
    else:
        df = pd.read_csv("%s/df_complete_dataset_ready_adapted_no_missing_values.csv" % dataset_dir,
                         index_col=0)

    num_lables = sum('_causal' in s for s in list(df.columns.values))
    print('num_lables :', num_lables)

    num_features = len(list(df.columns.values)) - num_lables
    print('num_features: ', num_features)

    return df, num_features


def run_trace_cv(job_name='', specified_tissues=SPECIFIED_TISSUES):

    results_dir = '../results/%s' % job_name
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    if not os.path.exists(results_dir + "/folds_gene_lists"):
        os.makedirs(results_dir + "/folds_gene_lists")

    # cross-validation - part 1

    # dataset, num_features = load_adapted_dataset(results_dir)
    dataset, num_features = load_adapted_dataset('')  # hard-coding for article pre-computed dataset

    # split data into X and y
    X = dataset.iloc[:, :num_features]
    Y = dataset.iloc[:, num_features:]

    X = X.values

    if not specified_tissues:
        tissues = Y.columns.values
    else:
        tissues = specified_tissues

    print('tissues for analysis:', tissues)

    clf1 = XGBClassifier(
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
    clf2 = RandomForestClassifier(class_weight={0: 0.01, 1: 1}, n_estimators=1000, random_state=None, max_depth=None)
    clf3 = LogisticRegression(class_weight={0: 0.01, 1: 1}, random_state=123, solver='lbfgs', C=1, max_iter=100000)
    clf4 = GradientBoostingClassifier(n_estimators=80, random_state=123,
                                      init=clf3, loss='deviance')
    clf5 = MLPClassifier(hidden_layer_sizes=(10, 10,), max_iter=1000, alpha=0.5, activation='relu',
                         tol=0.001, batch_size=200, solver='adam', verbose=10)

    meta_mlp8 = MLPClassifier(hidden_layer_sizes=(10, 10,), max_iter=1000, alpha=0.1, activation='relu',
                              solver='adam', verbose=10, learning_rate_init=0.01, learning_rate='invscaling')

    classifiers = [clf1, clf2, clf3, clf4, clf5]
    methods = ['XGB', 'RF', 'LR', 'LR+GB', 'MLP']

    meta_classifiers = [meta_mlp8]
    meta_methods = ['meta_MLP']

    classifiers = classifiers + meta_classifiers
    methods = methods + meta_methods

    results_df = pd.DataFrame(columns=tissues, index=methods)

    results_pr_df = pd.DataFrame(columns=tissues, index=methods)

    #########################################################

    # cross-validation - part 2

    # Run classifier with cross-validation and plot ROC curves
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=123)

    for tissue in tissues:

        print(tissue)

        print(Y[tissue].sum())
        if Y[tissue].sum() < 50:
            continue

        Y_pred_prob_df = pd.DataFrame(columns=methods, index=Y.index.values)

        mean_tpr_clfs = []
        mean_fpr_clfs = []
        mean_auc_clfs = []
        std_auc_clfs = []

        lsvc = LinearSVC(C=0.1, penalty="l1", dual=False, class_weight={0: 0.01, 1: 1}, tol=0.001,
                         max_iter=10000, random_state=0)

        features_names = dataset.iloc[:, :num_features].columns
        features_names = list(features_names)
        lsvc.fit(X, Y[tissue])
        coefficients = lsvc.coef_.tolist()[0]

        feature_coefficients = pd.DataFrame(
            {'features': features_names,
             'coefficients': coefficients})
        feature_coefficients.to_csv('%s/feature_coefficients_%s.csv' % (results_dir, tissue))

        for method_name, clf in zip(methods, classifiers):

            tprs = []
            aucs = []
            mean_fpr = np.linspace(0, 1, 100)

            mean_average_pr_score = []

            i = 0
            for train, test in cv.split(X, Y[tissue]):

                s = Y[tissue]
                s = s[s != 0]

                train_gene_list = Y_pred_prob_df.iloc[train].index.values.tolist()
                with open('%s/folds_gene_lists/%s_train_complete_gene_list_fold_%s.txt' % (results_dir, tissue, i),
                          'w') as f:
                    for gene in train_gene_list:
                        f.write("%s\n" % gene)

                positive_all_gene_list = s.index.values.tolist()
                positive_train_gene_list = [value for value in train_gene_list if value in positive_all_gene_list]
                with open('%s/folds_gene_lists/%s_train_positive_gene_list_fold_%s.txt' % (results_dir, tissue, i),
                          'w') as f:
                    for gene in positive_train_gene_list:
                        f.write("%s\n" % gene)

                test_gene_list = Y_pred_prob_df.iloc[test].index.values.tolist()
                with open('%s/folds_gene_lists/%s_test_gene_list_fold_%s.txt' % (results_dir, tissue, i), 'w') as f:
                    for gene in test_gene_list:
                        f.write("%s\n" % gene)

                if method_name in meta_methods:
                    print('meta learner: ', method_name)

                    meta_X = Y_pred_prob_df.iloc[:, :(len(methods) - len(meta_methods))]
                    meta_X = meta_X.values

                    X_train = meta_X[train]
                    X_test = meta_X[test]

                else:
                    X_train = X[train]
                    X_test = X[test]

                y_train = Y[tissue][train]
                y_test = Y[tissue][test]

                if method_name == 'XGB' or method_name == 'LR+GB':
                    X_train_new, X_test_new = select_features(X_train=X_train, Y=y_train, X_test=X_test, SEED=123)
                elif method_name == 'LR':
                    X_train_new, X_test_new = select_features(X_train=X_train, Y=y_train, X_test=X_test, SEED=123)
                else:
                    X_train_new, X_test_new = X_train, X_test

                if method_name == 'GB':
                    clf.fit(X=X_train_new, y=y_train, sample_weight=compute_sample_weight({0: 0.01, 1: 1}, y_train))
                else:
                    clf.fit(X_train_new, y_train)

                print('fitting')
                clf.fit(X_train_new, y_train)

                print('predicting')
                probas_ = clf.predict_proba(X_test_new)
                # print('clf.classes_:, ', clf.classes_)

                print('fold', i)
                print(clf)
                print('AUC: ', roc_auc_score(y_test, probas_[:, 1]))

                y_pro = probas_[:, 1].copy().reshape(-1, 1)
                scaler = preprocessing.MinMaxScaler(feature_range=(0, 10))

                y_pro_scaled = scaler.fit_transform(y_pro)

                print('AUC_scaled_scores: ', roc_auc_score(y_test, y_pro_scaled))

                # insert predicted probabilities to results dataframe
                k = 0
                for j in test:
                    Y_pred_prob_df[method_name].iloc[j] = y_pro_scaled.copy()[k][0]
                    k += 1

                ####################################################
                # calculate precision-recall curve
                precision, recall, thresholds = precision_recall_curve(y_test, y_pro_scaled)
                # calculate precision-recall AUC
                precision_recall_auc = auc(recall, precision)
                print('precision_recall_auc: ', precision_recall_auc)
                print('average_precision_score: ', average_precision_score(y_test, y_pro_scaled))
                mean_average_pr_score.append(average_precision_score(y_test, y_pro_scaled))
                ####################################################

                # Compute ROC curve and area the curve
                fpr, tpr, thresholds = roc_curve(y_test, y_pro_scaled)
                tprs.append(interp(mean_fpr, fpr, tpr))
                tprs[-1][0] = 0.0
                roc_auc = auc(fpr, tpr)
                aucs.append(roc_auc)

                plt.plot(fpr, tpr, lw=1, alpha=0.3,
                         label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))
                i += 1

            plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
                     label='Chance', alpha=.8)

            mean_tpr = np.mean(tprs, axis=0)
            mean_tpr[-1] = 1.0
            mean_auc = auc(mean_fpr, mean_tpr)
            std_auc = np.std(aucs)
            plt.plot(mean_fpr, mean_tpr, color='b',
                     label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
                     lw=2, alpha=.8)

            std_tpr = np.std(tprs, axis=0)
            tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
            tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
            plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                             label=r'$\pm$ 1 std. dev.')

            plt.xlim([-0.05, 1.05])
            plt.ylim([-0.05, 1.05])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title('Receiver operating characteristic plot: %s_%s' % (tissue, method_name))
            plt.legend(loc="lower right")
            plt.savefig("%s/ROC_%s_%s.png" % (results_dir, tissue, method_name))
            plt.close()

            results_df.loc[method_name, tissue] = mean_auc

            mean_pr_score = np.mean(mean_average_pr_score)
            print('mean_pr_score: ', mean_pr_score)
            results_pr_df.loc[method_name, tissue] = mean_pr_score

            mean_fpr_clfs.append(mean_fpr)
            mean_tpr_clfs.append(mean_tpr)
            mean_auc_clfs.append(mean_auc)
            std_auc_clfs.append(std_auc)

        for i in range(len(methods)):
            plt.plot(mean_fpr_clfs[i], mean_tpr_clfs[i], lw=2, alpha=0.5,
                     label=r'%s (AUC = %0.2f $\pm$ %0.2f)' % (methods[i], mean_auc_clfs[i], std_auc_clfs[i]))
        plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='black',
                 label='Chance', alpha=.8)
        plt.xlim([-0.05, 1.05])
        plt.ylim([-0.05, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver operating characteristic plot - %s' % tissue)
        plt.legend(loc="lower right")
        plt.savefig("%s/ROC_%s_all_methods.png" % (results_dir, tissue))
        plt.close()

        Y_pred_prob_df = pd.concat([Y_pred_prob_df, Y[tissue]], axis=1, sort=False)
        Y_pred_prob_df.rename(columns={tissue: 'Class'}, inplace=True)
        Y_pred_prob_df.to_csv('%s/Y_pred_prob_%s.csv' % (results_dir, tissue))

    results_df.to_csv('%s/methods_tissues_auc.csv' % results_dir)
    results_pr_df.to_csv('%s/methods_tissues_pr_score.csv' % results_dir)
