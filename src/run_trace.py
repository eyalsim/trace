from create_trace_dataset import run_create_trace_dataset
from rename_features import rename_dataset_features
from prepare_trace_dataset import run_prepare_dataset
from infer_mechanisms import run_infer_mechanisms
from shap_feature_type_analysis import run_feature_type_analysis
from predict_cv import run_trace_cv
import os
import datetime

JOB_NAME = ''  # enter job name here (otherwise date and time)

CREATE_DATASET = False
PREPARE_DATASET = False
INFER_SELECTIVE_MECHANISMS = True
TRACE_CV = True

THRESHOLD = 7
FC_THRESHOLD = 2
NET_EMBEDDING = True


# TODO: (nice-to-haves for future)
#  - move all module constants here
#  - argparse
#  - More documentation where needed
#  - improve flow between modules


def run_trace(job_name=''):

    if job_name == '':
        job_name = str(datetime.datetime.now())

    results_dir = '../results/%s' % job_name
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    if CREATE_DATASET:
        run_create_trace_dataset(job_name, THRESHOLD, FC_THRESHOLD, NET_EMBEDDING)
        rename_dataset_features(job_name)

    if PREPARE_DATASET:
        run_prepare_dataset(job_name=job_name, features_to_keep=None)

    if INFER_SELECTIVE_MECHANISMS:
        # currently, hard-coded for Fig.2 example disease genes.
        run_infer_mechanisms(job_name=job_name, run_mode='infer_genes')
        run_infer_mechanisms(job_name=job_name, run_mode='infer_tissues')
        run_feature_type_analysis(job_name=job_name)

    if TRACE_CV:
        # currently, hard-coded to using article pre-computed dataset.
        run_trace_cv(job_name=job_name)

    return


run_trace(job_name=JOB_NAME)
