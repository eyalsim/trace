## TRACE: Predicting molecular mechanisms of hereditary diseases by using their tissue-selective manifestation

This repository contains data, code, and analysis as described in the paper "Predicting molecular mechanisms of hereditary diseases by using their tissue-selective manifestation" (*https://doi.org/10.15252/msb.202211407*).


**Webserver:** A catalog of the tissue-specific risks associated with each protein-coding gene in our dataset is available through the TRACE webserver (https://netbio.bgu.ac.il/trace/).

## Table of Contents
- [Article results](#article-results)
  - [Pipeline code](#pipeline-code)
  - [Plots and analyses](#plots-and-analyses)
- [Requirements](#requirements)
- [Data availability](#data-availability)

### Article results
Raw results are found in the respective "article_results/figureX/" folders.
#### Pipeline code
Trace main pipeline was written in python and can be reviewed and used by running "run_trace.py" from the src folder.
#### Plots and analyses
- Plots and analyses were mostly done in R.\
These are found in the "article_results" folder and can be re-created from article_results/src as long as the [Trace dataset](#data-availability) (both files) is downloaded to trace/article_results/dataset and the trace [data](#data-availability) folder is downloaded to trace/.
- Analyses that were performed as part of the main pipeline can be recreated by running the relevant module via "run_trace.py".

| Figure  | Panel | Script                                         |
|---------|-------|------------------------------------------------|
| Figure1 | B     | trace_dataset_statistics.R                     |
| Figure2 | A,B   | example_genes_model_interpratation.R           |
|         |       | example_process.R                              |
|         |       | infer_mechanisms.py (trace pipeline)           |
| Figure3 | A     | infer_mechanisms.py (trace pipeline)           |
|         | B,C   | features_ballon_plot.R                         |
|         |       | shap_feature_type_analysis.py (trace pipeline) |
| Figure4 | A     | all_tissues_methods_auc_barplot.R              |
|         | C     | trace_auc_barplot.R                            |
|         | D     | trace_boxplots.R                               |
|         | E     | trace_ridge_density.R                          |
| Figure5 | C     | trace_patients_percentile.R                    |
|         | D     | compare_trace_with_other_methods_ranks.R       |
| S3      | A     | xgb_auc_barplot.R                              |
|         | B     | predict_cv.py (trace pipeline)                 |
| S4      |       | infer_mechanisms.py (trace pipeline)           |
| S5      |       | predict_cv.py (trace pipeline)                 |
| S6      | A     | all_tissues_methods_auc_barplot.R              |
|         | B     | trace_auc_barplot.R                            |
| S7      | A     | trace_boxplots.R                               |
|         | B     | trace_ridge_density.R                          |
| S8      |       | widely_expressed_analysis.R                    |
| S9      | A     | trace_dataset_statistics.R                     |
|         | B     | trace_boxplots.R                               |
|         | C     | trace_ridge_density.R                          |

### Requirements
TRACE can be trained using a standard laptop (runtime ranges between ~few minutes/hours for specific prediction or analysis and ~days for the entire pipeline including recreating the features dataset.  
Required dependencies are detailed in **environment.yml**.

### Data availability
- A **data** folder and two **features datasets** are required to train all TRACE models and can be downloaded from https://sandbox.zenodo.org/record/1185590#.ZES8M3ZByUk.
