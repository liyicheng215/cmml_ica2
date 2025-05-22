# CMML-ICA2: Comparative Evaluation of scPred, SingleR, and scANVI for Single-Cell Annotation

## Overview

This repository provides a systematic benchmarking of three automated cell annotation methods—**scPred**, **SingleR**, and **scANVI**—evaluated on both cross-dataset and large-scale in-dataset scenarios using standardized PBMC datasets.

The study assesses each method's performance in terms of accuracy, F1 score, computational efficiency, and scalability, offering insights into their applicability across varying dataset complexities.

## Repository Structure

* **`data_process/`**: Scripts for data preprocessing and change.
* **`run_model/`**: Implementation of training and prediction pipelines for each annotation method.
* **`result_csv/`**: Generated prediction labels and corresponding ground truth for evaluation.
* **`visualization/`**: Scripts for generating performance metrics visualizations, such as bar plots and confusion matrices.
* **`requirement.txt`**: List of Python and R dependencies required to run the analysis.
* **`README.md`**: Project documentation and usage guidelines.

## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/liyicheng215/cmml_ica2.git
   cd cmml_ica2
   ```

2. **Set up the Python environment according to `requirement.txt`**


## Usage

### Data Preprocessing

Navigate to the `data_process/` directory and execute the preprocessing scripts to prepare the datasets for analysis. 

### Model Training and Prediction

In the `run_model/` directory, you'll find separate scripts for each annotation method:

* **scPred**: Utilizes PCA for dimensionality reduction followed by SVM for classification.
* **SingleR**: Employs correlation-based matching against reference datasets.
* **scANVI**: Leverages a semi-supervised variational autoencoder framework for label transfer.

Script is configured to handle cross-dataset and in-dataset evaluation scenarios.

### Evaluation and Visualization

Post-prediction, use the scripts in the `visualization/` directory to compute performance metrics and generate visual representations:

* **Bar Plots**: Depicting overall accuracy and F1 scores across methods.
* **Per-Cell-Type F1 Scores**: Highlighting method performance on individual cell types.
* **Confusion Matrices**: Illustrating true vs. predicted labels for detailed error analysis.

## Results 
Generated prediction labels and corresponding ground truth for evaluation. These tables are available in the `result_csv/` directory.
