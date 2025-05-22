# cmml_ica2
Based on the contents of the GitHub repository [liyicheng215/cmml\_ica2](https://github.com/liyicheng215/cmml_ica2) and our previous discussions, here's a comprehensive `README.md` file tailored for your project:

---

# CMML-ICA2: Benchmarking Automated Cell Annotation Methods in Single-Cell RNA-seq Analysis

## Overview

This repository provides a systematic benchmarking of three automated cell annotation methods—**scPred**, **SingleR**, and **scANVI**—evaluated on both cross-dataset and large-scale in-dataset scenarios using standardized PBMC datasets.

The study assesses each method's performance in terms of accuracy, F1 score, computational efficiency, and scalability, offering insights into their applicability across varying dataset complexities.

## Repository Structure

* **`data_process/`**: Scripts for data preprocessing, including normalization and feature selection.
* **`run_model/`**: Implementation of training and prediction pipelines for each annotation method.
* **`result_csv/`**: Generated prediction labels and corresponding ground truth for evaluation.
* **`visualization/`**: Scripts for generating performance metrics visualizations, such as bar plots and confusion matrices.
* **`requirement.txt`**: List of Python dependencies required to run the analysis.
* **`README.md`**: Project documentation and usage guidelines.

## Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/liyicheng215/cmml_ica2.git
   cd cmml_ica2
   ```



2. **Set up the Python environment:**

   ```bash
   pip install -r requirement.txt
   ```



*Note: Ensure you have Python 3.8 or higher installed.*

## Usage

### Data Preprocessing

Navigate to the `data_process/` directory and execute the preprocessing scripts to prepare the datasets for analysis. This includes normalization, selection of highly variable genes (HVGs), and formatting data structures compatible with each annotation method.

### Model Training and Prediction

In the `run_model/` directory, you'll find separate scripts for each annotation method:

* **scPred**: Utilizes PCA for dimensionality reduction followed by SVM for classification.
* **SingleR**: Employs correlation-based matching against reference datasets.
* **scANVI**: Leverages a semi-supervised variational autoencoder framework for label transfer.

Each script is configured to handle both cross-dataset and in-dataset evaluation scenarios.

### Evaluation and Visualization

Post-prediction, use the scripts in the `visualization/` directory to compute performance metrics and generate visual representations:

* **Bar Plots**: Depicting overall accuracy and F1 scores across methods.
* **Per-Cell-Type F1 Scores**: Highlighting method performance on individual cell types.
* **Confusion Matrices**: Illustrating true vs. predicted labels for detailed error analysis.

## Results 
Generated prediction labels and corresponding ground truth for evaluation. These tables are available in the `result_csv/` directory.
