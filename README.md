## Contents
- [Overview](#overview)
- [Files](#files)
- [Usage](#usage)
- [Publication](#publication)


## Overview
Protein activity is a significant characteristic for recombinant proteins which can be used as biocatalysts. Since protein activity and solubility are correlated for some proteins, the publicly available solubility dataset may be adopted to develop models that can predict protein solubility from sequence. In literature, predicting protein solubility from sequence has been intensively explored, but the predicted solubility represented in binary values from all the developed models was not suitable for guiding experimental designs to improve protein solubility. We first implemented a novel approach that predicted protein solubility in continuous numerical values instead of binary ones. After combining it with various machine learning algorithms, we achieved a R2 of 0.4115 when Support Vector Machine (SVM) algorithm was used. 


## Files
* seqinfo.xlsx: The protein sequence we downloaded from NCBI database
* all_data.xlsx: The original eSol database including protein name, protein solubility and other information
* solubility2.csv: data after matching the protein sequence and solubility by name
* cleandata_1name.csv: The amino acid composition exacting from sequence and solubility
* solubility prediction-final.R: The model training and evaluation of machine learning models except ANN(ANNs-solubility-final.ipynb) and XGboost(XGboost-final.ipynb)


## Usage
* Run prediction-final.R to clean, analyse data and predict protein solubility with logistic regression, decision tree, SVM, Naïve Bayes, random forest
* Run ANNs-solubility-final.ipynb to predict protein solubility with ANN
* Run XGboost-final.ipynb to predict protein solubility with XGBoost


## Publication
Han, X., Wang, X., & Zhou, K. (2019). Develop machine learning-based regression predictive models for engineering protein solubility. Bioinformatics, 35(22), 4640-4646.




