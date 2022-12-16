# DRUG RESPONSE PREDICTION
Predicting the clinical response of therapeutic agents is a major challenge in cancer treatment. 
Unfortunately, establishing accurate and robust methods in such a context is difficult. 
Most approaches exploit only clinical data and gene expression information without considering other resources, such as pathways, capable of shedding light on the comprehension of such phenomena. 
We present a new method for drug response prediction which exploits pathway alteration information retrieved by PHENSIM.

![alt text](https://github.com/Hela06/Drug-Response-Prediction/blob/main/docs/images/Workflow-drug-prediction-1-6step.png)

We propose a workflow consisting of 6 main steps. 
First, TCGA-curated drug response clinical data [1] are retrieved and preprocessed as described in [2].  
Secondly, each cancer-drug dataset is divided into train/test sets. Then, we identify differentially expressed genes in each sample by computing the 25th and 75th percentile of the gene distribution in normal samples.
Next, for each gene, we build a contingency matrix including the number of responders and non-responders samples having that gene greater than 75 percentile or less than 25 percentile. Finally, a chi-square test is applied to select the enrichment genes list. Such a list is used as input to the PHENSIM [3] simulator, together with expression data for each patient, to get the perturbation of genes. 
Finally, the PHENSIM simulator results are used to train a learning model capable of predicting drug response.
