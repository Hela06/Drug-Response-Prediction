# DRUG RESPONSE PREDICTION

## INTRODUCTION
Predicting the clinical response of therapeutic agents is a major challenge in cancer treatment. 
Unfortunately, establishing accurate and robust methods in such a context is difficult. 
Most approaches exploit only clinical data and gene expression information without considering other resources, such as pathways, capable of shedding light on the comprehension of such phenomena. 
We present a new method for drug response prediction which exploits pathway alteration information retrieved by PHENSIM.

![alt text](https://github.com/Hela06/Drug-Response-Prediction/blob/main/docs/images/Workflow-drug-prediction-1-6step.png)

We propose a workflow consisting of 6 main steps. 
1.  TCGA-curated drug response clinical data [1][1] are retrieved and preprocessed as described in [2][2].  
2.  Each cancer-drug dataset is divided into train/test sets. 
3.  The differentially expressed genes are identified in each sample by computing the 25th and 75th percentile of the gene distribution in normal samples.
4.  For each gene, contingency matrix is built including the number of responders and non-responders samples having that gene greater than 75 percentile or less than 25 percentile.
5. A chi-square test is applied to select the enrichment genes list. Such a list is used as input to the PHENSIM [3][3] simulator, together with expression data for each patient, to get the perturbation of genes. 
6. The PHENSIM simulator results are used to train a learning model capable of predicting drug response.

## USAGE
TODO

# References
[1] Ding Z, Zu S, Gu J. Evaluating the molecule-based prediction of clinical drug responses in cancer. Bioinformatics. 2016;32: 2891–2895.
[2] Hermida LC, Gertz EM, Ruppin E. Predicting cancer prognosis and drug response from the tumor microbiome. Nat Commun. 2022;13: 2896.
[3] Alaimo S, Rapicavoli RV, Marceca GP, La Ferlita A, Serebrennikova OB, Tsichlis PN, et al. PHENSIM: Phenotype Simulator. PLoS Comput Biol. 2021;17: e1009069.


[1]: https://pubmed.ncbi.nlm.nih.gov/27354694/ "Ding Z, Zu S, Gu J. Evaluating the molecule-based prediction of clinical drug responses in cancer. Bioinformatics. 2016;32: 2891–2895."
[2]: https://www.nature.com/articles/s41467-022-30512-3 "Hermida LC, Gertz EM, Ruppin E. Predicting cancer prognosis and drug response from the tumor microbiome. Nat Commun. 2022;13: 2896."
[3]: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009069  "Alaimo S, Rapicavoli RV, Marceca GP, La Ferlita A, Serebrennikova OB, Tsichlis PN, et al. PHENSIM: Phenotype Simulator. PLoS Comput Biol. 2021;17: e1009069."
