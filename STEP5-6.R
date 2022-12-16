#STEP5-6:TRAIN THE MODEL AND TEST RESPONSE
library(data.table)
library(dplyr)
library(e1071)
library(ROCR)
s_vec_mach <- function(EP,EA,PP,PA,C,sample_info, prj){
  EP <- readRDS(EP)
  EA <- readRDS(EA)
  PP <- readRDS(PP)
  PA <- readRDS(PA)
  C <- readRDS(C)
  info <- readRDS(sample_info)
  
  for (d in unique(info$drug_name[which(info$PROJECT == prj)])) {
    matr <- c("END-PE","END-AC","PATH-PE","PATH-AC","CLIN", "EXPR-SB", "END-CLIN")
    df <- data.frame(MODELS=matr, MEAN_ACC=rep(NA,length(matr)), TPR=rep(NA,length(matr)),FPR=rep(NA,length(matr)),
                     TP=rep(NA,length(matr)), FP=rep(NA,length(matr)), TN=rep(NA,length(matr)), FN=rep(NA,length(matr)))
    for (i in matr) {
      if(i == "END-PE"){
        print(paste(prj,i,d))
        m <- EP
        m[,1] <- NULL
        m <- m[complete.cases(m),]
        m_t <- t(m)
      }
      else if (i == "END-AC"){
        print(paste(prj,i,d))
        m <- EA
        m[,1] <- NULL
        m <- m[complete.cases(m),]
        m_t <- t(m)
      }
      else if (i == "PATH-PE"){
        print(paste(prj,i,d))
        m <- PP
        m[,c(1:3)] <- NULL
        m <- m[complete.cases(m),]
        m_t <- t(m)
      }
      else if (i == "PATH-AC"){
        print(paste(prj,i,d))
        m <- PA
        m[,c(1:3)] <- NULL
        m <- m[complete.cases(m),]
        m_t <- t(m)
      } 
      else if (i == "CLIN"){
        print(paste(prj,i,d))
        m_t <- C
        row.names(m_t) <- paste(m_t$PROJECT, m_t$drug_name, m_t$CONDITION, m_t$'TR/TE',m_t$SAMPLES, sep="_")
        m_t <- m_t[,c(6:8)]
        m_t$tumor_stage <- as.numeric(m_t$tumor_stage)
      } 
      #else if(i == "EXPR"){
      #print(paste(prj,i,d))
      #m <- readRDS(paste0("./expression_matrix/TCGA-",prj,"_",d,"_expr_matrix_tumor.rds"))
      #colnames(m) <- substring(colnames(m),1,12)
      #for (c in 1:ncol(m)) {
      #r <- which(info$SAMPLES == colnames(m)[c] & info$drug_name == d)
      #new <- paste(prj,d,info$CONDITION[r],info$`TR/TE`[r],info$SAMPLES[r],sep="_")
      #colnames(m)[c] <- new
      #}
      #m_t <- t(m)
      #}
      else if(i == "EXPR-SB"){
        print(paste(prj,i,d))
        m <- readRDS(paste0("./expression_matrix/TCGA-",prj,"_",d,"_expr_matrix_tumor.rds"))
        colnames(m) <- substring(colnames(m),1,12)
        row.names(m) <- gsub("\\..*","",row.names(m))
        str <- readRDS(paste0("./gene_stratified/TCGA-",prj,"_",d,"_genes_stratification_pval.rds"))
        str <- subset(str, chiTest_pv < 0.05)
        gene_sb <- row.names(str)
        m <- subset(m,row.names(m) %in% gene_sb)
        for (c in 1:ncol(m)) {
          r <- which(info$SAMPLES == colnames(m)[c] & info$drug_name == d)
          new <- paste(prj,d,info$CONDITION[r],info$`TR/TE`[r],info$SAMPLES[r],sep="_")
          colnames(m)[c] <- new
        }
        m_t <- t(m)
      }
      else if (i == "END-CLIN"){
        print(paste(prj,i,d))
        m <- EP
        m[,1] <- NULL
        m <- m[complete.cases(m),]
        m <- t(m)
        m_t <- C
        row.names(m_t) <- paste(m_t$PROJECT, m_t$drug_name, m_t$CONDITION, m_t$'TR/TE',m_t$SAMPLES, sep="_")
        m_t <- m_t[,c(6:8)]
        m_t$tumor_stage <- as.numeric(m_t$tumor_stage)
        m_t <- cbind(m, m_t)
      }
      
      print("split train test")
      m_t <- as.data.frame(m_t)
      m_t$CONDITION <- ifelse(row.names(m_t) %like% "NR", "NR","R")
      m_t$CONDITION = factor(m_t$CONDITION, levels = c("NR", "R"))
      
      m_train <- subset(m_t, row.names(m_t) %like% d & row.names(m_t) %like% "TRAIN" )
      #m_test <- subset(m_t, row.names(m_t)%like%d & row.names(m_t)%like%"TEST")
      m_test <- subset(m_t, row.names(m_t)%like%d)
      set.seed(246)
      print("svm analysis")
      #tune.out <- tune(svm, CONDITION~ ., data = m_train, kernel = "radial",
      #ranges = list(cost = c(0.1, 1, 10, 100, 1000),
      #gamma = c(0.5, 1, 2, 3, 4)))
      #tune.out$
      classifier = svm(formula = CONDITION ~ .,
                       data = m_train,
                       cost=10,
                       gamma=0.5,
                       type="C-classification",
                       kernel = 'radial',
                       cross = 3)
      df$MEAN_ACC[which(df$MODELS == i)] <- classifier$tot.accuracy
      
      print("Compute Prediction")
      m_test$P <- predict(classifier, newdata = m_test[-ncol(m_test)], type="response")
      print("Confusion matrix")
      #cm = table(m_test[,"CONDITION"], m_test$P)
      print("Select prediction")
      prdc <- m_test[,c(ncol(m_test)-1,ncol(m_test))]
      prdc <- ifelse(prdc == "NR", 0,1)
      prdc <- as.data.frame(prdc)
      
      print("Identify TP, FP, TN, FN")
      pr<-prediction(prdc$P, prdc$CONDITION)
      print("Identify TPR and FPR")
      pref <- performance(pr, "tpr", "fpr")
      print("Fill the dataframe")
      
      
      print("Insert TPR and FPR information")
      df$TPR[which(df$MODELS == i)] <- pref@y.values
      df$FPR[which(df$MODELS == i)] <- pref@x.values
      print("Insert TP, FP, TN, FN information")
      df$TP[which(df$MODELS == i)] <- pr@tp
      print("prova d")
      df$FP[which(df$MODELS == i)] <- pr@fp
      print("prova e")
      df$TN[which(df$MODELS == i)] <- pr@tn
      print("prova f")
      df$FN[which(df$MODELS == i)] <- pr@fn
      df <- tidyr::unnest_longer(df, c("TPR","FPR","TP","FP","TN","FN"))
    }
    saveRDS(df, paste0("./svm_perf/performance_TESTallDT_svm_",d,"_",prj,".rds"))
  }
}

s_vec_mach(EP="./matrix_Perturbation_BLCA.rds", 
           EA="./matrix_ActivationScore_BLCA.rds",
           PP="./matrix_each_pathway_Perturbation_BLCA.rds",
           PA="./matrix_each_pathway_ActivationScore_BLCA.rds",
           C="./clinical_data/BLCA_info_data.rds",
           sample_info = "./ALLsamples_info_TrainTest.rds",
           prj="BLCA")


