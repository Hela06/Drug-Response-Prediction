#STEP2-3:TRAIN/TEST SPLIT AND PHENSIM INPUT GENERATION 
library(TCGAbiolinks)
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(data.table)
library(plyr)
library(stringr)

drug_sim <- function(sample_info, m_expr, output, m_comb, strat){
  
  input <- readRDS(sample_info)
  input$PROJECT <- unlist(input$PROJECT)
  input <- input[!duplicated(input[c(1,2,3,4)]),]
  print("input file loaded!")
  
  input$CONDITION <- ifelse(input$CONDITION == "NON_RESPONDER", "NR","R")
  
  input$"TR/TE" <- NA
  saveRDS(input,"./ALLsamples_info_TrainTest.rds")
  
  project <- unique(input$PROJECT)
  
  for (p in project) {
    p_name <- paste0("TCGA-", p)
    drug <- unique(input$drug_name[which(input$PROJECT == p)])
    for (d in drug) {
      if(!file.exists(paste0(strat,p_name,"_",d,"_genes_stratification.rds"))){
        if(!file.exists(paste0(m_comb,p_name,"_",d,"_expr_matrix_comb_normalized.rds"))){
          print(paste("Load expression matrixs for project", p_name, d))
          ex_mat <- readRDS(paste0(m_expr, p_name, "expr_matrix_normal.rds"))
          p_mat <- readRDS(paste0(m_expr, p_name,"_",d, "_expr_matrix_tumor.rds"))
          
          comb <- TCGAanalyze_Normalization(tabDF = cbind(ex_mat,p_mat),
                                            geneInfoHT,
                                            method = "gcContent")
          comb <- na.omit(comb)
          saveRDS(comb, paste0(m_comb,p_name,"_",d,"_expr_matrix_comb_normalized.rds"))
          
          n <- ncol(p_mat)
          group_case <- rep("CASE",n)
          m <- ncol(ex_mat)
          group_ctrl <- rep("CONTROL", m)
          group_total <- c(group_case, group_ctrl) 
          d0 <- DGEList(comb)
          mm <- model.matrix(~0 + group_total)
          keep <- filterByExpr(d0, design = mm)
          d0 <- d0[keep,keep.lib.sizes=FALSE]
          d0 <- calcNormFactors(d0)
          y <- voom(d0, mm,plot = F)
          barcode <-colnames(p_mat)
          new_mat <- y$E
          new_mat <- new_mat[ , !(colnames(new_mat) %in% barcode)]
          new_p_mat <- y$E
          new_p_mat <- new_p_mat[,(colnames(new_p_mat) %in% barcode)]
          
          print("Compute 75th and 25th")
          df_normal <- matrix(nrow = nrow(new_mat), ncol = (3 + length(barcode)))
          row.names(df_normal) <- row.names(new_mat)
          colnames(df_normal) <- c("25th", "50th", "75th", barcode)
          
          for (n in row.names(df_normal)) {
            df_normal[n,c(1:3)] <- quantile(new_mat[n,], c(0.25,0.5,0.75))
          }
          
          print("Identify upregulated and downregulated genes")
          for (s in colnames(df_normal)[4:ncol(df_normal)]) {
            for (n in row.names(new_p_mat)) {
              if (new_p_mat[n,s] < df_normal[n,1]){
                df_normal[n,s] <- "DOWN_REG"
              }
              else if (new_p_mat[n,s] > df_normal[n,3]){
                df_normal[n,s] <- "UP_REG"
              }
              else { df_normal[n,s] <- "NO_CHANGE" }
            }
          }
          
          colnames(df_normal)[4:ncol(df_normal)]<- substring(colnames(df_normal)[4:ncol(df_normal)],1,12)
          for(c in 4:ncol(df_normal)){
            if(input$CONDITION[which(input$SAMPLES == colnames(df_normal)[c] & input$drug_name == d)] == "R"){
              colnames(df_normal)[c] <- paste0("RESPONDER-", colnames(df_normal)[c])
            }
            else if(input$CONDITION[which(input$SAMPLES == colnames(df_normal)[c]& input$drug_name == d)] == "NR") {
              colnames(df_normal)[c] <- paste0("PROGRESSION-", colnames(df_normal)[c])
            }
          }
          
          
          saveRDS(df_normal, paste0(strat, p_name,"_",d,"_genes_stratification.rds"))
          
          df_normal <- as.data.frame(df_normal)
          df_normal$chiTest_pv <- NA
          
          print("Generate contingency matrixs")
          for (g in row.names(df_normal)) {
            if(length(which(df_normal[g,] == "NO_CHANGE")) == (ncol(df_normal) - 4)){
              df_normal[g,"chiTest_pv"] <- "NO_CHANGE"
            } else {
              cont_m <- matrix(nrow = 2, ncol = 2)
              row.names(cont_m) <- c("R", "NR")
              colnames(cont_m) <- c("DOWN", "UP")
              cont_m[1,1] <- length(which(colnames(df_normal) %like% "RESPONDER" & df_normal[g,] == "DOWN_REG" ))
              cont_m[2,1] <- length(which(colnames(df_normal) %like% "PROGRESSION" & df_normal[g,] == "DOWN_REG" ))
              cont_m[1,2] <- length(which(colnames(df_normal) %like% "RESPONDER" & df_normal[g,] == "UP_REG" ))
              cont_m[2,2] <- length(which(colnames(df_normal) %like% "PROGRESSION" & df_normal[g,] == "UP_REG" ))
              ctst <- chisq.test(cont_m) 
              pv <- ctst$p.value
              df_normal[g,"chiTest_pv"] <- pv
            }
          }
          
          saveRDS(df_normal, paste0(strat, p_name,"_",d,"_genes_stratification_pval.rds"))
          
        }
        else {
          print(paste("Expression matrix combined and normalized exists for project", p_name, d))
          comb <- readRDS(paste0(m_comb, p_name, "_", d,"_expr_matrix_comb_normalized.rds"))
          ex_mat <- readRDS(paste0(m_expr, p_name, "expr_matrix_normal.rds"))
          p_mat <- readRDS(paste0(m_expr, p_name,"_",d, "_expr_matrix_tumor.rds"))
          n <- ncol(p_mat)
          group_case <- rep("CASE",n)
          m <- ncol(ex_mat)
          group_ctrl <- rep("CONTROL", m)
          group_total <- c(group_case, group_ctrl) 
          d0 <- DGEList(comb)
          mm <- model.matrix(~0 + group_total)
          keep <- filterByExpr(d0, design = mm)
          d0 <- d0[keep,keep.lib.sizes=FALSE]
          d0 <- calcNormFactors(d0)
          y <- voom(d0, mm,plot = F)
          barcode <-colnames(p_mat)
          new_mat <- y$E
          new_mat <- new_mat[ , !(colnames(new_mat) %in% barcode)]
          new_p_mat <- y$E
          new_p_mat <- new_p_mat[,(colnames(new_p_mat) %in% barcode)]
          
          print("Compute 75th and 25th")
          df_normal <- matrix(nrow = nrow(new_mat), ncol = (3 + length(barcode)))
          row.names(df_normal) <- row.names(new_mat)
          colnames(df_normal) <- c("25th", "50th", "75th", colnames(new_p_mat))
          
          
          for (n in row.names(df_normal)) {
            df_normal[n,c(1:3)] <- quantile(new_mat[n,], c(0.25,0.5,0.75))
          }
          
          print("Identify upregulated and downregulated genes")
          for (s in colnames(df_normal)[4:ncol(df_normal)]) {
            for (n in row.names(new_p_mat)) {
              if (new_p_mat[n,s] < df_normal[n,1]){
                df_normal[n,s] <- "DOWN_REG"
              }
              else if (new_p_mat[n,s] > df_normal[n,3]){
                df_normal[n,s] <- "UP_REG"
              }
              else { df_normal[n,s] <- "NO_CHANGE" }
              saveRDS(df_normal, paste0(strat, p_name,"_genes_stratification.rds"))
            }
          }
          
          colnames(df_normal)[4:ncol(df_normal)]<- substring(colnames(df_normal)[4:ncol(df_normal)],1,12)
          for(c in 4:ncol(df_normal)){
            if(input$CONDITION[which(input$SAMPLES == colnames(df_normal)[c] & input$drug_name == d)] == "R"){
              colnames(df_normal)[c] <- paste0("RESPONDER-", colnames(df_normal)[c])
            }
            else if(input$CONDITION[which(input$SAMPLES == colnames(df_normal)[c] & input$drug_name == d)] == "NR") {
              colnames(df_normal)[c] <- paste0("PROGRESSION-", colnames(df_normal)[c])
            }
          }
          
          saveRDS(df_normal, paste0(strat, p_name,"_",d,"_genes_stratification.rds"))
          df_normal <- as.data.frame(df_normal)
          df_normal$chiTest_pv <- NA
          
          print("Generate contingency matrixs")
          for (g in row.names(df_normal)) {
            if(length(which(df_normal[g,] == "NO_CHANGE")) == (ncol(df_normal) - 4)){
              df_normal[g,"chiTest_pv"] <- "NO_CHANGE"
            } else {
              cont_m <- matrix(nrow = 2, ncol = 2)
              row.names(cont_m) <- c("R", "NR")
              colnames(cont_m) <- c("DOWN", "UP")
              cont_m[1,1] <- length(which(colnames(df_normal) %like% "RESPONDER" & df_normal[g,] == "DOWN_REG" ))
              cont_m[2,1] <- length(which(colnames(df_normal) %like% "PROGRESSION" & df_normal[g,] == "DOWN_REG" ))
              cont_m[1,2] <- length(which(colnames(df_normal) %like% "RESPONDER" & df_normal[g,] == "UP_REG" ))
              cont_m[2,2] <- length(which(colnames(df_normal) %like% "PROGRESSION" & df_normal[g,] == "UP_REG" ))
              ctst <- chisq.test(cont_m) 
              pv <- ctst$p.value
              df_normal[g,"chiTest_pv"] <- pv
            }
          }
          
          saveRDS(df_normal, paste0(strat, p_name,"_",d,"_genes_stratification_pval.rds"))
          
        } 
      }
      else { 
        print(paste("Expression matrix combined and normalized exists for project", p_name, d))
        ex_mat <- readRDS(paste0(m_expr, p_name, "expr_matrix_normal.rds"))
        p_mat <- readRDS(paste0(m_expr, p_name,"_",d, "_expr_matrix_tumor.rds"))
        if(!file.exists(paste0(strat,p_name,"_",d,"_genes_stratification_pval.rds"))){
          print(paste("Matrix with genes stratification exists for project", p_name, d))
          df_normal <- readRDS(paste0(strat, p_name, "_", d,"_genes_stratification.rds"))
          df_normal <- as.data.frame(df_normal)
          df_normal$chiTest_pv <- NA
          
          print(paste(" Generate contingengy matrixs"))
          for (g in row.names(df_normal)) {
            if(length(which(df_normal[g,] == "NO_CHANGE")) == (ncol(df_normal) - 4)){
              df_normal[g,"chiTest_pv"] <- "NO_CHANGE"
            } else {
              cont_m <- matrix(nrow = 2, ncol = 2)
              row.names(cont_m) <- c("R", "NR")
              colnames(cont_m) <- c("DOWN", "UP")
              cont_m[1,1] <- length(which(colnames(df_normal) %like% "RESPONDER" & df_normal[g,] == "DOWN_REG" ))
              cont_m[2,1] <- length(which(colnames(df_normal) %like% "PROGRESSION" & df_normal[g,] == "DOWN_REG" ))
              cont_m[1,2] <- length(which(colnames(df_normal) %like% "RESPONDER" & df_normal[g,] == "UP_REG" ))
              cont_m[2,2] <- length(which(colnames(df_normal) %like% "PROGRESSION" & df_normal[g,] == "UP_REG" ))
              ctst <- chisq.test(cont_m) 
              pv <- ctst$p.value
              df_normal[g,"chiTest_pv"] <- pv
            }
          }
          
          saveRDS(df_normal, paste0(strat, p_name,"_",d,"_genes_stratification_pval.rds"))
        } else {
          print(paste("Stratification gene matrix with pvalue exists for project", p_name, d))
          df_normal <- readRDS(paste0(strat, p_name, "_", d,"_genes_stratification_pval.rds"))
        }
      }
      
      #Rimuovere tutti i geni che hanno un pvalue non significativo
      print("Remove no significant genes!")
      df_normal_sb <- subset(df_normal, chiTest_pv < 0.05)
      df_normal_sb[,c(1:3,ncol(df_normal_sb))] <- NULL
      
      #Train/test set split
      print("Divide the samples in training and test set!")
      info <- readRDS("./ALLsamples_info_TrainTest.rds")
      set.seed(246)
      sampling <- sample.int(n = ncol(df_normal_sb), size = floor(.75*ncol(df_normal)), replace = FALSE)
      for (z in 1:ncol(df_normal_sb)) {
        nm <- sub(".+?-", "",colnames(df_normal_sb)[z])
        q <- which(info$SAMPLES == nm & info$PROJECT == p & info$drug_name == d)
        if(z %in% sampling){
          info$"TR/TE"[q] <- "TRAIN"
        } else {
          info$"TR/TE"[q] <- "TEST"
        }
      }
      
      saveRDS(info, "./ALLsamples_info_TrainTest.rds")
      colnames(df_normal_sb) <- sub(".+?-", "",colnames(df_normal_sb))
      
      
      p_mat_mg <- p_mat
      row.names(p_mat_mg) <- gsub("\\..*","",row.names(p_mat_mg))
      print("Generate Annotation information.")
      exp_info <- AnnotationDbi::select(org.Hs.eg.db, keys = row.names(p_mat_mg), columns = c("ENTREZID", "SYMBOL"), keytype = "ENSEMBL")
      exp_info <- dplyr::distinct(exp_info)
      exp_info <- ddply(exp_info, .(ENSEMBL, ENTREZID), summarize,
                        SYMBOL=paste(SYMBOL,collapse="/"))
      df_normal_sb <- merge(df_normal_sb, exp_info, by.x=0, by.y="ENSEMBL")
      
      
      p_mat_mg <- merge(p_mat_mg, exp_info, by.x=0, by.y="ENSEMBL" )
      colnames(p_mat_mg) <- substring(colnames(p_mat_mg),1,12)
      
      
      print(paste("Generate input file for Phensim for project", p_name, d))
      for (s in colnames(df_normal_sb)[2:(ncol(df_normal_sb)-2)]) { 
        name_sb <- paste0(p, "_",d,"_", input$CONDITION[which(input$SAMPLES == s & input$drug_name ==d)],"_",info$`TR/TE`[which(info$SAMPLES == s & info$drug_name == d)],"_",s)
        
        if(file.exists(paste0(output, name_sb,"_","enrich_file.txt"))){
          print(paste("Input files for", s, "exist!"))
        } else{
          print(paste("Start - input file for Phensim:", s))
          colonne <- c("Row.names", s, "ENTREZID", "SYMBOL")
          no_exp_i <- p_mat_mg[,colonne]
          no_exp_i <- subset(no_exp_i, no_exp_i[,s] < 10)
          sb <- df_normal_sb[,colonne]
          sb <- subset(sb, sb[,s] == "UP_REG" | sb[,s] == "DOWN_REG" )
          sb[,5] <- ifelse(sb[,s] == "UP_REG", "ACTIVATION", "INHIBITION")
          
          print(paste("###",s, "generation enrich patient file"))
          enrich_file <- matrix(nrow = nrow(sb), ncol = 8)
          enrich_file[,1] <- "P"
          enrich_file[,2] <- s
          enrich_file[,3] <- "PATIENT"
          enrich_file[,4] <- sb[,3]
          enrich_file[,5] <- sb[,4]
          enrich_file[,6] <- "GENE"
          enrich_file[,7] <- "PATIENT_EDGE"
          enrich_file[,8] <- sb[,5]
          enrich_file <- enrich_file[!is.na(enrich_file[,4]),]
          colnames(enrich_file) <- c("X1", "X2", "X3", "X4", "X5", "X6","X7", "X8")
          enrich_file <- as.data.frame(enrich_file)
          enrich_file$X9 <- unique(input$CONDITION[which(input$SAMPLES == s & input$drug_name == d)])
          enrich_file <- enrich_file[!is.na(enrich_file$X4),]
          
          print(paste("###",s, "generation no-expression genes file"))
          no_exp <- no_exp_i[!(no_exp_i$ENTREZID %in% enrich_file$X4),]
          no_exp <- no_exp$ENTREZID
          no_exp <- no_exp[complete.cases(no_exp)]
          
          #name_sb <- paste0(p, "_", enrich_file[1,9],"_",s)
          write.table(enrich_file, file = paste0(output, name_sb,"_","enrich_file.txt"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
          write.table(no_exp, file = paste0(output,name_sb,"_","no_exp_genes.txt"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
          print(paste("Finish - input file for Phensim:", s))
        }
      }
    }
  }
  print("Generation simulation parameters file")
  sim_par <- matrix(ncol = 2, nrow = 1 )
  sim_par[,1] <- "P"
  sim_par[,2] <- "OVEREXPRESSION"
  
  print("Generation custom node file")
  cust_nd <- matrix(ncol = 2, nrow = 1)
  cust_nd[,1] <- "PATIENT"
  cust_nd[,2] <- "0"
  
  print("Generation custom edge file")
  cust_ed <- matrix(ncol = 1, nrow = 1)
  cust_ed[,1] <- c("PATIENT_EDGE")
  
  write.table(sim_par, file = paste0(output,"sim_par.txt"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(cust_nd, file = paste0(output, "custom_node_file.txt"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(cust_ed, file = paste0(output,"custom_edge_file.txt"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
}

drug_sim(sample_info="./ALLsamples_info_downloaded.rds",m_expr = "./expression_matrix/", output="./input_PHENSIM/", m_comb = "./comb_matrix/", strat = "./gene_stratified/")
sample_info="./ALLsamples_info_downloaded.rds"
m_expr = "./expression_matrix/"
output="./input_PHENSIM/"
m_comb = "./comb_matrix/"
strat = "./gene_stratified/"