#STEP1:DATAPREPROCESSING
library(TCGAbiolinks)
library(limma)
library(edgeR)
library(org.Hs.eg.db)
library(data.table)
library(plyr)
library(stringr)

normal_ctrl <- function(table,output){
  sample_info <- readRDS(table)
  project <- unique(sample_info$PROJECT)
  
  for (p in project) {
    p_name <- paste0("TCGA-", p)
    
    if(file.exists(paste0(output,p_name,"expr_matrix_normal.rds"))){
      print(paste("Expression Solid Normal Tissue matrix exists for", p_name))
    } else {
      print(paste("Download Solid Normal Tissue expression matrix for", p_name))
      query <- GDCquery(
        project = p_name,
        data.category = "Transcriptome Profiling", 
        access = "open", 
        legacy = FALSE, 
        data.type = "Gene Expression Quantification", 
        sample.type = "Solid Tissue Normal"
      )
      # Setting the working directory is useless since the directory structure of the repository starts at the root
      setwd("../../../../../")
      GDCdownload(query) 
      expr <- GDCprepare(query)
      
      # Same as above, no working directory should be set
      setwd("./PhD Sistemi Complessi/Multilayer_microbiome_analysis_Ferro/Phensim_microbiome/Pipeline_TCGA_COAD_Phensim/ALLsamples_forUMATH/")
      ex_mat <- TCGAanalyze_Preprocessing(object = expr, cor.cut = 0.6)
      
      saveRDS(ex_mat, paste0(output, p_name, "expr_matrix_normal.rds")) 
      print("Next project!")
    }
  }
}


tumoral_samples <- function(table, output){
  input <- readRDS(table)
  remove <- c()
  project <- as.vector(unique(input$PROJECT))
  drug <- unique(input$drug_name)
  for (p in project) {
    p_name <- paste0("TCGA-", p)
    for (d in drug) {
      code <- unique(input$SAMPLES[which(input$PROJECT == p & input$drug_name == d)],)
      print(paste("the number of tumor sample for project", p_name, d, "is:", length(code)))
      if(length(code) < 18){
        print(paste("NO ENOUGH SAMPLES FOR ", p_name, d))
        n <- which(input$PROJECT == p & input$drug_name == d)
        remove <- c(remove, n)
        
      } 
      else if (length(code) >= 18) {
        resp <- input$CONDITION[which(input$PROJECT == p & input$drug_name == d)]
        if(length(resp[which(resp== "RESPONDER")]) >= 4 & length(resp[which(resp== "NON_RESPONDER")]) >= 4){
          print(paste("ENOUGH NUMBER OF SAMPLES FOR", p_name, d))
          if(file.exists(paste0(output,p_name,"_",d,"_expr_matrix_tumor.rds"))){
            print(paste("Expression Tumoral matrix exists for", p_name, d))
          } else {
            print(paste("Download Tumoral expression matrix for", p_name, d))
            query <- GDCquery(
              project = p_name, 
              data.category = "Transcriptome Profiling", 
              access = "open", 
              legacy = FALSE, 
              data.type = "Gene Expression Quantification", 
              sample.type = "Primary Tumor",
              barcode= code
            )
            setwd("../../../../../")
            GDCdownload(query) 
            expr_p <- GDCprepare(query)
            ##CONTROLLARE NOMI DIRECTORY
            setwd("./PhD Sistemi Complessi/Multilayer_microbiome_analysis_Ferro/Phensim_Ruppin/Pipeline_TCGA_COAD_Phensim/ALLsamples_forUMAP_v2/")
            p_mat <- TCGAanalyze_Preprocessing(object = expr_p,
                                               cor.cut = 0.6   
            )
            #s <- unique(substring(colnames(p_mat),1,12))
            #r <- pmat
            
            saveRDS(p_mat, paste0(output,p_name,"_",d,"_expr_matrix_tumor.rds"))
            print("Next project!") 
          }
        } else { 
          print(paste("NO ENOUGH SAMPLES FOR ", p_name, d))
          m <- which(input$PROJECT == p & input$drug_name == d)
          remove <- c(remove, m)
        }
      } 
      else if (p == "STAD" & d == "Oxaliplatin"){
        print(paste("EXCEPTION FOR", p_name, d)) 
        if(file.exists(paste0(output,p_name,"_",d,"_expr_matrix_tumor.rds"))){
          print(paste("Expression Tumoral matrix exists for", p_name, d))
        } else {
          print(paste("Download Tumoral expression matrix for", p_name, d))
          query <- GDCquery(
            project = p_name, 
            data.category = "Transcriptome Profiling", 
            access = "open", 
            legacy = FALSE, 
            data.type = "Gene Expression Quantification", 
            sample.type = "Primary Tumor",
            barcode= code
          )
          setwd("../../../../../")
          GDCdownload(query) 
          expr_p <- GDCprepare(query)
          ##CONTROLLARE NOMI DIRECTORY
          setwd("./PhD Sistemi Complessi/Multilayer_microbiome_analysis_Ferro/Phensim_Ruppin/Pipeline_TCGA_COAD_Phensim/ALLsamples_forUMAP_v2/")
          p_mat <- TCGAanalyze_Preprocessing(object = expr_p,
                                             cor.cut = 0.6   
          )
          
          saveRDS(p_mat, paste0(output,p_name,"_",d,"_expr_matrix_tumor.rds"))
          print("Next project!") 
        }  
      } else {
        print("There is something wrong!!")
      }
    }
  }
  df <- input[-(remove),]
  saveRDS(df, "./ALLsamples_info_downloaded.rds")
}

normal_ctrl(table = "./ALLsamples_info.rds", output = "./expression_matrix/")
tumoral_samples(table = "./ALLsamples_info.rds", output = "./expression_matrix/")
