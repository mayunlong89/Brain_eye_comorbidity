

### This is for H-MGAMA-based permutation analysis

###---------------10000 times of Monte Carlo permutation analysis-------------------------

setwd("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\04-metasoft\\09-H-MAGMA-S-MultiXcan\\MTAG GWAS result_H-MAGMA\\")

eye_brain_plei <- read.table("91_eye_brain_pleiotropic_genes.txt",header = TRUE)
plei_genes <- eye_brain_plei$ID 


##Function for reading all files in a directory
getwd_cd <- getwd()

#read file name
first_file_name <- list.files()
first_file_name <- first_file_name[-c(1,5)]
dir <- paste(getwd_cd,"/",first_file_name,sep = "")
n <- length(dir)

n_sub <- rep(0,n)
n_sub <- as.data.frame(n_sub)
n_sub <- t(n_sub)

#read file
b <- list.files(dir[1])
m <- length(b)

#diseas name file
name_col <- strsplit(b,"_Adult_brain.genes.out")
name_col <- unlist(name_col)


#--read all files to generate a list file
a <- vector("list",length = n)

for (i in 1:n){
 
 
   b <- list.files(dir[i])
   m <- length(b)
   
  temp <- as.data.frame(matrix(nrow = 0, ncol = 0))
  
  if(ncol(temp)==0){
    
    file = paste(dir[i],"/",b[1],sep="")
    all_data <- read.csv(file,sep = "\t")
    temp <- as.data.frame(all_data$GENE)
    colnames(temp) <- "GENE"
  } 
    
    for (j in 1:m){
      
      file = paste(dir[i],"/",b[j],sep="")
      
      all_data <- read.csv(file,sep = "\t")
      
      temp1 <- all_data[,c(1,9)]
      
      ###utilizing merge function to combine various datasets with different columns
      temp <- merge(temp,temp1,by="GENE",all=F)
    }
  
  #name of the data.frame for each cell type or tissue type
  colnames(temp) <- c("gene_name",name_col) 
  a[[i]] <- temp
}
##Name of this list file
names(a) <- first_file_name
head(a[[1]])


a[[1]]$gene_name


##---Take all genes from one trait as background genes
#data_adult_brain_ADHD <- read.csv("Adult_brain/ADHD_Adult_brain.genes.out",sep = "\t",header=T)
#Among 91 genes, one gene (GNL3L	ENSG00000130119) did not have annotations
#So we only used 90 genes as input for permutation analysis
##removing these 90 eye-brain pleiotropic genes from the background
#gene_background <- data_adult_brain_ADHD$GENE[-which(data_adult_brain_ADHD$GENE %in% plei_genes)]


library(dplyr)
names(a)  ###The list format datasets for all single-cell data
z=6  ### z = 1:length(a) element number indicating which dataset in a is planned to be used 
gene_background <- a[[z]]$gene_name[-which(a[[z]]$gene_name %in% plei_genes)]

set.seed(12345)

##construct a function for permutation

#MC permutation analysis
MC_test <- function(all_bkgenes,num_ran,threshold){
  
rand_selected_genes <- sample(gene_background,num_ran)

temp_ran_data <- a[[z]][which(a[[z]]$gene_name %in% rand_selected_genes),]

temp_ran_data1 <- temp_ran_data[,-1]
rownames(temp_ran_data1)<-temp_ran_data[,1]
sum_count <- as.data.frame(rowSums(as.matrix(temp_ran_data1 < threshold)))
num <- count(sum_count,sum_count[,1]>1) #counting the number of genes with two significant hits
#num <- hist(sum_count[,1],plot=FALSE)$counts
num <- num[2,2]

  return(num)
  
}


# conducting 1000 times permutation analysis for each tissue or cell type
all_bkgenes = gene_background
num_ran = 90
threshold = 0.00055  #(0.05/90)
results_permut <- replicate(10000,MC_test(all_bkgenes,num_ran,threshold))



#Ploting function

count_obs <- c(45,44,44,41,44,48)
names(a)
count_obs2 <- data.frame(names(a),count_obs)

#set the name for each plot
namexx <- paste(names(a)[z],"Monte Carlo test",sep = ": ")

#P-value
P_value = length(results_permut[results_permut>count_obs2[z,2]])/length(results_permut)
P_value
#Visualization
hist(results_permut, col="#D2691E",xlab="Counts of overlapped genes",
     xlim = c(0,50),
     main= namexx)
     #main=NULL)
abline(v=count_obs2[z,2],col="blue",lty="longdash")


