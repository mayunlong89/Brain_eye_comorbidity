###----------COSG method identifies cell type-specific genes

#COSG_top100_general <- read.table("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\05-singlecell\\COSG_celltype_top100.txt",header=TRUE)
COSG_top500_general <- read.table("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\05-singlecell\\COSG_celltype_top500.txt",header=TRUE)
pleiotropic_91genes <- read.table("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\05-singlecell\\91_pleiotropic_genes.txt",header=TRUE)
nonpleiotropic_58genes <- read.table("E:\\00-新课题-总结\\15-眼脑交叉学科建设\\00-Solid paper-for 定平\\01-Data-analysis\\03-All_analysis\\05-singlecell\\58_nonpleiotropic_genes.txt",header=TRUE)

#n1=length(COSG_top100_general$names.Microglia)
n2=length(COSG_top500_general$names.Microglia)
m <- length(pleiotropic_91genes$eye_brain_pleiotropic_genes)
m2 <-length(nonpleiotropic_58genes$nonpleiotropic_genes)
names(COSG_top500_general)

##pleiotropic genes
j <- 7
All_bkgenes <- 25000
common <- intersect(pleiotropic_91genes[,1],COSG_top500_general[,j])
n <- length(common)
phyper(n-1,m-n,All_bkgenes-(m-n),n2,lower.tail = F)
common
names(COSG_top500_general)[j]

#non-pleiotropic genes
j <- 8
All_bkgenes <- 25000
common <- intersect(nonpleiotropic_58genes[,1],COSG_top500_general[,j])
n <- length(common)
phyper(n-1,m2-n,All_bkgenes-(m2-n),n2,lower.tail = F)
common
names(COSG_top500_general)[j]


