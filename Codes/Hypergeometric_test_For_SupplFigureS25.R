###-------------------------Hypergeometric test for gene overlap

#Hypergeomtric tests between pairs of strategies

##------------------------------------------------------------------------------
#Finally, we identified 615 genes that were mapped from the 267 loci, 
#including 287 genes involved through positional mapping, 
#123 genes involved through eQTL mapping, 
#and 345 genes involved through chromatin interaction mapping. 
##------------------------------------------------------------------------------

##Assuming there were 25000 background genes

#--Position vs eQTL
phyper(48-1,180+59,25000-(180+59),116,lower.tail = F)
#phyper(48-1,116,20000-(116),287,lower.tail = F)

#--Position vs Hi-C
phyper(66-1,180+41,25000-(180+41),345, lower.tail = F)

#--eQTL vs Hi-C
phyper(17-1,58+41,25000-(58+41),345, lower.tail = F)

