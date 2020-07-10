setwd("~/Downloads")
library(ggplot2)
library(reshape2)
library(amap)

table1=read.table("DE_CC_vs_HC.csv",header=T,row.names = 1)
table2=read.table("DE_FLU_vs_CC.csv",header=T,row.names = 1)
table3=read.table("DE_FLU_vs_HC.csv",header=T,row.names = 1)
table4=read.table("EM-2.csv",header=T,row.names = 1)
table5=read.table("sample_sheet.csv",header=T,row.names = 1)
table6=read.table("mart_export-2.txt",header=T,row.names = 1,sep = '\t')
annotated_expression_matrix=merge(table4,table6,by=0)

## Rename column Gene.name
names(annotated_expression_matrix)[names(annotated_expression_matrix) == "Gene.name"] <- "gene_symbol"

##mean
annotated_expression_matrix$mean = apply(annotated_expression_matrix[,c(2:43)],1,mean)

## annotated DE tables    
annotated_table1=merge(table1,table6,by=0)
annotated_table2=merge(table2,table6,by=0)
annotated_table3=merge(table3,table6,by=0)
names(annotated_table1)[names(annotated_table1) == "Gene.name"] <- "gene_symbol"
names(annotated_table2)[names(annotated_table2) == "Gene.name"] <- "gene_symbol"
names(annotated_table3)[names(annotated_table3) == "Gene.name"] <- "gene_symbol"
annotated_table1_CC_vs_HC=annotated_table1[,c(1,2,3,4,5)]
annotated_table2_FLU_vs_CC=annotated_table2[,c(1,2,3,4,5)]
annotated_table3_FLU_vs_HC=annotated_table3[,c(1,2,3,4,5)]

##-log10p calculation
annotated_table1_CC_vs_HC$log10p=-log10(annotated_table1_CC_vs_HC$p)
names(annotated_table1_CC_vs_HC)[names(annotated_table1_CC_vs_HC) == "log10p"] <- "-log10p"
annotated_table2_FLU_vs_CC$log10p=-log10(annotated_table2_FLU_vs_CC$p)
names(annotated_table2_FLU_vs_CC)[names(annotated_table2_FLU_vs_CC) == "log10p"] <- "-log10p"
annotated_table3_FLU_vs_HC$log10p=-log10(annotated_table3_FLU_vs_HC$p)
names(annotated_table3_FLU_vs_HC)[names(annotated_table3_FLU_vs_HC) == "log10p"] <- "-log10p"

##SIGNIFICANCE 
annotated_table1_CC_vs_HC$significance= as.factor(annotated_table1_CC_vs_HC$p.adj < 0.05 & abs(annotated_table1_CC_vs_HC$log2fold) > 1.0)
annotated_table2_FLU_vs_CC$significance= as.factor(annotated_table2_FLU_vs_CC$p.adj < 0.05 & abs(annotated_table2_FLU_vs_CC$log2fold) > 1.0)
annotated_table3_FLU_vs_HC$significance= as.factor(annotated_table3_FLU_vs_HC$p.adj < 0.05 & abs(annotated_table3_FLU_vs_HC$log2fold) > 1.0)

########################################################################################################################################################################################################################################################################################################################################################
#organising all the information in a master table

##merging all tables
merged_table1=merge(annotated_table1_CC_vs_HC,annotated_table2_FLU_vs_CC,by=0)
merged_table1=merged_table1[c(2:8,10:15)]

names(merged_table1)[names(merged_table1) == "Row.names.x"] <- "Row.names" 
names(merged_table1)[names(merged_table1) == "gene_symbol.x"] <- "gene_symbol"
names(merged_table1)[names(merged_table1) == "log2fold.x"] <- "log2fold.CC_vs_HC"
names(merged_table1)[names(merged_table1) == "p.x"] <- "p.CC_vs_HC"
names(merged_table1)[names(merged_table1) == "p.adj.x"] <- "p.adj.CC_vs_HC"
names(merged_table1)[names(merged_table1) == "-log10p.x"] <- "-log10p.CC_vs_HC"
names(merged_table1)[names(merged_table1) == "significance.x"] <- "significance.CC_vs_HC"
names(merged_table1)[names(merged_table1) == "log2fold.y"] <- "log2fold.FLU_vs_CC"
names(merged_table1)[names(merged_table1) == "p.y"] <- "p.FLU_vs_CC"
names(merged_table1)[names(merged_table1) == "-log10p.y"] <- "-log10p.FLU_vs_CC"
names(merged_table1)[names(merged_table1) == "significance.y"] <- "significance.FLU_vs_CC"
names(merged_table1)[names(merged_table1) == "p.adj.y"] <- "p.adj.FLU_vs_CC"

merged_table2=merge(merged_table1,annotated_table3_FLU_vs_HC,by=0)
names(merged_table2)[names(merged_table2) == "log2fold"] <- "log2fold.FLU_vs_HC"
names(merged_table2)[names(merged_table2) == "p"] <- "p.FLU_vs_HC"
names(merged_table2)[names(merged_table2) == "-log10p"] <- "-log10p.FLU_vs_HC"
names(merged_table2)[names(merged_table2) == "significance"] <- "significance.FLU_vs_HC"
names(merged_table2)[names(merged_table2) == "p.adj"] <- "p.adj.FLU_vs_HC"
merged_table2=merged_table2[c(2:11,13:14,16:18,20:21)]
names(merged_table2)[names(merged_table2) == "Row.names.x"] <- "Row.names"
names(merged_table2)[names(merged_table2) == "gene_symbol.x"] <- "gene_symbol"

##making a master table
demo_master=merge(merged_table2,annotated_expression_matrix,by=0)
demo_master=demo_master[c(2:18,20:61,63:67)]
names(demo_master)[names(demo_master) == "Row.names.x"] <- "Row.names"
names(demo_master)[names(demo_master) == "gene_symbol.x"] <- "gene_symbol"
master_table=na.omit(demo_master)
master_table=master_table[c(1,5,2:4,6:59,64,60:63)]
row.names(master_table)=master_table$Row.names
master_table=master_table[,c(2:64)]
########################################################################################################################################################################################################################################################################################################################################################

##DIFFERENTIAL EXPRESSION TABLE
Diffexp.CC_vs_HC <- subset(master_table, significance.CC_vs_HC==TRUE | log2fold.CC_vs_HC >0 ,select=c(gene_symbol,log2fold.CC_vs_HC,p.CC_vs_HC,p.adj.CC_vs_HC,CC_1,CC_2,CC_3,CC_4,CC_5,CC_6,CC_7,CC_8,CC_9,CC_10,CC_11,CC_12,CC_13,CC_14,HC_1,HC_2,HC_3,HC_4,HC_5,HC_6,HC_7,HC_8,HC_9,HC_10,HC_11,HC_12,HC_13,HC_14,mean))
Diffexp.FLU_vs_CC <- subset(master_table, significance.FLU_vs_CC==TRUE | log2fold.FLU_vs_CC >0 ,select=c(gene_symbol,log2fold.FLU_vs_CC,p.FLU_vs_CC,p.adj.FLU_vs_CC,CC_1,CC_2,CC_3,CC_4,CC_5,CC_6,CC_7,CC_8,CC_9,CC_10,CC_11,CC_12,CC_13,CC_14,FLU_1,FLU_2,FLU_3,FLU_4,FLU_5,FLU_6,FLU_7,FLU_8,FLU_9,FLU_10,FLU_11,FLU_12,FLU_13,FLU_14,mean))
Diffexp.FLU_vs_HC <- subset(master_table, significance.FLU_vs_HC==TRUE | log2fold.FLU_vs_HC >0 ,select=c(gene_symbol,log2fold.FLU_vs_HC,p.FLU_vs_HC,p.adj.FLU_vs_HC,HC_1,HC_2,HC_3,HC_4,HC_5,HC_6,HC_7,HC_8,HC_9,HC_10,HC_11,HC_12,HC_13,HC_14,FLU_1,FLU_2,FLU_3,FLU_4,FLU_5,FLU_6,FLU_7,FLU_8,FLU_9,FLU_10,FLU_11,FLU_12,FLU_13,FLU_14,mean))                           
##save the tables
write.csv(master_table,file = "master_table.csv")
write.csv(Diffexp.CC_vs_HC,file = "Diffexp.CC_vs_HC.csv")
write.csv(Diffexp.FLU_vs_CC,file = "Diffexp.FLU_vs_CC.csv")
write.csv(Diffexp.FLU_vs_HC,file = "Diffexp.FLU_vs_HC.csv")

########################################################################################################################################################################################################################################################################################################################################################
#PLOTS

##MA plot
ggplot(data=master_table) + geom_point(aes(x=log10(mean), y= `log2fold.CC_vs_HC`,color=significance.CC_vs_HC))
ggplot(data=master_table) + geom_point(aes(x=log10(mean), y= `log2fold.FLU_vs_CC`,color=significance.FLU_vs_CC))
ggplot(data=master_table) + geom_point(aes(x=log10(mean), y= `log2fold.FLU_vs_HC`,color=significance.FLU_vs_HC))

##density plot
ggplot(master_table,aes(x=log10(mean)))+geom_density(colour="black",fill="brown")+labs(title="Expression Density",x="Mean Log10",y="Density")

##volcano plot
ggplot(data=master_table) + geom_point(aes(x=`log2fold.CC_vs_HC`, y= `-log10p.CC_vs_HC`,color=significance.CC_vs_HC))
ggplot(data=master_table) + geom_point(aes(x=`log2fold.FLU_vs_CC`, y= `-log10p.FLU_vs_CC`,color=significance.FLU_vs_CC))
ggplot(data=master_table) + geom_point(aes(x=`log2fold.FLU_vs_HC`, y= `-log10p.FLU_vs_HC`,color=significance.FLU_vs_HC))

##SORTED P  VALUE
sorted_CC_vs_HC=master_table[order(master_table$p.CC_vs_HC),]
sorted_FLU_vs_CC=master_table[order(master_table$p.FLU_vs_CC),]
sorted_FLU_vs_HC=master_table[order(master_table$p.FLU_vs_HC),]
##PCA
my_data=as.matrix(sapply(table4, as.numeric))
pca = prcomp(t(my_data))
pca_coordinates = data.frame(pca$x)

########################################################################################################################################################################################################################################################################################################################################################
#Heatmap

##zscore table
zscore=scale(table4)
##correlation
cors = cor(zscore, method="spearman")
melted = melt(cors)
ggplot(melted, aes(x = Var1, y = Var2, fill = value)) + geom_tile()+ scale_fill_gradientn(limits=c(-1, 1),colours = colorRampPalette(colours)(100))
plot(pca$x[,1], pca$x[,2],xlab ='PC1',ylab = 'PC2',col=table5$GROUP)


##diff heatmap
CC_vs_HC=Diffexp.CC_vs_HC[c(2,6:32)]
FLU_vs_CC=Diffexp.FLU_vs_CC[c(2,6:32)]
FLU_vs_HC=Diffexp.FLU_vs_HC[c(2,6:32)]
CC_vs_HC=CC_vs_HC[c(2:28)]
FLU_vs_CC=FLU_vs_CC[c(2:28)]
FLU_vs_HC=FLU_vs_HC[c(2:28)]
Scaled_FLU_vs_HC=scale(FLU_vs_HC)
Scaled_FLU_vs_CC=scale(FLU_vs_CC)
Scaled_CC_vs_HC=scale(CC_vs_HC)
##for CC_vs_HC
y.dist1 = Dist(Scaled_CC_vs_HC, method="spearman")
y.cluster1 = hclust(y.dist1, method="average")
y.dd1 = as.dendrogram(y.cluster1)
y.dd.reordered1 = reorder(y.dd1,0,FUN="average")
y.order1 = order.dendrogram(y.dd.reordered1)
CC_vs_HC_ordered = Scaled_CC_vs_HC[y.order1,]
melted1 = melt(as.matrix(CC_vs_HC_ordered))
heatmap(CC_vs_HC_ordered,col = terrain.colors(256))
##for FLU_vs_CC
y.dist2 = Dist(Scaled_FLU_vs_CC, method="spearman")
y.cluster2 = hclust(y.dist2, method="average")
y.dd2 = as.dendrogram(y.cluster2)
y.dd.reordered2 = reorder(y.dd2,0,FUN="average")
y.order2 = order.dendrogram(y.dd.reordered2)
FLU_vs_CC_ordered = Scaled_FLU_vs_CC[y.order2,]
melted2 = melt(as.matrix(FLU_vs_CC_ordered))
heatmap(FLU_vs_CC_ordered,col = terrain.colors(256))
##for FLU_vs_HC
y.dist3 = Dist(Scaled_FLU_vs_HC, method="spearman")
y.cluster3 = hclust(y.dist3, method="average")
y.dd3 = as.dendrogram(y.cluster3)
y.dd.reordered3 = reorder(y.dd3,0,FUN="average")
y.order3 = order.dendrogram(y.dd.reordered3)
FLU_vs_HC_ordered = Scaled_FLU_vs_HC[y.order3,]
melted3 = melt(as.matrix(FLU_vs_HC_ordered))
heatmap(FLU_vs_HC_ordered,col = terrain.colors(256))

########################################################################################################################################################################################################################################################################################################################################################



