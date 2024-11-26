
##install DESEq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2", force = TRUE)
install.packages("xlsx")

##load library
library(DESeq2)
library(phyloseq)
library("xlsx")

##load phyloseq object non-rarified
ITS_filtered

##subset data 1
(da1<-subset_samples(ITS_filtered, G_T_C%in% c("0_B", "1_SP_H_RH")))
metadata_da1 <- as(sample_data(da1), "data.frame")

# create deseq2 object 1
deseq_da1 <- phyloseq_to_deseq2(subset_samples(da1, G_T_C%in% c("0_B", "1_SP_H_RH")), ~ G_T_C)
deseq_da1<- DESeq(deseq_da1, test="Wald", sfType="poscounts",  fitType="local")

# make pairwise comparsions 1
res_da1 = results(deseq_da1, pAdjustMethod = "BH", contrast=c("G_T_C", "0_B", "1_SP_H_RH"), cooksCutoff = FALSE)
alpha = 0.05

sigtab_da1 = res_da1[which(res_da1$padj< alpha), ]
sigtab_da1 = sigtab_da1[which(sigtab_da1$log2FoldChange < 0), ]
sigtab_da1 = cbind(as(sigtab_da1, "data.frame"), as(phyloseq::tax_table(da1)[rownames(sigtab_da1), ], "matrix"))
view(sigtab_da1)
write.csv(sigtab_da1, "D:/Botrytis stress/ITS amplicon/sigtab_da1.csv", row.names=TRUE)

##subset data 2
(da2<-subset_samples(ITS_filtered, G_T_C%in% c("0_B", "2_SL_H_RH")))
metadata_da2 <- as(sample_data(da2), "data.frame")

# create deseq2 object 2
deseq_da2 <- phyloseq_to_deseq2(subset_samples(da2, G_T_C%in% c("0_B", "2_SL_H_RH")), ~ G_T_C)
deseq_da2<- DESeq(deseq_da2, test="Wald", sfType="poscounts",  fitType="local")

# make pairwise comparsions 2
res_da2 = results(deseq_da2, pAdjustMethod = "BH", contrast=c("G_T_C", "0_B", "2_SL_H_RH"), cooksCutoff = FALSE)
alpha = 0.05
sigtab_da2 = res_da2[which(res_da2$padj< alpha), ]
sigtab_da2 = sigtab_da2[which(sigtab_da2$log2FoldChange < 0), ]
sigtab_da2 = cbind(as(sigtab_da2, "data.frame"), as(phyloseq::tax_table(da2)[rownames(sigtab_da2), ], "matrix"))
view(sigtab_da2)
write.csv(sigtab_da2, "D:/Botrytis stress/ITS amplicon/sigtab_da2.csv", row.names=TRUE)

##subset data 3
(da3<-subset_samples(ITS_filtered, G_T_C%in% c("0_B", "3_SP_H_RP")))
metadata_da3 <- as(sample_data(da3), "data.frame")

# create deseq2 object 3
deseq_da3 <- phyloseq_to_deseq2(subset_samples(da3, G_T_C%in% c("0_B", "3_SP_H_RP")), ~ G_T_C)
deseq_da3 <- DESeq(deseq_da3, test="Wald", sfType="poscounts",  fitType="local")

# make pairwise comparsions 3
res_da3 = results(deseq_da3, pAdjustMethod = "BH", contrast=c("G_T_C", "0_B", "3_SP_H_RP"), cooksCutoff = FALSE)
alpha = 0.05
sigtab_da3 = res_da3[which(res_da3$padj< alpha), ]
sigtab_da3 = sigtab_da3[which(sigtab_da3$log2FoldChange < 0), ]
sigtab_da3 = cbind(as(sigtab_da3, "data.frame"), as(phyloseq::tax_table(da3)[rownames(sigtab_da3), ], "matrix"))
view(sigtab_da3)
write.csv(sigtab_da3, "D:/Botrytis stress/ITS amplicon/sigtab_da3.csv", row.names=TRUE)

##subset data 4
(da4<-subset_samples(ITS_filtered, G_T_C%in% c("0_B", "4_SL_H_RP")))
metadata_da4 <- as(sample_data(da4), "data.frame")

# create deseq2 object 4
deseq_da4 <- phyloseq_to_deseq2(subset_samples(da4, G_T_C%in% c("0_B", "4_SL_H_RP")), ~ G_T_C)
deseq_da4 <- DESeq(deseq_da4, test="Wald", sfType="poscounts",  fitType="local")

# make pairwise comparsions 4
res_da4 = results(deseq_da4, pAdjustMethod = "BH", contrast=c("G_T_C", "0_B", "4_SL_H_RP"), cooksCutoff = FALSE)
alpha = 0.05
sigtab_da4 = res_da4[which(res_da4$padj< alpha), ]
sigtab_da4 = sigtab_da4[which(sigtab_da4$log2FoldChange < 0), ]
sigtab_da4 = cbind(as(sigtab_da4, "data.frame"), as(phyloseq::tax_table(da4)[rownames(sigtab_da4), ], "matrix"))
view(sigtab_da4)
write.csv(sigtab_da4, "D:/Botrytis stress/ITS amplicon/sigtab_da4.csv", row.names=TRUE)

##subset data 5
(da5<-subset_samples(ITS_filtered, G_T_C%in% c("0_B", "5_SP_S_RH")))
metadata_da5 <- as(sample_data(da5), "data.frame")

# create deseq2 object 5
deseq_da5 <- phyloseq_to_deseq2(subset_samples(da5, G_T_C%in% c("0_B", "5_SP_S_RH")), ~ G_T_C)
deseq_da5 <- DESeq(deseq_da5, test="Wald", sfType="poscounts",  fitType="local")

# make pairwise comparsions 5
res_da5 = results(deseq_da5, pAdjustMethod = "BH", contrast=c("G_T_C", "0_B", "5_SP_S_RH"), cooksCutoff = FALSE)
alpha = 0.05
sigtab_da5 = res_da5[which(res_da5$padj< alpha), ]
sigtab_da5 = sigtab_da5[which(sigtab_da5$log2FoldChange < 0), ]
sigtab_da5 = cbind(as(sigtab_da5, "data.frame"), as(phyloseq::tax_table(da5)[rownames(sigtab_da5), ], "matrix"))
view(sigtab_da5)
write.csv(sigtab_da5, "D:/Botrytis stress/ITS amplicon/sigtab_da5.csv", row.names=TRUE)

##subset data 6
(da6<-subset_samples(ITS_filtered, G_T_C%in% c("0_B", "6_SL_S_RH")))
metadata_da6 <- as(sample_data(da6), "data.frame")

# create deseq2 object 6
deseq_da6 <- phyloseq_to_deseq2(subset_samples(da6, G_T_C%in% c("0_B", "6_SL_S_RH")), ~ G_T_C)
deseq_da6 <- DESeq(deseq_da6, test="Wald", sfType="poscounts",  fitType="local")

# make pairwise comparison 6
res_da6 = results(deseq_da6, pAdjustMethod = "BH", contrast=c("G_T_C", "0_B", "6_SL_S_RH"), cooksCutoff = FALSE)
alpha = 0.05
sigtab_da6 = res_da6[which(res_da6$padj< alpha), ]
sigtab_da6 = sigtab_da6[which(sigtab_da6$log2FoldChange < 0), ]
sigtab_da6 = cbind(as(sigtab_da6, "data.frame"), as(phyloseq::tax_table(da6)[rownames(sigtab_da6), ], "matrix"))
view(sigtab_da6)
write.csv(sigtab_da6, "D:/Botrytis stress/ITS amplicon/sigtab_da6.csv", row.names=TRUE)

##subset data 7
(da7<-subset_samples(ITS_filtered, G_T_C%in% c("0_B", "7_SP_S_RP")))
metadata_da7 <- as(sample_data(da7), "data.frame")

# create deseq2 object 7
deseq_da7 <- phyloseq_to_deseq2(subset_samples(da7, G_T_C%in% c("0_B", "7_SP_S_RP")), ~ G_T_C)
deseq_da7 <- DESeq(deseq_da7, test="Wald", sfType="poscounts",  fitType="local")

# make pairwise comparsions 7
res_da7 = results(deseq_da7, pAdjustMethod = "BH", contrast=c("G_T_C", "0_B", "7_SP_S_RP"), cooksCutoff = FALSE)
alpha = 0.05
sigtab_da7 = res_da7[which(res_da7$padj< alpha), ]
sigtab_da7 = sigtab_da7[which(sigtab_da7$log2FoldChange < 0), ]
sigtab_da7 = cbind(as(sigtab_da7, "data.frame"), as(phyloseq::tax_table(da7)[rownames(sigtab_da7), ], "matrix"))
view(sigtab_da7)
write.csv(sigtab_da7, "D:/Botrytis stress/ITS amplicon/sigtab_da7.csv", row.names=TRUE)

##subset data 8
(da8<-subset_samples(ITS_filtered, G_T_C%in% c("0_B", "8_SL_S_RP")))
metadata_da8 <- as(sample_data(da8), "data.frame")

# create deseq2 object 8
deseq_da8 <- phyloseq_to_deseq2(subset_samples(da8, G_T_C%in% c("0_B", "8_SL_S_RP")), ~ G_T_C)
deseq_da8 <- DESeq(deseq_da8, test="Wald", sfType="poscounts",  fitType="local")

# make pairwise comparsions 7
res_da8 = results(deseq_da8, pAdjustMethod = "BH", contrast=c("G_T_C", "0_B", "8_SL_S_RP"), cooksCutoff = FALSE)
alpha = 0.05
sigtab_da8 = res_da8[which(res_da8$padj< alpha), ]
sigtab_da8 = sigtab_da8[which(sigtab_da8$log2FoldChange < 0), ]
sigtab_da8 = cbind(as(sigtab_da8, "data.frame"), as(phyloseq::tax_table(da8)[rownames(sigtab_da8), ], "matrix"))
view(sigtab_da8)

# Extract row names from each data frame
row_names_da1 <- row.names(sigtab_da1)
row_names_da2 <- row.names(sigtab_da2)
row_names_da3 <- row.names(sigtab_da3)
row_names_da4 <- row.names(sigtab_da4)
row_names_da5 <- row.names(sigtab_da5)
row_names_da6 <- row.names(sigtab_da6)
row_names_da7 <- row.names(sigtab_da7)
row_names_da8 <- row.names(sigtab_da8)

# Combine row names into a single vector and get unique names
all_row_da <- unique(c(row_names_da1, row_names_da2, row_names_da3, row_names_da4, row_names_da5,
                       row_names_da6, row_names_da7, row_names_da8))

##prune taxa not present/significantly higher than in bulk (unplanted soil)
significant_asv <- c("ASV4" , "ASV9" , "ASV21" , "ASV27" , "ASV39" , "ASV45"  ,"ASV59"  ,"ASV61",  "ASV62"  ,"ASV72" , "ASV87" , "ASV88",
                     "ASV116", "ASV127" ,"ASV134", "ASV141", "ASV162" ,"ASV164" ,"ASV171" ,"ASV186", "ASV192", "ASV195" ,"ASV201" ,"ASV210",
                     "ASV214" , "ASV219" , "ASV226" , "ASV394" , "ASV444", "ASV68" , "ASV91",  "ASV93" , "ASV101" ,"ASV150", "ASV187" ,"ASV193",
                     "ASV198", "ASV249", "ASV289", "ASV42" , "ASV89" , "ASV94" , "ASV96" , "ASV103", "ASV108" , "ASV135" , "ASV143" ,"ASV153",
                     "ASV173", "ASV191", "ASV199" ,"ASV253", "ASV258", "ASV275", "ASV276", "ASV290","ASV319", "ASV322" ,"ASV359", "ASV374",
                     "ASV406", "ASV446", "ASV52",  "ASV69" , "ASV86" , "ASV128", "ASV133", "ASV138", "ASV174", "ASV183" ,"ASV224" ,"ASV267",
                     "ASV301", "ASV344", "ASV427", "ASV470" ,"ASV73" , "ASV188", "ASV235" ,"ASV251",  "ASV259" , "ASV280", "ASV15" , "ASV47",
                     "ASV67",  "ASV74" , "ASV145", "ASV157", "ASV202", "ASV204", "ASV231","ASV269" ,"ASV287", "ASV34" , "ASV248" ,"ASV11", 
                     "ASV84" , "ASV381" , "ASV423")
ITS_significant <-prune_taxa(significant_asv,ITS_filtered )
extract_otu_ITS<-otu_table(ITS_significant)
write.csv(extract_otu_ITS, file = "OTU_ITS.csv")
extract_tax_ITS<- phyloseq::tax_table (ITS_significant)
write.csv(extract_tax_ITS, file = "TAX_ITS.csv")

