
setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/3-lfmm")

library(gradientForest)
library(MaizePal)

gfData <- read.table("2200.GF",header =T,sep="\t",row.names = "Popu")  # 加载环境数据、等位基因频率
candidate <- gfData[,grep("chr",names(gfData))] #提取包含candSNPs的列，等位基因频率
present <- gfData[,c(1,2,grep("bio",names(gfData)))]  # 提取第一列、第二列和包含bio的列，即坐标和生物气候数据
bioclimatic <- paste("bio",1:19,sep = "")  # 生成向量(bio_1, ..., bio_19)
maxLevel <- log2(0.368*nrow(candidate)/2)# 固定公式 log2（0.368*居群数/2）
gf_candidate <- gradientForest(cbind(present[,bioclimatic], candidate),  predictor.vars=colnames(present[,bioclimatic]),response.vars=colnames(candidate), ntree=500,maxLevel=maxLevel, trace=T, corr.threshold=0.50)  #决策树默认500，相关性阈值0.5/0.7
#plot(gf_candidate, plot.type = "Overall.Importance", col=c(rep("grey",15),MaizePal::maize_pal("HighlandMAGIC", 4) ),las=2,cex.names=0.8)
bio_cand <- gf_candidate$overall.imp[order(gf_candidate$overall.imp,decreasing = T)]  

###bio4
most_cand <- names(bio_cand[1]) 

barplot(bio_cand,las=2,cex.names=0.8,col=c(MaizePal::maize_pal("JimmyRed",4),rep("grey",15)),ylab="Weighted importance")
barplot(bio_cand,las=2,cex.names=0.8,col=c(MaizePal::maize_pal("JimmyRed",5),rep("grey",14)),ylab="Weighted importance")


temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand, type=c("Overall"),standardize = T) # al candidate SNPs 
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand, type=c("Species"),standardize = T) #each individual candidate allele


pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))]) 
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) 
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand])))
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) 
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm])
plot(temp_cand_overall,type="n",ylim=c(0,0.25),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= paste("Temperature Seasonality (standard deviation ×100) (BIO4)",sep="")) 
for(j in 1:length(temp_cand_SNP)){   lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6)) } 

lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4) 
warm_col=MaizePal::maize_pal("JimmyRed",4) 
cold_col=MaizePal::maize_pal("MaizAzul",4) 
id_c <- order(temp$bio[cold]) 

id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col)) 

id_w <- order(temp$bio[warm]) 

id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col)) 

#points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg=rev(id_wcol),cex=1.5) 

#points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg=id_ccol,cex=1.5) 

points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg="#185A56",cex=1.5) 

points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg="#B89076",cex=1.5) 

rownames(pop_turn)[cold]
[1] "GER_Bonn"        "UK_Cam"          "AUS_Ehrenhausen" "AUS_Grafenegg"   "IRE_ARDARA"      "IRE_MAGHERA"     "IRE_LETTERKENNY" "IRE_Dublin"      "AUS_MiBa"        "UK_Norwich"      "UK_Oxford"       "AUS_Schubert"    "SWZ_Zurich"     
id_c
[1]  5  7  6  8 10 11  2  1 13 12  4  3  9
rownames(pop_turn)[warm]
[1] "JPN_HakubaMura"  "JPN_Suwa2"       "JPN_Guha"        "JPN_MiyaOh"      "JPN_Taimadera"   "JPN_Suwa1"       "JPN_Hiro"        "JPN_Mats"        "JPN_Takaragaike"
id_w
[1] 7 5 9 3 2 6 8 4 1




#Bio18
most_cand <- names(bio_cand[2]) 

temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand, type=c("Overall"),standardize = T) # al candidate SNPs 
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand, type=c("Species"),standardize = T) #each individual candidate allele

pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))]) 
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) 
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand])))
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) 
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm])
plot(temp_cand_overall,type="n",ylim=c(0,0.25),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= paste("Precipitation of Warmest Quarter (BIO18)",sep="")) 
for(j in 1:length(temp_cand_SNP)){   lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6)) } 

lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4) 
warm_col=MaizePal::maize_pal("JimmyRed",4) 
cold_col=MaizePal::maize_pal("MaizAzul",4) 
id_c <- order(temp$bio[cold]) 

id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col)) 

id_w <- order(temp$bio[warm]) 

id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col)) 

points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg="#185A56",cex=1.5) 

points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg="#B89076",cex=1.5) 

rownames(pop_turn)[cold]
[1] "GER_Bonn"        "UK_Cam"          "AUS_Ehrenhausen" "AUS_Grafenegg"   "IRE_ARDARA"      "IRE_MAGHERA"     "IRE_LETTERKENNY" "IRE_Dublin"      "AUS_MiBa"        "UK_Norwich"      "UK_Oxford"       "AUS_Schubert"    "SWZ_Zurich"     

id_c
2 11 10  8 12  1  9  4  6  7  3 13  5

rownames(pop_turn)[warm]
[1] "JPN_HakubaMura"  "JPN_Suwa2"       "JPN_Guha"        "JPN_MiyaOh"      "JPN_Taimadera"   "JPN_Suwa1"       "JPN_Hiro"        "JPN_Mats"        "JPN_Takaragaike"


id_w
2 6 1 4 8 5 9 7 3




## BIO5

most_cand <- names(bio_cand[3]) 
temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand, type=c("Overall"),standardize = T) # al candidate SNPs 
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand, type=c("Species"),standardize = T) #each individual candidate allele

pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))]) 
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) 
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand])))
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) 
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm])
plot(temp_cand_overall,type="n",ylim=c(0,0.25),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= paste("Max Temperature of Warmest Month (BIO5)",sep="")) 
for(j in 1:length(temp_cand_SNP)){   lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6)) } 

lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4) 
warm_col=MaizePal::maize_pal("JimmyRed",4) 
cold_col=MaizePal::maize_pal("MaizAzul",4) 
id_c <- order(temp$bio[cold]) 

id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col)) 

id_w <- order(temp$bio[warm]) 

id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col)) 

points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg="#185A56",cex=1.5) 

points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg="#B89076",cex=1.5) 

rownames(pop_turn)[cold]
[1] "GER_Bonn"        "UK_Cam"          "AUS_Ehrenhausen" "AUS_Grafenegg"   "IRE_ARDARA"      "IRE_MAGHERA"     "IRE_LETTERKENNY" "IRE_Dublin"      "AUS_MiBa"        "UK_Norwich"      "UK_Oxford"       "AUS_Schubert"    "SWZ_Zurich"     

id_c
5  7  6  8 10 11  2 13  1  3 12  4  9

rownames(pop_turn)[warm]
[1] "JPN_HakubaMura"  "JPN_Suwa2"       "JPN_Guha"        "JPN_MiyaOh"      "JPN_Taimadera"   "JPN_Suwa1"       "JPN_Hiro"        "JPN_Mats"        "JPN_Takaragaike"


id_w
8 4 1 2 6 3 7 5 9

#bio10
most_cand <- names(bio_cand[4]) 

## Allele turnover functions across the landscape. 

temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand, type=c("Overall"),standardize = T) # al candidate SNPs 
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand, type=c("Species"),standardize = T) #each individual candidate allele

pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))]) 
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) 
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand])))
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) 
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm])
plot(temp_cand_overall,type="n",ylim=c(0,0.25),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= paste("Mean Temperature of Warmest Quarter (BIO10)",sep="")) 
for(j in 1:length(temp_cand_SNP)){   lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6)) } 

lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4) 
warm_col=MaizePal::maize_pal("JimmyRed",4) 
cold_col=MaizePal::maize_pal("MaizAzul",4) 
id_c <- order(temp$bio[cold]) 

id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col)) 

id_w <- order(temp$bio[warm]) 

id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col)) 

points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg="#185A56",cex=1.5) 

points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg="#B89076",cex=1.5) 

rownames(pop_turn)[cold]
[1] "GER_Bonn"        "UK_Cam"          "AUS_Ehrenhausen" "AUS_Grafenegg"   "IRE_ARDARA"      "IRE_MAGHERA"     "IRE_LETTERKENNY" "IRE_Dublin"      "AUS_MiBa"        "UK_Norwich"      "UK_Oxford"       "AUS_Schubert"    "SWZ_Zurich"     

id_c
5  7  6  8 10 11  2 13  1  3  4  9 12

rownames(pop_turn)[warm]
[1] "JPN_HakubaMura"  "JPN_Suwa2"       "JPN_Guha"        "JPN_MiyaOh"      "JPN_Taimadera"   "JPN_Suwa1"       "JPN_Hiro"        "JPN_Mats"        "JPN_Takaragaike"


id_w
8 4 1 2 6 3 7 5 9

#BIO8
most_cand <- names(bio_cand[5]) 
temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand, type=c("Overall"),standardize = T) # al candidate SNPs 
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand, type=c("Species"),standardize = T) #each individual candidate allele

pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))]) 
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) 
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand])))
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) 
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm])
plot(temp_cand_overall,type="n",ylim=c(0,0.25),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= paste("Mean Temperature of Wettest Quarter (BIO8)",sep="")) 
for(j in 1:length(temp_cand_SNP)){   lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6)) } 

lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4) 
warm_col=MaizePal::maize_pal("JimmyRed",4) 
cold_col=MaizePal::maize_pal("MaizAzul",4) 
id_c <- order(temp$bio[cold]) 

id_ccol <- as.character(cut(1:length(id_c),length(cold_col),labels=cold_col)) 

id_w <- order(temp$bio[warm]) 

id_wcol <- as.character(cut(1:length(id_w),length(warm_col),labels=warm_col)) 

points(temp$bio[warm][id_w],temp$imp[warm][id_w],pch=21,bg="#185A56",cex=1.5) 

points(temp$bio[cold][id_c],temp$imp[cold][id_c],pch=21,bg="#B89076",cex=1.5) 

rownames(pop_turn)[cold]
[1] "GER_Bonn"        "UK_Cam"          "AUS_Ehrenhausen" "AUS_Grafenegg"   "IRE_ARDARA"      "IRE_MAGHERA"     "IRE_LETTERKENNY" "IRE_Dublin"      "AUS_MiBa"        "UK_Norwich"      "UK_Oxford"       "AUS_Schubert"    "SWZ_Zurich"     

id_c
11  6 10  7  8  5  2  1 13 12  3  4  9

rownames(pop_turn)[warm]
[1] "JPN_HakubaMura"  "JPN_Suwa2"       "JPN_Guha"        "JPN_MiyaOh"      "JPN_Taimadera"   "JPN_Suwa1"       "JPN_Hiro"        "JPN_Mats"        "JPN_Takaragaike"


id_w
8 4 1 2 6 5 3 7 9

######LAT LON corr with imp IMPORTANT#####
setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/3-lfmm/LAT_LON_corrIMP")
bio5_euro=read.table("LAT_LON_corrIMP/BIO5-Euro.tsv",header = T)

ggscatter(bio5_euro, x = "Lat", y = "y",
                 size = 2.5,
                 add = "reg.line",add.params = list(color="#9c9c9c", fill="#d6d6d6", size = 1),conf.int = TRUE)+stat_cor(method = "spearman")

df=read.table("LAT_LON_corrIMP/test2.tsv",header = T)

LAT=ggscatter(df, x = "Lat", y = "y",
                    add = "reg.line",                         
                     conf.int = TRUE,  ellipse = TRUE,                        
                     color = "Anno",palette = "lancet",facet.by = "Bio"        
           )+stat_cor(method = "spearman",aes(color = Anno), label.x = 40)

LON=ggscatter(df, x = "Lon", y = "y",
              add = "reg.line",                         
              conf.int = TRUE,  ellipse = TRUE,                        
              color = "Anno",palette = "lancet",facet.by = "Bio"        
)+stat_cor(method = "spearman",aes(color = Anno), label.x = 40)

LAT+LON+plot_layout(ncol=1)

ggscatter(df, x = "x", y = "y",
          add = "reg.line",                         
          conf.int = TRUE,  ellipse = TRUE,                        
          color = "Bio",facet.by = "Bio"        
)+stat_cor(method = "spearman",aes(color = Bio), label.x = 40)

###LAT LON corr with value
df_Bio4=read.table("LAT_LON_corrIMP/Bio4.tsv",header = T)
Bio4_Lat=ggscatter(df_Bio4, x = "Latitude", y = "Value",title = "Bio4",
          add = "reg.line",                         
          conf.int = TRUE,  ellipse = TRUE,                        
          color = "Anno",palette = "lancet"      
)+stat_cor(method = "spearman",aes(color = Anno), label.x = 40)
Bio4_Lon=ggscatter(df_Bio4, x = "Longitude", y = "Value",title = "Bio4",
                   add = "reg.line",                         
                   conf.int = TRUE,  ellipse = TRUE,                        
                   color = "Anno",palette = "lancet"      
)+stat_cor(method = "spearman",aes(color = Anno), label.x = 40)



df_Bio5=read.table("LAT_LON_corrIMP/Bio5.tsv",header = T)
Bio5_Lat=ggscatter(df_Bio5, x = "Latitude", y = "Value",title = "Bio5",
                   add = "reg.line",                         
                   conf.int = TRUE,  ellipse = TRUE,                        
                   color = "Anno",palette = "lancet"      
)+stat_cor(method = "spearman",aes(color = Anno), label.x = 40)
Bio5_Lon=ggscatter(df_Bio5, x = "Longitude", y = "Value",title = "Bio5",
                   add = "reg.line",                         
                   conf.int = TRUE,  ellipse = TRUE,                        
                   color = "Anno",palette = "lancet"      
)+stat_cor(method = "spearman",aes(color = Anno), label.x = 40)



df_Bio10=read.table("LAT_LON_corrIMP/Bio10.tsv",header = T)
Bio10_Lat=ggscatter(df_Bio10, x = "Latitude", y = "Value",title = "Bio10",
                    add = "reg.line",                         
                    conf.int = TRUE,  ellipse = TRUE,                        
                    color = "Anno",palette = "lancet"      
)+stat_cor(method = "spearman",aes(color = Anno), label.x = 40)
Bio10_Lon=ggscatter(df_Bio10, x = "Longitude", y = "Value",title = "Bio10",
                    add = "reg.line",                         
                    conf.int = TRUE,  ellipse = TRUE,                        
                    color = "Anno",palette = "lancet"      
)+stat_cor(method = "spearman",aes(color = Anno), label.x = 40)



Bio5_Lat+Bio10_Lat+Bio4_Lat+Bio5_Lon+Bio10_Lon+Bio4_Lon+plot_layout(ncol=3)



# ggscatter(df, x = "Lat", y = "y",
#           add = "reg.line",                         # Add regression line
#           conf.int = TRUE,                          # Add confidence interval
#           color = "Anno", palette = "jco",           # Color by groups "cyl"
#           shape = "Anno"                             # Change point shape by groups "cyl"
# )+stat_cor(aes(color = Anno), label.x = 30) 
# 
# facet.by = "dataset"
# ggscatter(df, x = "value", y = "imp",
#           add = "reg.line",                         
#           conf.int = TRUE,  ellipse = TRUE,                        
#           color = "Anno",palette = "lancet",facet.by = "Type"        
# )+stat_cor(method = "spearman",aes(color = Anno), label.x = 50)  











#heritability
/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/heritability
ls|grep csv|awk '{print "sed \0471d\047 "$1"|awk -F\",\" \047{print \"0\\t\"\$1\"\\t\"\$2}\047 > "$1".gcta"}'|sed 's/csv.gcta/gcta/'|sh
plink2 --vcf ../vcf2gwas/test.vcf.gz -chr-set -8 --make-bed --out file3
gcta64 --bfile file3 --autosome --maf 0.01 --make-grm --out test --thread-num 10
gcta64 --grm test --pheno test.pheno --reml --out test --thread-num 10
ls|grep gcta|awk '{print "gcta64 --grm test --pheno "$1" --reml --out "$1".hert --thread-num 10"}'|sed 's/gcta.hert/herit/'
more *herit.hsq|egrep 'herit|V\(G\)/'|sed 'N;s/\n/\t/'|awk '{print $1"\t"$3"\t"1-$3}'|awk '$2>0.6'












