
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







####gwas####

#Bio10
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/3-lfmm")
gwasResults<-read.delim("bio10.pvalue.gwas.result",header=TRUE,sep="\t")


chr_len <- gwasResults %>%
  group_by(CHR) %>%
  summarise(chr_len=max(BP))

chr_pos <- chr_len  %>%
  mutate(total = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)

Snp_pos <- chr_pos %>%
  left_join(gwasResults, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)


X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
#colors=rep(c("#3B6EFF", "#A4C8FE"), 8 )
colors=c("#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF")

max_p_Snp=-log10(0.05/3479056)

ggplot(Snp_pos, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=CHR))+
  scale_color_manual(values = colors) +
  scale_x_continuous( label = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"), breaks= X_axis$center ) +
  labs(x="",y=paste("-log10(P):BIO10")) +
  scale_y_continuous(limits=c(0,-log10(min(Snp_pos$P))+0.5),expand = c(0, 0) ) +
  geom_hline(yintercept =10.985, color = "grey", size = 2, alpha=0.8, linetype = "twodash")+
  theme_bw() +
  theme(
    plot.title=element_text(size = 15,vjust = -0.2,hjust=0.9,face = 'bold'),
    axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
    axis.text.x=element_text(vjust=1,hjust=0,size=20,colour="black"),
    axis.text.y=element_text(vjust=1,hjust=0.9,size=20,colour="black"),
    axis.title.y = element_text(size =20),
    axis.line.x=element_line(colour="black"),
    axis.line.y=element_line(colour="black"),
    panel.border=element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 11),
    panel.background=element_rect(fill="white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
    plot.margin=unit(c(0.5,0.5,0.5,0.5),"mm"))
ggsave(filename="my.bio10.gwas.3.png",width=24,height = 14)


#Bio4
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/3-lfmm")
gwasResults<-read.delim("bio4.pvalue.gwas.result",header=TRUE,sep="\t")


chr_len <- gwasResults %>%
  group_by(CHR) %>%
  summarise(chr_len=max(BP))

chr_pos <- chr_len  %>%
  mutate(total = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)

Snp_pos <- chr_pos %>%
  left_join(gwasResults, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)


X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
#colors=rep(c("#3B6EFF", "#A4C8FE"), 8 )
colors=c("#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF")


ggplot(Snp_pos, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=CHR))+
  scale_color_manual(values = colors) +
  scale_x_continuous( label = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"), breaks= X_axis$center ) +
  labs(x="",y=paste("-log10(P):BIO4")) +
  scale_y_continuous(limits=c(0,-log10(min(Snp_pos$P))+0.5),expand = c(0, 0) ) +
  geom_hline(yintercept =9.3, color = "grey", size = 2, alpha=0.8, linetype = "twodash")+
  theme_bw() +
  theme(
    plot.title=element_text(size = 15,vjust = -0.2,hjust=0.9,face = 'bold'),
    axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
    axis.text.x=element_text(vjust=1,hjust=0,size=20,colour="black"),
    axis.text.y=element_text(vjust=1,hjust=0.9,size=20,colour="black"),
    axis.title.y = element_text(size =20),
    axis.line.x=element_line(colour="black"),
    axis.line.y=element_line(colour="black"),
    panel.border=element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 11),
    panel.background=element_rect(fill="white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
    plot.margin=unit(c(0.5,0.5,0.5,0.5),"mm"))
ggsave(filename="my.bio4.gwas.3.png",width=24,height = 14)
ggsave(filename="my.bio4.gwas.3.pdf",dpi=300)

#Bio18
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/3-lfmm")
gwasResults<-read.delim("bio18.pvalue.gwas.result",header=TRUE,sep="\t")


chr_len <- gwasResults %>%
  group_by(CHR) %>%
  summarise(chr_len=max(BP))

chr_pos <- chr_len  %>%
  mutate(total = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)

Snp_pos <- chr_pos %>%
  left_join(gwasResults, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)


X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
#colors=rep(c("#3B6EFF", "#A4C8FE"), 8 )
colors=c("#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF")


ggplot(Snp_pos, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=CHR))+
  scale_color_manual(values = colors) +
  scale_x_continuous( label = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"), breaks= X_axis$center ) +
  labs(x="",y=paste("-log10(P):BIO18")) +
  scale_y_continuous(limits=c(0,-log10(min(Snp_pos$P))+0.5),expand = c(0, 0) ) +
  geom_hline(yintercept =7.84, color = "grey", size = 2, alpha=0.8, linetype = "twodash")+
  theme_bw() +
  theme(
    plot.title=element_text(size = 15,vjust = -0.2,hjust=0.9,face = 'bold'),
    axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
    axis.text.x=element_text(vjust=1,hjust=0,size=20,colour="black"),
    axis.text.y=element_text(vjust=1,hjust=0.9,size=20,colour="black"),
    axis.title.y = element_text(size =20),
    axis.line.x=element_line(colour="black"),
    axis.line.y=element_line(colour="black"),
    panel.border=element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 11),
    panel.background=element_rect(fill="white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
    plot.margin=unit(c(0.5,0.5,0.5,0.5),"mm"))
ggsave(filename="my.bio18.gwas.3.png",width=24,height = 14)

#Bio5
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/3-lfmm")
gwasResults<-read.delim("bio5.pvalue.gwas.result",header=TRUE,sep="\t")


chr_len <- gwasResults %>%
  group_by(CHR) %>%
  summarise(chr_len=max(BP))

chr_pos <- chr_len  %>%
  mutate(total = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)

Snp_pos <- chr_pos %>%
  left_join(gwasResults, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)


X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
#colors=rep(c("#3B6EFF", "#A4C8FE"), 8 )
colors=c("#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF")


ggplot(Snp_pos, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=CHR))+
  scale_color_manual(values = colors) +
  scale_x_continuous( label = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"), breaks= X_axis$center ) +
  labs(x="",y=paste("-log10(P):BIO5")) +
  scale_y_continuous(limits=c(0,-log10(min(Snp_pos$P))+0.5),expand = c(0, 0) ) +
  geom_hline(yintercept =9.447, color = "grey", size = 2, alpha=0.8, linetype = "twodash")+
  theme_bw() +
  theme(
    plot.title=element_text(size = 15,vjust = -0.2,hjust=0.9,face = 'bold'),
    axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
    axis.text.x=element_text(vjust=1,hjust=0,size=20,colour="black"),
    axis.text.y=element_text(vjust=1,hjust=0.9,size=20,colour="black"),
    axis.title.y = element_text(size =20),
    axis.line.x=element_line(colour="black"),
    axis.line.y=element_line(colour="black"),
    panel.border=element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 11),
    panel.background=element_rect(fill="white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
    plot.margin=unit(c(0.5,0.5,0.5,0.5),"mm"))
ggsave(filename="my.bio5.gwas.3.png",width=24,height = 14)

#Bio8
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/3-lfmm")
gwasResults<-read.delim("bio8.pvalue.gwas.result",header=TRUE,sep="\t")


chr_len <- gwasResults %>%
  group_by(CHR) %>%
  summarise(chr_len=max(BP))

chr_pos <- chr_len  %>%
  mutate(total = cumsum(chr_len) - chr_len) %>%
  select(-chr_len)

Snp_pos <- chr_pos %>%
  left_join(gwasResults, ., by="CHR") %>%
  arrange(CHR, BP) %>%
  mutate( BPcum = BP + total)


X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
#colors=rep(c("#3B6EFF", "#A4C8FE"), 8 )
colors=c("#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF", "#A4C8FE","#3B6EFF")


ggplot(Snp_pos, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=CHR))+
  scale_color_manual(values = colors) +
  scale_x_continuous( label = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8"), breaks= X_axis$center ) +
  labs(x="",y=paste("-log10(P):BIO8")) +
  scale_y_continuous(limits=c(0,-log10(min(Snp_pos$P))+0.5),expand = c(0, 0) ) +
  geom_hline(yintercept =7.84, color = "grey", size = 2, alpha=0.8, linetype = "twodash")+
  theme_bw() +
  theme(
    plot.title=element_text(size = 15,vjust = -0.2,hjust=0.9,face = 'bold'),
    axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
    axis.text.x=element_text(vjust=1,hjust=0,size=20,colour="black"),
    axis.text.y=element_text(vjust=1,hjust=0.9,size=20,colour="black"),
    axis.title.y = element_text(size =20),
    axis.line.x=element_line(colour="black"),
    axis.line.y=element_line(colour="black"),
    panel.border=element_blank(),
    legend.position = "none",
    axis.title.x = element_text(size = 11),
    panel.background=element_rect(fill="white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
    plot.margin=unit(c(0.5,0.5,0.5,0.5),"mm"))
ggsave(filename="my.bio8.gwas.3.png",width=24,height = 14)



##




#heritability
/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/heritability
ls|grep csv|awk '{print "sed \0471d\047 "$1"|awk -F\",\" \047{print \"0\\t\"\$1\"\\t\"\$2}\047 > "$1".gcta"}'|sed 's/csv.gcta/gcta/'|sh
plink2 --vcf ../vcf2gwas/test.vcf.gz -chr-set -8 --make-bed --out file3
gcta64 --bfile file3 --autosome --maf 0.01 --make-grm --out test --thread-num 10
gcta64 --grm test --pheno test.pheno --reml --out test --thread-num 10
ls|grep gcta|awk '{print "gcta64 --grm test --pheno "$1" --reml --out "$1".hert --thread-num 10"}'|sed 's/gcta.hert/herit/'
more *herit.hsq|egrep 'herit|V\(G\)/'|sed 'N;s/\n/\t/'|awk '{print $1"\t"$3"\t"1-$3}'|awk '$2>0.6'












