######Figure 1A geographic map####
setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/2-non-clones-gvcf/3-pi-fst/RDA/miss09/new-lfmm/pvalue/gwas-pvalue/GF")

data <- read.delim("map.txt", sep="\t", header=T)
world_map <- map_data("world")


data$Longitude <- as.numeric(data$Longitude)


# make a map
#ggplot() +geom_polygon(data = world_map, aes(x = long, y =lat, group = group), fill="grey90") +geom_point(data = data, aes(x = Longitude, y = Latitude,color="blue"), size=0.8) +scale_color_aaas()+theme_void() +ylim(-55,85) +theme(legend.position="none")

Figure_1_geo_map=ggplot() +geom_polygon(data = world_map, aes(x = long, y =lat, group = group), fill="grey90") +geom_point(data = data, aes(x = Longitude, y = Latitude,color="blue"), size=0.8) +scale_color_aaas()+theme_bw()+ylim(-55,85)+theme(legend.position="none")+xlab("Longitude")+ylab("Latitude")
zoom1=Figure_1_geo_map +
  coord_cartesian(xlim = c(-10, 18), ylim = c(45, 57))
zoom2=Figure_1_geo_map +
  coord_cartesian(xlim = c(130, 145), ylim = c(30, 50))



######Figure 1B Genome distance#####
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-7/rudealis/All-data/new/gvcf-209/1-GD-remove-clones")
vcftools --vcf autosome.filter.snp.vcf --plink --out test
plink --file test  --distance-matrix

cp plink.mdist gd.tsv
sed -i 's/ /\t/g' gd.tsv
less -S plink.mdist.id|cut -f1 > id
paste -d"\t" id gd.tsv > GenomeDistance.209.tsv


a=read.table("GenomeDistance.209.tsv",header = T,row.names=1)
a[lower.tri(a)] <- NA
a$y <- rownames(a)
da <- melt(data = a) %>% na.omit()

da$variable <- factor(da$variable,levels = unique(da$variable))
da$y <- factor(da$y,levels = unique(da$y))
head(da,3)

write.table(da,"./autosome.pass.mat.simplify",sep="\t")

less -S autosome.pass.mat.simplify |awk '$NF>0'|awk '{print $NF}' > autosome.pass.mat.simplify.value

a=read.table("autosome.pass.mat.simplify.value",header = F)
ggplot(a,aes(x=V1))+geom_histogram(bins=9500,color = "grey",fill='transparent')+ylim(0,50)+geom_density(color="#03B0AB",alpha=0.2,cex=1,fill='transparent')+theme_classic()+xlab("Genome distance")+ylab("Frequency")+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=12),axis.text.x=element_text(size=12))+geom_vline(aes(xintercept=0.02), colour="#03B0AB", linetype="dotted",size=0.5)

ggplot(a,aes(x=V1))+geom_histogram(bins=9500,aes(y=..ncount..),binwidth = 0.001)+ylim(0,50)+geom_density(color="#03B0AB",alpha=0.2,cex=1,fill='transparent')+theme_classic()+xlab("Genome distance")+ylab("Frequency")+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=12),axis.text.x=element_text(size=12))+geom_vline(aes(xintercept=0.02), colour="#03B0AB", linetype="dotted",size=0.5)

Figure_S_all209_GD=ggplot(a, aes(x = V1)) +
  geom_histogram(bins = 1000, aes(y = ..density.. / 100), binwidth = 0.0001, color = "grey",fill='transparent') +
  geom_density(aes(y = ..density.. / 100), color = "#03B0AB") +geom_vline(aes(xintercept=0.01), colour="blue", linetype="dotted",size=0.7)+
  xlab("Hamming distance")+ylab("Frequency")+
  theme_classic()



/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/1-GD-remove-clones
setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/1-GD-remove-clones")

less -S Table.S.GD |grep Bonn|awk '$1~/Bonn/ && $2~/Bonn/'|cut -f3 > Bonn.GD.value
a=read.table("Bonn.GD.value",header = F)




ggplot(a,aes(x=V1))+geom_histogram(aes(y=..ncount..),binwidth = 0.001)+theme_classic()



#Test Bonn


#One FigureS for all

less -S list |awk '{print $1"\t"$9"-"$10}'|grep -v '\-N' > 12.pops.name

perl /groups/dolan/user/shuangyang.wu/script/merge_files.pl -r Table.S.GD -n 1 -i 12.pops.name -c 1 -o out.field1
perl /groups/dolan/user/shuangyang.wu/script/merge_files.pl -r out.field1 -n 2 -i 12.pops.name -c 1 -o out.field2

less -S out.field2 |awk '{print $5"\t"$7"\t"$3}'|grep -v '^-'|awk '!($2~/^-/)'|awk '$1==$2'|sort -k1,1|sed 's/-Y//g' > 12.pops.GD.value


GD_12=read.table("12.pops.GD.value",header = F)

ggplot(GD_12,aes(x=V3))+geom_histogram(aes(y=..ncount..),binwidth = 0.001)+geom_vline(aes(xintercept=0.01), colour="blue", linetype="dotted",size=0.7)+xlab("Hamming distance")+ylab("Frequency")+theme(axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.y=element_text(size=12),axis.text.x=element_text(size=12))+theme_bw()+facet_wrap(vars(V1), nrow = 4)

#####FIGURE 1 combine Figure_1_geo_map+Figure_S_all209_GD+plot_annotation(tag_levels = 'A')###
######Figure 2A PCA#####
setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/1-PCA")
table=read.table("78.info.pca",header = T,sep = "\t")
country_colours <-
  c("Austria" = "#00A087",
    "Canada" = "#3C5488",
    "Croatia" = "#4DBBD5",
    "Czech" = "#9DAAC4",
    "Denmark" = "#E64B35",
    "UK" = "#8c2d04",
    "France"= "#762a83",
    "Germany"="#ec7014",
    "Hungary"="#fee090",
    "Ireland"="#c51b7d",
    "Japan"="black",
    "Switzerland"="#f4a582",
    "US" = "#5e4fa2"
  )
all_pca=ggplot(table,aes(x=PC1,y=PC2,color=Country,label =RecordingName ),alpha=1)+geom_point(size=3, alpha=1)+theme_bw()+scale_colour_manual(values = country_colours)+xlab("PCA1 8.0%")+ylab("PCA2 3.0%")
#ggplot(table,aes(x=PC1,y=PC2,color=Note,label =RecordingName ),alpha=1)+geom_point(size=3, alpha=1)+theme_bw()+scale_colour_manual(values = country_colours)+xlab("PCA1 8.0%")+ylab("PCA2 3.0%")+geom_text_repel(label=table[,11],max.overlaps = Inf)

setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/1-PCA/onlyEurope")
table=read.table("54.info.pca",header = T,sep = "\t")
europe_pca=ggplot(table,aes(x=PC1,y=PC2,color=Country,label =RecordingName ),alpha=1)+geom_point(size=3, alpha=1)+theme_bw()+scale_colour_manual(values = country_colours)+xlab("PCA1 4.6%")+ylab("PCA2 3.5%")+theme(legend.position = "NA")
all_pca|europe_pca


######Figure 2B Admixture####
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-7/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/2-Admixture")
a=read.table("K.dis",header = F,sep="\t")
ggplot(a,aes(x =V1,y = V2,group = 1))+geom_line(linetype= "dashed",color="navyblue",size=1)+geom_point(size=3)+theme_classic()+xlab("K value")+ylab("CV error")

setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-7/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/2-Admixture/Q")

library(pophelper)
library(ggplot2)
require(gridExtra)

grep 'CHROM' *.vcf > id
perl /groups/dolan/user/shuangyang.wu/script/merge_files.pl -r id -n 1 -i 78.info -c 1 -o out
less -S out |awk '{print $10"\t"$1}' > new.id
###only Q
sfiles <- list.files(path="/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/2-Admixture/Q", full.names=T)
labels <- read.csv("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/2-Admixture/new.1.id", header=F, stringsAsFactors=F, sep = "\t")
#labels <- read.csv("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/2-Admixture/new.id", header=F, stringsAsFactors=F, sep = "\t")
slist <- readQ(files=sfiles)
tabulateQ(qlist=readQ(sfiles))
summariseQ(tabulateQ(qlist=readQ(sfiles)))
rownames(slist[[1]]) <- labels$V2
rownames(slist[[2]]) <- labels$V2
rownames(slist[[3]]) <- labels$V2
rownames(slist[[4]]) <- labels$V2
rownames(slist[[5]]) <- labels$V2
rownames(slist[[6]]) <- labels$V2
rownames(slist[[7]]) <- labels$V2
rownames(slist[[8]]) <- labels$V2
rownames(slist[[9]]) <- labels$V2

p1<-plotQ(slist,returnplot=T,exportplot=F,basesize=11,showindlab=T,useindlab=TRUE,showyaxis=T,showticks=T,sortind="all")

grid.arrange(p1$plot[[1]],p1$plot[[2]],p1$plot[[3]],p1$plot[[4]],p1$plot[[5]],p1$plot[[6]],p1$plot[[7]],p1$plot[[8]],p1$plot[[9]],ncol=1)

###test without id###
p1<-plotQ(slist,returnplot=T,exportplot=F,basesize=11,showindlab=F,useindlab=TRUE,showyaxis=T,showticks=T,sortind="all")

grid.arrange(p1$plot[[1]],p1$plot[[2]],p1$plot[[3]],p1$plot[[4]],p1$plot[[5]],p1$plot[[6]],p1$plot[[7]],p1$plot[[8]],p1$plot[[9]],ncol=1)

###test withgroupID

p1<-plotQ(slist,returnplot=T,exportplot=F,basesize=11,grplab=data.frame(lab1=labels$V1,stringsAsFactors=F),showindlab=TRUE,useindlab=TRUE,showyaxis=T,showticks=T)
###test withgroupID test no id###
p1<-plotQ(slist,returnplot=T,exportplot=F,basesize=11,grplab=data.frame(lab1=labels$V1,stringsAsFactors=F),showindlab=F,useindlab=TRUE,showyaxis=T,showticks=T)

grid.arrange(p1$plot[[2]],p1$plot[[3]],p1$plot[[4]],p1$plot[[5]],p1$plot[[6]],p1$plot[[7]],p1$plot[[8]],p1$plot[[9]],ncol=1)

###GROUPID final
less -S new.id|sed 's/^/E-/'|sed 's/E-JPN/J-JPN/' |sed 's/E-AUS/E-AUT/'|sed 's/E-CRO/E-HRV/'|sed 's/E-DEN/E-DNK/'|sed 's/E-GER/E-DEU/'|sed 's/E-IRE/E-IRL/'|sed 's/E-SWZ/E-CHE/'|sed 's/E-US/E-USA/'|sed 's/E-UK/E-GBR/' > new.1.id
p1<-plotQ(slist,returnplot=T,exportplot=F,basesize=11,grplab=data.frame(lab1=labels$V1,stringsAsFactors=F),showindlab=F,useindlab=TRUE,showyaxis=T,showticks=T,ordergrp = TRUE)
grid.arrange(p1$plot[[2]],p1$plot[[3]],p1$plot[[4]],p1$plot[[5]],p1$plot[[6]],p1$plot[[7]],p1$plot[[8]],p1$plot[[9]],ncol=1)


######Figure 3B fst heatmap####
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-7/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst/Fst-vcftools")
df<-read.table("fst.group.pair",head=T)
ggplot(df,aes(x=region1, y=region2))+
  geom_tile(aes(fill=fst), color = 'white',alpha = 0.6) +geom_text(aes(label=round(fst,2)))+scale_fill_gradient(low = "white",high = "#2c7fb8")+theme_grey(base_size = 10)+
  labs(x = NULL,y = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 8,angle=0,hjust = 1),
        axis.text = element_text(color='black'),
        panel.background = element_blank())

#new fst matrix

/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/new-7-fst-IBD-IBE

pixy --stats pi fst dxy --vcf test_filtered1.vcf.gz --populations 7.sample.info --window_size 20000  --n_cores 10 --output_prefix filtered_7_pixy_20K  --bypass_invariant_check yes

#cat list | awk '{print "less -S all_pixy_20K_fst.txt|grep "$1"|grep "$2" |cut -f6|grep -v NA|awk \047{if(\$1<0){print \"0\"}else{print \$0}}\047|datamash mean 1 >"$1"-"$2".fst"}' | sh &
  
less -S filtered_7_pixy_20K_fst.txt|sed 's/\t/_/'|cut -f1,2,5|sed '1d'|egrep -v 'chrU|chrV'|sort -k1,1 -k2,2|grep -v NA|awk '{if($3<0){print $1"\t"$2"\t0"}else{print $0}}'|datamash -g 1 mean 3 >fst.pair.tsv

less -S filtered_7_pixy_20K_dxy.txt|sed 's/\t/_/'|cut -f1,2,5|sed '1d'|egrep -v 'chrU|chrV'|sort -k1,1 -k2,2|grep -v NA|awk '{if($3<0){print $1"\t"$2"\t0"}else{print $0}}'|datamash -g 1 mean 3|sed 's/_/\t/2' > dxy.pair.tsv

setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/new-7-fst-IBD-IBE")
df<-read.table("fst.pair.tsv",head=T)
Fst_7=ggplot(df,aes(x=region1, y=region2))+
  geom_tile(aes(fill=fst), color = 'white',alpha = 0.6)+geom_text(aes(label=round(fst,5)))+scale_fill_gradient(low = "white",high = "#0897B4")+theme_grey(base_size = 10)+
  labs(x = NULL,y = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 8,angle=0,hjust = 1),
        axis.text = element_text(color='black'),
        panel.background = element_blank())

Use the information from last time, update the plot data with new code

/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/new-7-fst-IBD-IBE
setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/new-7-fst-IBD-IBE")
plot_data=read.csv("plot.table.final",head=T)
##manully remove 0 line
library(ggpubr)
p_dis<-ggscatter(plot_data, x = "geo_dit", y = "fst",
                 size = 2.5,
                 add = "reg.line",  # 添加回归线
                 add.params = list(color="#9c9c9c", fill="#d6d6d6", size = 1),  # 自定义回归线的颜色
                 conf.int = TRUE  # 添加置信区间
)+stat_cor(method = "spearman",label.sep = "\n") +
  labs(x = "Geographical Distance (m)",y = expression(italic(F)[italic(ST)]/(1-italic(F)[italic(ST)])),size = 3)

p_env<-ggscatter(plot_data, x = "env_dist", y = "fst",
                 size = 2.5,
                 add = "reg.line",  # 添加回归线
                 add.params = list(color="#9c9c9c", fill="#d6d6d6", size = 1),  # 自定义回归线的颜色
                 conf.int = TRUE  # 添加置信区间
)+stat_cor(method = "spearman", label.sep = "\n") +
  labs(x = "Environment distance",y = expression(italic(F)[italic(ST)]/(1-italic(F)[italic(ST)])),size = 3)

p_dis|p_env


pi_data <- read.table("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/new-7-fst-IBD-IBE/filtered_7_pixy_20K_pi.txt", header=T)
Pi_7=ggplot(pi_data, aes(pop,avg_pi,color=pop),guide="none")+geom_jitter()+geom_boxplot(fill=NA, col="black")+scale_color_aaas()+labs(x = "Population" ,y = "Nucleotide Diversity (Pi)") +theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position = "none")


df1<-read.table("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/new-7-fst-IBD-IBE/dxy.pair.tsv",head=T)
dxy_7=ggplot(df1,aes(x=region1, y=region2))+
  geom_tile(aes(fill=dxy), color = 'white',alpha = 0.6) +geom_text(aes(label=round(dxy,5)))+scale_fill_gradient(low = "white",high = "#7BCDC3")+theme_grey(base_size = 10)+
  labs(x = NULL,y = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 8,angle=0,hjust = 1),
        axis.text = element_text(color='black'),
        panel.background = element_blank())

FIGURE 3 combine Pi_7+dxy_7+Fst_7+p_dis+p_env+plot_annotation(tag_levels = 'A')




#bootstrap for equal number of ind from 10 to 21

for n in {10..22}; do
# Inner loop from 1 to 100
seq 1 100 | awk -v n="$n" '{
        print "less -S sample.info | grep Europe | shuf -n "n" > "n".n" n ".sample_" $1 ".europe.japan.sample.info";
        print "less -S sample.info | grep Japan | shuf -n "n" >> "n".n" n ".sample_" $1 ".europe.japan.sample.info"
    }' | bash
done

/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/bootstrap/n10

sbatch -c 10 -J fst -a 1-1200 array.sh

ls|grep pixy_20K_pi.txt|awk '{print "cat "$1"|egrep -v \047chrU|chrV\047|sed 1d|sort -k1,1 -k2,2 -k3,3n|datamash -g 1 median 5 |awk \047{print \""$1"\\t\"\$0}\047 > "$1".result"}'|sed 's/_pixy_20K_pi.txt//2'|sed 's/_pixy_20K_pi.txt//2'|sh &
  
setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/bootstrap/n10")
cat *.result|less -S|sed -r 's/\.\S+//' > bootstrap.pi.result

#Figure BOOTSTRAP
BS_pi=read.table("bootstrap.pi.result",header = F)

ggplot(BS_pi, aes(V2,V3,color=V1),guide="none")+geom_jitter()+geom_boxplot(fill=NA, col="black")+scale_color_viridis()+labs(x = "Population" ,y = "Nucleotide Diversity (Pi)") +theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position = "none")+facet_wrap(~V1)

  



######Figure 3C IBD IBE####



p_dis<-ggscatter(plot_data, x = "geo_dit", y = "fst",
          size = 2.5,
          add = "reg.line",  # 添加回归线
          add.params = list(color="#9c9c9c", fill="#d6d6d6", size = 1),  # 自定义回归线的颜色
          conf.int = TRUE  # 添加置信区间
)+stat_cor(method = "spearman",label.x = 10000, label.y = 0.063, label.sep = "\n") +
  labs(x = "Geographical Distance (m)",y = expression(italic(F)[italic(ST)]/(1-italic(F)[italic(ST)])),size = 3)

p_env<-ggscatter(plot_data, x = "env_dist", y = "fst",
                 size = 2.5,
                 add = "reg.line",  # 添加回归线
                 add.params = list(color="#9c9c9c", fill="#d6d6d6", size = 1),  # 自定义回归线的颜色
                 conf.int = TRUE  # 添加置信区间
)+stat_cor(method = "spearman",label.x = 1, label.y = 0.063, label.sep = "\n") +
  labs(x = "Environment distance",y = expression(italic(F)[italic(ST)]/(1-italic(F)[italic(ST)])),size = 3)

p_dis|p_env

#####Figure 3 A 7 pi######
pi_data <- read.table("7.1.pi", header=T)
ggplot(pi_data, aes(pop,avg_pi,color=pop),guide="none")+geom_jitter()+geom_boxplot(fill=NA, col="black")+scale_color_aaas()+labs(x = "Population" ,y = "Nucleotide Diversity (Pi)") +theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1))


df1<-read.table("7.dxy.value",head=T)
ggplot(df1,aes(x=region1, y=region2))+
  geom_tile(aes(fill=dxy), color = 'white',alpha = 0.6) +geom_text(aes(label=round(dxy,4)))+scale_fill_gradient(low = "white",high = "#2c7fb8")+theme_grey(base_size = 10)+
  labs(x = NULL,y = NULL) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + 
  theme(legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_text(size = 8,angle=0,hjust = 1),
        axis.text = element_text(color='black'),
        panel.background = element_blank())


#####Figure S3 and 4A pi main figure####
setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan")
pi_data <- read.table("default_pixy_noUV_20K_pi.txt", header=T)
pi_data_euro <- pi_data %>%
  dplyr::filter(pop=="Europe")
pi_data_asia <- pi_data %>%
  dplyr::filter(pop=="Japan")
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

pi_data %>%
  group_by(chromosome) %>%
  summarise(max = max(position, na.rm = TRUE))

pi_data %>%
  group_by(pop) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))

pi_data %>% ggplot()+aes(position*20000, avg_pi, color = factor(pop, levels = c("Europe", "Japan"), labels = c("Europe", "Japan")),fill = factor(pop, levels = c("Europe", "Japan"), labels = c("Europe", "Japan")))+geom_point(alpha = 0.5, shape = 21, color = "black") +facet_grid(.~chromosome, scales = "free")+geom_smooth(span = 0.4, se = F, method = "loess")+scale_color_manual(values = c("#FDD876", "#87c9c3"))+scale_fill_manual(values = c("#FDD876", "#87c9c3"))+theme_bw()+
  theme(legend.box = NULL,
        strip.background = element_blank(),
        legend.background = element_rect(colour = NA),
        legend.key = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "Pi", color = "Population", fill = "Population")->Pi_Figure  

###For each pi accross the chromosome
ggplot(pi_data, aes(chromosome,avg_pi,color=pop),guide="none")+geom_boxplot()+scale_color_aaas()+facet_grid(.~chromosome,scales = "free")+stat_compare_means()+theme_bw()+labs(y = "Pi")

###For all pi in the genome
ggplot(pi_data, aes(pop,avg_pi,color=pop),guide="none")+geom_jitter()+geom_boxplot(fill=NA, col="black")+labs(y = "Pi") +theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1))+scale_color_aaas()+stat_compare_means() 

#####Figure 4 fst####
library(tidyverse)
library(ggridges)

# load data
dxy_data <- read.table("default_pixy_noUV_20K_dxy.txt", header=T)
fst_data <- read.table("default_pixy_noUV_20K_fst_over0.txt", header=T)

# add some columns
dxy_data$data_type <- "Dxy"
dxy_data <- mutate(dxy_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))

fst_data$data_type <- "Fst"
fst_data <- mutate(fst_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))

# subset and merge dataframes
dxy_data_sub <- dxy_data %>% dplyr::select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

fst_data_sub <- fst_data %>% dplyr::select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
colnames(fst_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")

# cheeky fix to get matching rows from both datasets
tmp_dxy <- semi_join(dxy_data_sub, fst_data_sub, by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))
tmp_fst <- semi_join(fst_data_sub, dxy_data_sub,  by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))
data <- full_join(tmp_dxy, tmp_fst)
data <- data %>%
  group_by(comparison, data_type) %>%
  mutate(position = 1:n())
data %>%
  group_by(chromosome) %>%
  summarise(max = max(position, na.rm = TRUE))

data %>%group_by(comparison,data_type) %>%summarise(median = median(value, na.rm = TRUE))
view(data %>%group_by(comparison,data_type) %>%summarise(median = median(value, na.rm = TRUE)))

tmp_fst <- data %>% dplyr::filter(data_type=="Fst")
head(tmp_fst)


tmp_fst%>%ggplot() +
  aes(x = position*20000, y = value)+
  geom_point(alpha = 0.5, color = "black", fill = "#999999", shape = 21)+
  geom_smooth(span = 0.4, se = F, method = "loess")+
  facet_grid(.~chromosome, scales = "free")+
  theme_bw()+
  theme(legend.box = NULL,
        strip.background = element_blank(),
        legend.background = element_rect(colour = NA),
        legend.key = element_blank(),
        strip.text.x = element_blank(),
        axis.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "Fst")->FST_Figure



#####Figure 3 S ABCDE####

Pi_chr<-ggplot(pi_data, aes(chromosome,avg_pi,color=pop),guide="none")+geom_boxplot()+scale_color_aaas()+facet_grid(.~chromosome,scales = "free")+theme_bw()+labs(y = "Pi")


Pi_all<-ggplot(pi_data, aes(pop,avg_pi,color=pop),guide="none")+geom_jitter()+geom_boxplot(fill=NA, col="black")+labs(x = "Population" , y = "Nucleotide Diversity (Pi)") +theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1))+scale_color_aaas()+stat_compare_means() 

Fst_all=ggplot(tmp_fst, aes(value,comparison, fill=comparison), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +theme_bw() + theme(legend.position = "none", axis.text.y=element_blank(),axis.text = element_text(size = 5)) +xlim(0,1) +scale_fill_brewer(palette = "Set2") +labs(x="Fst value of autosomes", y="Density")

Fst_chr=ggplot(tmp_fst, aes(value,chromosome, fill=chromosome), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +theme_bw() + theme(legend.position = "none", axis.text.y=element_blank(),axis.text = element_text(size = 5)) +xlim(0,1) +geom_vline(aes(xintercept=0.137), colour="navyblue", linetype="dashed")+scale_fill_brewer(palette = "Set1") +labs(x="Fst value of each autosome", y="Density")

tmp_fst %>%group_by(chromosome) %>%summarise(median = median(value, na.rm = TRUE))
chromosome median
<chr>       <dbl>
1 chr1      0.112
2 chr2        0.129
3 chr3        0.128
4 chr4        0.120
5 chr5        0.194
6 chr6        0.119
7 chr7        0.156
8 chr8        0.168


Pi_all+Pi_chr+Fst_all+Fst_chr+plot_layout(ncol=2)

######Figure 4 B centromere######
region_fst<-read.table("refine.fst.txt",header = T)
library(ggpubr)
ggplot(region_fst, aes(region,avg_wc_fst,color=region),guide="none")+geom_jitter()+geom_boxplot(fill=NA, col="black")+labs(x = "Region" , y = "Fst") +theme_bw() +theme(axis.text.x = element_text(angle = 45, hjust=1))+facet_grid(chromosome~.)


library('ggthemr')
ggthemr("dust")
ggthemr("flat")
compaired <- list(c("centromere", "long"), 
                  c("centromere","short"), 
                  c("short","long"))
ggplot(region_fst, aes(region,avg_wc_fst,fill=region),guide="none")+geom_violin(draw_quantiles = c(0.5))+labs(y = "Fst") +geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test)+scale_fill_brewer(palette = "Set1")+facet_wrap(chromosome~.,ncol = 4,nrow = 2)

####Figure 5 Sex#####
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/UV-specific")
pi_data <- read.table("chrVU.20K.pi", header=T)

pi_data_euro <- pi_data %>%
  dplyr::filter(pop=="Europe")
pi_data_asia <- pi_data %>%
  dplyr::filter(pop=="Japan")
pi_data <- pi_data %>%
  group_by(pop) %>%
  mutate(position = 1:n())

pi_data %>%
  group_by(chromosome) %>%
  summarise(max = max(position, na.rm = TRUE))

pi_data %>%
  group_by(pop) %>%
  summarise(median = median(avg_pi, na.rm = TRUE))

pi_data %>% ggplot()+aes(position*20000, avg_pi, color = factor(pop, levels = c("Europe", "Japan"), labels = c("Europe", "Japan")),fill = factor(pop, levels = c("Europe", "Japan"), labels = c("Europe", "Japan")))+geom_point(alpha = 0.5, shape = 21, color = "black") +facet_grid(.~chromosome, scales = "free")+geom_smooth(span = 0.4, se = F, method = "loess")+scale_color_manual(values = c("#FDD876", "#87c9c3"))+scale_fill_manual(values = c("#FDD876", "#87c9c3"))+theme_bw()+
  theme(legend.box = NULL,
        strip.background = element_blank(),
        legend.background = element_rect(colour = NA),
        legend.key = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "Pi", color = "Population", fill = "Population")->Pi_Figure_UV

###For each pi accross the autosome
Pi_chr_UV<-ggplot(pi_data, aes(chromosome,avg_pi,color=pop),guide="none")+geom_boxplot()+scale_color_manual(values = c("#FDD876", "#87c9c3"))+facet_grid(.~chromosome,scales = "free")+theme_bw()+labs(y = "Pi")
Pi_Figure_UV|Pi_chr_UV

setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/UV-specific")
dxy_data <- read.table("chrUV.20K.dxy", header=T)
fst_data <- read.table("chrVU.20K.fst", header=T)
dxy_data$data_type <- "Dxy"
dxy_data <- mutate(dxy_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))
fst_data$data_type <- "Fst"
fst_data <- mutate(fst_data,
                   comparison = paste(pop1, pop2, sep = '_v_'))
# subset and merge dataframes
dxy_data_sub <- dxy_data %>% dplyr::select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_dxy)
colnames(dxy_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")
fst_data_sub <- fst_data %>% dplyr::select(comparison,data_type,chromosome,window_pos_1,window_pos_2,avg_wc_fst)
colnames(fst_data_sub) <- c("comparison","data_type","chromosome","window_pos_1","window_pos_2","value")
# cheeky fix to get matching rows from both datasets
tmp_dxy <- semi_join(dxy_data_sub, fst_data_sub, by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))
tmp_fst <- semi_join(fst_data_sub, dxy_data_sub,  by=c("comparison", "chromosome", "window_pos_1", "window_pos_2"))
data <- full_join(tmp_dxy, tmp_fst)
data <- data %>%
  group_by(comparison, data_type) %>%
  mutate(position = 1:n())
data %>%
  group_by(chromosome) %>%
  summarise(max = max(position, na.rm = TRUE))
data %>%group_by(comparison,data_type) %>%summarise(median = median(value, na.rm = TRUE))
view(data %>%group_by(comparison,data_type) %>%summarise(median = median(value, na.rm = TRUE)))
data %>%group_by(chromosome,data_type) %>%summarise(median = median(value, na.rm = TRUE))






Fst_chr_UV=ggplot(tmp_fst, aes(value,chromosome, fill=chromosome), guide="none") +
  geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +theme_bw() + theme(legend.position = "none", axis.text.y=element_blank(),axis.text = element_text(size = 5)) +xlim(0,0.5)+scale_fill_manual(values = c("#FDD876", "#87c9c3")) +labs(x="Fst value of each autosome", y="Density")
Fst_chr_UV
Pi_Figure_UV+Pi_chr_UV+FST_UV+Fst_chr_UV+plot_layout(ncol=2)
###Each chr median

data %>%group_by(chromosome,data_type) %>%summarise(median = median(value, na.rm = TRUE))

ggplot(pi_data, aes(chromosome,avg_pi,color=pop),guide="none")+geom_boxplot()+stat_compare_means()+scale_color_manual(values = c("#FDD876", "#87c9c3"))+facet_grid(.~chromosome,scales = "free")+theme_bw()+labs(y = "Pi")


####try tjmd#####
tjmd_data<-read.table("all.V.tjmd", header=T)
tjmd_data %>% ggplot()+aes(end, tjmd, color = factor(pop, levels = c("Europe", "Japan"), labels = c("Europe", "Japan")),fill = factor(pop, levels = c("Europe", "Japan"), labels = c("Europe", "Japan")))+geom_point(alpha = 0.5, shape = 21, color = "black") +facet_grid(.~chr, scales = "free")+geom_smooth(span = 0.4, se = F, method = "loess")+scale_color_manual(values = c("#FDD876", "#87c9c3"))+scale_fill_manual(values = c("#FDD876", "#87c9c3"))+theme_bw()+
  theme(legend.box = NULL,
        strip.background = element_blank(),
        legend.background = element_rect(colour = NA),
        legend.key = element_blank(),
        legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_blank())+
  labs(x = "Genomic Position (Mb)", y = "TajimaD", color = "Population", fill = "Population")

#####LD and SMC refer to all code#####
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/5-LD")
read.table("test.all.final.Japanese")->EJapanese;
plot(EJapanese[,1]/1000,EJapanese[,2],type="l",col="#00AFBB",main="LD decay",xlab="Distance(Kb)",xlim=c(0,300),ylim=c(0,0.521414665371257),ylab=expression(r^{2}),bty="n",lwd=2)
read.table("test.all.final.European")->EEuropean;
lines(EEuropean[,1]/1000,EEuropean[,2],col="#E7B800",lwd=2)
read.table("test.all.final.Japanese")->EJapanese;
plot(EJapanese[,1]/1000,EJapanese[,2],type="l",col="#00AFBB",main="LD decay",xlab="Distance(Kb)",xlim=c(0,300),ylim=c(0,0.521414665371257),ylab=expression(r^{2}),bty="n",lwd=2)
read.table("test.all.final.European")->EEuropean;
lines(EEuropean[,1]/1000,EEuropean[,2],col="#E7B800",lwd=2)
#Japan
abline(h=0.3262,v=4.3,lwd=1.5,col="grey",lty=3)
abline(h=0.1694,v=1.6,lwd=1.5,col="grey",lty=3)
abline(h=0.3262,v=4.3,lwd=1.5,col="#00AFBB",lty=3)
#European
abline(h=0.1694,v=1.6,lwd=1.5,col="#E7B800",lty=3)
legend("topright",c("Japanese","European"),col=c("#00AFBB","#E7B800"),cex=1,lty=c(1,1),bty="n",lwd=2);
#####GWAS#####

#1 climate factor
snp id with pheno id
##climate /scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/worldclim
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/worldclim")

library(rgdal)
library(sp)
library(raster)

longla<-read.csv("update.new.csv") #读取经纬度文件
coordinates(longla)=c("Lon","Lat") #设置你的经纬度
#worldclim1<-raster("wc2.1_2.5m_bio_4.tif")  #加载栅格数据
#data1<-extract(x= worldclim1,y=longla) #按照坐标提取值

n=c(1:19)
result <- data.frame(row.names = c(1:78)) #1:XX数字取决于你的坐标数据
for (i in n) {
  file=paste0("wc2.1_2.5m_bio_",i,".tif")
  print(file) #把变量打印到屏幕中
  worldclim1<-raster(file)  #加载栅格数据
  data1<-extract(x= worldclim1,y=longla) #按照坐标提取值
  data_all<-data.frame(longla,data1) #合并
  write.csv(data_all,paste0("78-bio",i,".csv"))#导出气候数据
  bio <- read.csv(paste0("78-bio",i,".csv"), as.is = T ) #导入每个点的气候数据
  bio <- bio$data1 #提取指定列data1
  result <- cbind(result,bio) #合并数据（向量）
  colnames(result)[i] <- paste0("bio",i)}
write.csv(result, "78-bio1-19.csv") #转化为csv格式导出

paste -d"," new.csv 78-bio1-19.csv|cut -d"," -f1,5-23|sed 's/\"//g' > pheno.csv

less -S pheno.csv |sed 's/,/\t/g' > pheno1
perl /groups/dolan/user/shuangyang.wu/script/merge_files.pl -r pheno1 -n 1 -i all.info -c 1 -o out
less -S out |cut -f1-20,27-29|sed 's/,/\./g'|awk '{OFS="\t"}{print $1,$23,$21,$22,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}'> 78.geo.bio1-19.tsv

grep -f id 78.geo.bio1-19.tsv|awk '!a[$2]++'|cut -f2-23 > 22.pops.geo.bio1-19.tsv


# correalation

setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/1-pheno")
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(psych)
library(corrplot)


raw_data<-read.csv("78-bio1-19.csv",header=TRUE)
raw_data=raw_data[,2:20]
head(raw_data)

#raw_data<-raw_data[!duplicated(raw_data,fromLast=TRUE),]
col3 <- colorRampPalette(c("#0072B5FF", "white", "#BC3C29FF")) 
corrmatrix = cor(raw_data, method = "spearman")
rownames(corrmatrix) = as.character(names(raw_data))
colnames(corrmatrix) = as.character(names(raw_data))

testRes = cor.mtest(corrmatrix, conf.level = 0.95)

corrplot(corrmatrix, p.mat = testRes$p, method = 'circle', diag = TRUE, type = 'upper', sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, insig = 'label_sig', pch.col = 'grey20',col=col3(20))



#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --partition=c
#SBATCH --qos="long"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=10G
#SBATCH --time=5-00:00:00
#SBATCH --output=my.stdout.merge1

module load vcftools/0.1.16-foss-2018b-perl-5.28.0
module load bcftools/1.9-foss-2018b
module load tabix/0.2.6-gcccore-7.3.0

cp merged.nuclear_variants.miss09.final.snps.vcf test.vcf|bgzip|tabix
vcf2gwas -v test.vcf.gz -pf pheno.csv -o test_all -ap -lmm

#lfmm
/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/3-lfmm

less -S ../vcf2gwas/test.vcf.gz |egrep --color=auto 'CHROM|^chr' | cut -f1-2,10-87 | less -S | sed 's/\t/_/' | sed -r 's/\t1:(\S+)/\t1/g' | sed -r 's/\t0:(\S+)/\t0/g' | sed -r 's/\t\.(\S+)/\t9/g' > 78.miss09.geno &
  
less -S 78.miss09.geno|sed '1d'|cut -f2-79 > 78.miss09.geno1.1

library(data.table)
Y<-fread('78.miss09.geno1.1',header = FALSE)
b<-t(Y)
write.table(b,"78.miss09.geno1.format",sep="\t")

less -S 78.miss09.geno1.format | sed '1d' | awk '{OFS="\t"}{$1="";print $0}' | sed -r 's/^\t//' > 78.miss09.geno1.format.1

mv 78.miss09.geno1.format.1 78.miss09.geno1.lfmm

#id correspond to sample id#

Rscript all.lfmm.1.r

ls|grep pvalue|awk '{print "paste -d\"\\t\" snp.id "$1" | sed -r \047s/\\\".+,//\047 | sed \047s/\\\"V1\\\"/P/\047 > "$1".gwas.result"}'

top 0.00005
ls --color=auto | grep --color=auto gwas.result | awk '{print "less -S "$1"|awk \047{printf(\"%.12f\\t\",\$NF);print \$0}\047|sed 1d|sort -k1,1n|head -n 174|awk \047{print \$2\":\"\$3}\047 > "$1".174.id"}' | sed 's/pvalue.gwas.result.//2' | sh &
top 0.0001
ls --color=auto | grep --color=auto gwas.result$ | awk '{print "less -S "$1"|awk \047{printf(\"%.12f\\t\",\$NF);print \$0}\047|sed 1d|sort -k1,1n|head -n 348|awk \047{print \$2\":\"\$3}\047 > "$1".348.id"}' | sed 's/pvalue.gwas.result.//2' | sh &

  
perl /groups/dolan/user/shuangyang.wu/script/merge_files.pl -r 2200.id -n 1 -i ../../worldclim/t1.1 -c 1 -o out
less -S out |cut -f2-24|grep -v 'NA' > adaptive.loc.frq
head ../../worldclim/t1.1 
sed -i 's/_alt//g' adaptive.loc.frq

#2191
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/3-lfmm")
a=read.table("adaptive.loc.frq",header = T,sep="\t",row.names=1)
b<-t(a)
write.table(b,"adaptive.loc.frq.t",sep="\t")

sed -i 's/\"//g' adaptive.loc.frq.t

perl /groups/dolan/user/shuangyang.wu/script/merge_files.pl -r ../1-pheno/22.pops.geo.bio1-19.tsv -n 1 -i adaptive.loc.frq.t -c 1 -o out1

mv out1 2200.GF


perl /groups/dolan/user/shuangyang.wu/script/merge_files.pl -r 4058.id -n 1 -i ../../worldclim/t1.1 -c 1 -o out
less -S out |cut -f2-24|grep -v 'NA' > adaptive.loc.frq
head ../../worldclim/t1.1 
sed -i 's/_alt//g' adaptive.loc.frq

#4047
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/3-lfmm")
a=read.table("adaptive.loc.frq",header = T,sep="\t",row.names=1)
b<-t(a)
write.table(b,"adaptive.loc.frq.t",sep="\t")

sed -i 's/\"//g' adaptive.loc.frq.t

perl /groups/dolan/user/shuangyang.wu/script/merge_files.pl -r ../1-pheno/22.pops.geo.bio1-19.tsv -n 1 -i adaptive.loc.frq.t -c 1 -o out1

mv out1 2200.GF









# 2 22 pops over 2 Alt frequency

ls|grep list|awk '{print "vcftools --gzvcf test.vcf.gz --freq --keep "$1" --out "$1}'|sed 's/.list//2'|sh &
  
ls|grep frq|awk '{print "sed -i \047s/{ALLELE:FREQ}/"$1"_ref\\t"$1"_alt/\047 "$1}'|sed 's/.frq_ref/_ref/'|sed 's/.frq_alt/_alt/'|sh

#id
less -S AUS_Grafenegg.frq|cut -f1-2|sed 's/\t/:/' > snp
#Frq
paste -d"\t" AUS_Ehrenhausen.frq AUS_Grafenegg.frq AUS_MiBa.frq AUS_Schubert.frq GER_Bonn.frq IRE_ARDARA.frq IRE_Dublin.frq IRE_LETTERKENNY.frq IRE_MAGHERA.frq JPN_Guha.frq JPN_HakubaMura.frq JPN_Hiro.frq JPN_Mats.frq JPN_MiyaOh.frq JPN_Suwa1.frq JPN_Suwa2.frq JPN_Taimadera.frq JPN_Takaragaike.frq SWZ_Zurich.frq UK_Cam.frq UK_Norwich.frq UK_Oxford.frq |awk '{for (i = 6; i <= 132; i+=6) printf("%s\t", $i); printf("\n")}'|sed -r 's/[A-Z]://g'|sed 's/-nan/NA/g'|sed 's/\t$//' > all.freq.1 &
# all info of alt freq
paste -d"\t" snp all.freq.1 | sed 's/CHROM:POS/snp/' > t1 
mv t1 t1.1

t1.1 is the final table for all

# adaptive loc
/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/vcf2gwas
cat ./test_pheno.part15/Output/Linear_Mixed_Model/bio15/bio15_20231126_154624/best_p-values/p_wald_bio15_mod_sub_pheno.part15_test_top0.0001.csv ./test_pheno.part14/Output/Linear_Mixed_Model/bio14/bio14_20231126_152157/best_p-values/p_wald_bio14_mod_sub_pheno.part14_test_top0.0001.csv ./test_pheno.part19/Output/Linear_Mixed_Model/bio19/bio19_20231127_030350/best_p-values/p_wald_bio19_mod_sub_pheno.part19_test_top0.0001.csv ./test_pheno.part6/Output/Linear_Mixed_Model/bio6/bio6_20231127_064539/best_p-values/p_wald_bio6_mod_sub_pheno.part6_test_top0.0001.csv ./test_pheno.part11/Output/Linear_Mixed_Model/bio11/bio11_20231126_140417/best_p-values/p_wald_bio11_mod_sub_pheno.part11_test_top0.0001.csv ./test_pheno.part13/Output/Linear_Mixed_Model/bio13/bio13_20231126_145244/best_p-values/p_wald_bio13_mod_sub_pheno.part13_test_top0.0001.csv ./test_pheno.part12/Output/Linear_Mixed_Model/bio12/bio12_20231126_142807/best_p-values/p_wald_bio12_mod_sub_pheno.part12_test_top0.0001.csv ./test_pheno.part7/Output/Linear_Mixed_Model/bio7/bio7_20231127_074624/best_p-values/p_wald_bio7_mod_sub_pheno.part7_test_top0.0001.csv ./test_pheno.part5/Output/Linear_Mixed_Model/bio5/bio5_20231127_062115/best_p-values/p_wald_bio5_mod_sub_pheno.part5_test_top0.0001.csv ./test_pheno.part8/Output/Linear_Mixed_Model/bio8/bio8_20231127_080639/best_p-values/p_wald_bio8_mod_sub_pheno.part8_test_top0.0001.csv ./test_pheno.part17/Output/Linear_Mixed_Model/bio17/bio17_20231126_165348/best_p-values/p_wald_bio17_mod_sub_pheno.part17_test_top0.0001.csv ./test_pheno.part2/Output/Linear_Mixed_Model/bio2/bio2_20231127_050203/best_p-values/p_wald_bio2_mod_sub_pheno.part2_test_top0.0001.csv ./test_pheno.part16/Output/Linear_Mixed_Model/bio16/bio16_20231126_162520/best_p-values/p_wald_bio16_mod_sub_pheno.part16_test_top0.0001.csv ./test_pheno.part18/Output/Linear_Mixed_Model/bio18/bio18_20231127_035232/best_p-values/p_wald_bio18_mod_sub_pheno.part18_test_top0.0001.csv ./test_pheno.part4/Output/Linear_Mixed_Model/bio4/bio4_20231127_060601/best_p-values/p_wald_bio4_mod_sub_pheno.part4_test_top0.0001.csv ./test_pheno.part9/Output/Linear_Mixed_Model/bio9/bio9_20231127_083406/best_p-values/p_wald_bio9_mod_sub_pheno.part9_test_top0.0001.csv ./test_pheno.part1/Output/Linear_Mixed_Model/bio1/bio1_20231127_042029/best_p-values/p_wald_bio1_mod_sub_pheno.part1_test_top0.0001.csv ./test_pheno.part10/Output/Linear_Mixed_Model/bio10/bio10_20231126_134222/best_p-values/p_wald_bio10_mod_sub_pheno.part10_test_top0.0001.csv ./test_pheno.part3/Output/Linear_Mixed_Model/bio3/bio3_20231127_054148/best_p-values/p_wald_bio3_mod_sub_pheno.part3_test_top0.0001.csv |sed 's/,/\t/g'|sed 's/\"//g'|grep -v 'allele1'|cut -f3|sort|uniq|sed 's/\\.*//'|sed 's/^/chr/' > 3662.id

perl /groups/dolan/user/shuangyang.wu/script/merge_files.pl -r 3662.id -n 1 -i ../worldclim/t1.1 -c 1 -o out
less -S out |cut -f2-24|grep -v 'NA' > adaptive.loc.frq
head ../worldclim/t1.1 
sed -i 's/_alt//g' adaptive.loc.frq



setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/2-vcf2gwas")
a=read.table("adaptive.loc.frq",header = T,sep="\t",row.names=1)
b<-t(a)
write.table(b,"adaptive.loc.frq.t",sep="\t")

sed -i 's/\"//g' adaptive.loc.frq.t

perl /groups/dolan/user/shuangyang.wu/script/merge_files.pl -r ../worldclim/22.pops.geo.bio1-19.tsv -n 1 -i adaptive.loc.frq.t -c 1 -o out1

ln -s out1 3662.gf.data

#######GF######
setwd("/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/GF")
less -S 3662.gf.data |cut -f1-22 >22.geo
devtools::install_github("AndiKur4/MaizePal")
install.packages("gradientForest", repos="http://R-Forge.R-project.org")
library(gradientForest)
library(MaizePal)

gfData <- read.table("3662.gf.data",header =T,sep="\t",row.names = "Population")  # 加载环境数据、等位基因频率
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

#plot(gf_candidate, plot.type="S", imp.vars= bioclimatic, leg.posn="topright", cex.legend=0.4, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1)))
#plot(gf_candidate, plot.type="Cumulative.Importance", imp.vars= bioclimatic, show.overall=T, legend=T,common.scale=T,leg.posn="topleft", leg.nspecies=5, cex.lab=0.7, cex.legend=0.4,  cex.axis=0.6, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1),omi=c(0,0.3,0,0)))


## Allele turnover functions across the landscape. 

temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand, type=c("Overall"),standardize = T) # al candidate SNPs 
temp_cand_SNP <- cumimp(gf_candidate,predictor = most_cand, type=c("Species"),standardize = T) #each individual candidate allele

# ylim <- NULL
# for(j in 1:length(temp_cand_SNP)){ ylim <- c(ylim,max(temp_cand_SNP[[j]][[2]])) }
# 
# par(mfrow=c(1,2))
# par(mai=c(0.9,0.8,0.4,0))
# plot(temp_cand_overall,type="n",ylim=c(0,0.25),mgp=c(2,0.6,0), ylab="Cumulative importance",xlab= "bio4")
# 
# for(j in 1:length(temp_cand_SNP)){lines(temp_cand_SNP[[j]],col=adjustcolor(MaizePal::maize_pal("RubyGold")[5],alpha.f = 0.6)) }
# 
# lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4)


### Define populations locally adapted to contrasting environments 

pop_turn <- predict(gf_candidate,present[,grep("bio",names(present))]) 
temp <- data.frame(bio=present[,most_cand],imp=pop_turn[,most_cand]) 
warm <- which(pop_turn[,most_cand] >= (mean(pop_turn[,most_cand])))
cold <- which(pop_turn[,most_cand] < (mean(pop_turn[,most_cand]))) 
categories <- list(cold=rownames(pop_turn)[cold],warm=rownames(pop_turn)[warm])
plot(temp_cand_overall,type="n",ylim=c(0,0.25),mgp=c(2,0.6,0),ylab="Cumulative importance",xlab= paste("Temperature Seasonality (standard deviation ×100) (BIO4)",sep=""),main="Adaptive SNPs") 
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
5  7  6  8 10 11  2  1 13 12  4  3  9

rownames(pop_turn)[warm]
[1] "JPN_HakubaMura"  "JPN_Suwa2"       "JPN_Guha"        "JPN_MiyaOh"      "JPN_Taimadera"   "JPN_Suwa1"       "JPN_Hiro"        "JPN_Mats"        "JPN_Takaragaike"

id_w
7 5 9 3 2 6 8 4 1

####phenotype distribution#####

52.197  0.1308

df<-read.table("")
df <- mtcars[1:4]
head(df, n = 5)



library("ggplot2")
library("GGally")


setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/3-lfmm")
dis=read.table("22.pops.dis",header = T)
p_dis=ggpairs(dis, columns = c(6,7,10,12,20), aes(colour = State))

# ggscatter <- function(data, mapping, ...) {
#   x <- GGally::eval_data_col(data, mapping$x)
#   y <- GGally::eval_data_col(data, mapping$y)
#   df <- data.frame(x = x, y = y)
#   sp1 <- ggplot(df, aes(x=x, y=y)) +
#     geom_point() +
#     geom_abline(intercept = 0, slope = 1, col = 'darkred')
#   return(sp1)
# }
# 
# ggdehist <- function(data, mapping, ...) {
#   x <- GGally::eval_data_col(data, mapping$x)
#   df <- data.frame(x = x)
#   dh1 <- ggplot(df, aes(x=x)) +
#     geom_histogram(aes(y=..density..), bins = 50, fill = 'steelblue', color='black', alpha=.4) +
#     geom_density(aes(y=..density..)) + 
#     theme_minimal()
#   return(dh1)
# }
# 
# ggpairs(dis,columns = c(6,7,10,12,20),
#         lower = list(continuous = wrap(ggscatter)),
#         diag = list(continuous = wrap(ggdehist))) + 
#   theme_minimal() +
#   theme(panel.grid = element_blank(),
#         panel.border = element_rect(fill=NA),
#         axis.text =  element_text(color='black'))






###this is final command
ggpairs(dis, columns = c(6,7,10,12,20), aes(colour = State))+theme(panel.grid = element_blank(), panel.border = element_rect(fill=NA),axis.text =  element_text(color='black'))

######annovar#####
/scratch-cbe/users/shuangyang.wu/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/6-GWAS/update/5-annovar
wget https://marchantia.info/download/MpTak_v6.1/MpTak_v6.1r1.gff.gz
module load cufflinks/2.2.1-foss-2018b
gffread MpTak_v6.1r1.gff -T -o my.gtf
module load kent_tools/20190507-linux.x86_64
gtfToGenePred  -genePredExt my.gtf genome_refGene.txt
/groups/dolan/user/shuangyang.wu/software/annovar/retrieve_seq_from_fasta.pl --format refGene --seqfile Tak1.V6.withFullUV.fa genome_refGene.txt --out genome_refGeneMrna.fa
mkdir genome
mv genome* genome
/groups/dolan/user/shuangyang.wu/software/annovar/convert2annovar.pl -format vcf4old ../snpOnlymaxmiss09.recode.vcf  >vcf.annovar.inputi
/groups/dolan/user/shuangyang.wu/software/annovar/annotate_variation.pl --geneanno --neargene 2000 --dbtype refGene -out tes-1 --buildver genome vcf.annovar.input ./genome/


  
less -S bio4.pvalue.gwas.result|awk '{OFS="\t"}{print $0,-(log($NF)/log(10))}'|awk '$NF>9'|sort -k4,4nr|head -n 6
less -S bio5.pvalue.gwas.result|awk '{OFS="\t"}{print $0,-(log($NF)/log(10))}'|awk '$NF>8'|sort -k4,4nr|head -n 10
less -S bio10.pvalue.gwas.result|awk '{OFS="\t"}{print $0,-(log($NF)/log(10))}'|awk '$NF>7.8'|sort -k4,4nr|head -n 4
less -S bio8.pvalue.gwas.result|awk '{OFS="\t"}{print $0,-(log($NF)/log(10))}'|awk '$NF>7.8'|sort -k4,4nr|head -n 9  

egrep '15395471|17861387|15409213|17864867|20208658|4969164' tes-1.variant_function |egrep 'chr5|chr7'|cut -f3-5 > BIO4.select.bed
egrep '15395471|15409213|17861387|17864867|20208658|4969164|3871867|9442297|12356576|20242062' tes-1.variant_function |egrep 'chr5|chr7'|cut -f3-5 > BIO5.select.bed
egrep '17864867|20208658|4969164|17861387' tes-1.variant_function |egrep 'chr5|chr7'|cut -f3-5 > BIO10.select.bed
egrep '15395471|17685722|15409213|6755608|15394670|11067011|11068260|6641326|9653963' tes-1.variant_function |egrep 'chr5|chr7|chr8'|cut -f3-5 > BIO8.select.bed
  
######male female protein coding gene TE fst######
#autosome GENE and TE bed
GENE bed: /groups/dolan/user/shuangyang.wu/ref/
less -S MpTak_v6.1r1.gff.gz|grep 'gene\b'|grep -v 'locus_type'|grep '^chr'|egrep -v 'chrV|chrU'|cut -d";" -f1|awk '{OFS="\t"}{print $1,$4,$5,$NF,$7}'|sed 's/ID=//' > MpTak_v6.1r1.autosome.gene.bed
TE bed:/groups/dolan/user/shuangyang.wu/ref/EDTA-Tak1.V6.withFullUV
less -S TAK1.V6.TE.gtf|grep '^chr'|egrep -v 'chrU|chrV'|cut -d";" -f1|awk '{OFS="\t"}{print $1,$4,$5,$NF,$7}'|sed 's/\"//g' > Tak1.V6.autosome.TE.bed
less -S Tak1.V6.autosome.TE.bed|sort -k1,1 -k2,2n > Tak1.V6.autosome.TE.1.bed
#prepare female and male pops fst in 500 bp and intersect to get TE and gene differentiation
/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/UV-specific/autosome
#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --partition=c
#SBATCH --qos="long"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5G
#SBATCH --time=5-00:00:00
#SBATCH --output=my.stdout.pixy.femalePOP

pixy --stats pi fst dxy --vcf ../../../test_filtered1.vcf.gz --populations U.sample.info --window_size 500  --n_cores 20 --output_prefix female_pixy_500  --bypass_invariant_check yes

#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --partition=c
#SBATCH --qos="long"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5G
#SBATCH --time=5-00:00:00
#SBATCH --output=my.stdout.malePOP

pixy --stats pi fst dxy --vcf ../../../test_filtered1.vcf.gz --populations V.sample.info --window_size 500  --n_cores 20 --output_prefix male_pixy_500  --bypass_invariant_check yes

#post file preparation
less -S male_pixy_500_fst.txt|awk '{if($6<0){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t0\t"$7}else{print $0}}'|awk '{OFS="\t"}{print $3,$4,$5,$6,$6,"+"}'|sed '1d'|grep -v NA > male.pop.autosome.fst.bed
less -S female_pixy_500_fst.txt|awk '{if($6<0){print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t0\t"$7}else{print $0}}'|awk '{OFS="\t"}{print $3,$4,$5,$6,$6,"+"}'|sed '1d'|grep -v NA > female.pop.autosome.fst.bed
less -S female_pixy_500_pi.txt|grep 'Europe'|cut -f2-5 > female.autosome.Europe.pi
less -S female_pixy_500_pi.txt|grep 'Asia'|cut -f2-5 > female.autosome.Asia.pi
less -S male_pixy_500_pi.txt|grep 'Asia'|cut -f2-5 > male.autosome.Asia.pi
less -S male_pixy_500_pi.txt|grep 'Europe'|cut -f2-5 > male.autosome.Europe.pi

#male population
bedtools map -a MpTak_v6.1r1.autosome.gene.bed -b male.pop.autosome.fst.bed -c 5 -o mean|sed 's/\.$/0/'|cut -f1-4,6 > male.autosome.gene.fst.value
bedtools map -a Tak1.V6.autosome.TE.1.bed -b male.pop.autosome.fst.bed -c 5 -o mean|sed 's/\.$/0/'|cut -f1-4,6 > male.autosome.TE.fst.value

#female population
bedtools map -a MpTak_v6.1r1.autosome.gene.bed -b female.pop.autosome.fst.bed -c 5 -o mean|sed 's/\.$/0/'|cut -f1-4,6  > female.autosome.gene.fst.value
bedtools map -a Tak1.V6.autosome.TE.1.bed -b female.pop.autosome.fst.bed -c 5 -o mean|sed 's/\.$/0/'|cut -f1-4,6  > female.autosome.TE.fst.value

cat Gamtolog/V/V.gamte.fst.value Specific/V/V.specific.fst.value >autosome/male-population/male.sexchrome.gene.fst.value
cp V.TE.value ../autosome/male-population/male.sexchrome.TE.fst.value

cat Gamtolog/U/U.gamte.fst.value Specific/U/U.specific.fst.value >autosome/female-population/female.sexchrome.gene.fst.value
cp U.TE.value ../autosome/female-population/female.sexchrome.TE.fst.value

#male population
awk '{print $0"\tGene"}' male.autosome.gene.fst.value > male.autosome.gene.fst.1.value
awk '{print $0"\tTE"}'  male.autosome.TE.fst.value > male.autosome.TE.fst.1.value
awk '{print $0"\tTE"}'  male.sexchrome.TE.fst.value > male.sexchrome.TE.fst.1.value
awk '{print $0"\tGene"}'  male.sexchrome.gene.fst.value > male.sexchrome.gene.fst.1.value

cat male.autosome.gene.fst.1.value male.autosome.TE.fst.1.value male.sexchrome.gene.fst.1.value male.sexchrome.TE.fst.1.value > all.gene.TE.value
less -S all.gene.TE.fst.value |egrep 'chromo|chr1\b|chr2\b|chr3\b|chr4\b|chr5\b|chr6\b|chr7\b|chr8\b|chrV' > all.gene.TE.fst.1.value

/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/UV-specific/autosome/male-population
setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/UV-specific/autosome/male-population")

library(ggplot2)
library(RColorBrewer)
male_pop=read.table("all.gene.TE.fst.1.value",header = T)
compaired <- list(c("Gene", "TE"))
ggplot(male_pop, aes(Type,Fst,fill=Type),guide="none")+geom_violin(draw_quantiles = c(0.5))+labs(y = "Fst") +geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test)+scale_fill_brewer(palette = "Set2")+facet_wrap(chromosome~.,ncol = 3,nrow = 3)+theme_bw()

#female population

#male population
awk '{print $0"\tGene"}' female.autosome.gene.fst.value > male.autosome.gene.fst.1.value
awk '{print $0"\tTE"}'  female.autosome.TE.fst.value > male.autosome.TE.fst.1.value
awk '{print $0"\tTE"}'  female.sexchrome.TE.fst.value > male.sexchrome.TE.fst.1.value
awk '{print $0"\tGene"}'  female.sexchrome.gene.fst.value > male.sexchrome.gene.fst.1.value

cat male.autosome.gene.fst.1.value male.autosome.TE.fst.1.value male.sexchrome.gene.fst.1.value male.sexchrome.TE.fst.1.value |egrep 'chromo|chr1\b|chr2\b|chr3\b|chr4\b|chr5\b|chr6\b|chr7\b|chr8\b|chrU' > all.gene.TE.fst.1.value

setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/UV-specific/autosome/female-population")
female_pop=read.table("all.gene.TE.fst.1.value",header = T)
compaired <- list(c("Gene", "TE"))
library(ggsignif)
ggplot(female_pop, aes(Type,Fst,fill=Type),guide="none")+geom_violin(draw_quantiles = c(0.5))+labs(y = "Fst") +geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test)+scale_fill_brewer(palette = "Set2")+facet_wrap(chromosome~.,ncol = 3,nrow = 3)+theme_bw()

#male population
awk '{print $0"\tGene"}' male.autosome.gene.fst.value > male.autosome.gene.fst.1.value
awk '{print $0"\tTE"}'  male.autosome.TE.fst.value > male.autosome.TE.fst.1.value
awk '{print $0"\tTE"}'  male.sexchrome.TE.fst.value > male.sexchrome.TE.fst.1.value
awk '{print $0"\tGene"}'  male.sexchrome.gene.fst.value > male.sexchrome.gene.fst.1.value

cat male.autosome.gene.fst.1.value male.autosome.TE.fst.1.value male.sexchrome.gene.fst.1.value male.sexchrome.TE.fst.1.value > all.gene.TE.value
less -S all.gene.TE.fst.value |egrep 'chromo|chr1\b|chr2\b|chr3\b|chr4\b|chr5\b|chr6\b|chr7\b|chr8\b|chrV' > all.gene.TE.fst.1.value

/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/UV-specific/autosome/male-population
setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/UV-specific/autosome/male-population")

library(ggplot2)
library(RColorBrewer)
male_pop=read.table("all.gene.TE.fst.1.value",header = T)
compaired <- list(c("Gene", "TE"))
ggplot(male_pop, aes(Type,Fst,fill=Type),guide="none")+geom_violin(draw_quantiles = c(0.5))+labs(y = "Fst") +geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test)+scale_fill_brewer(palette = "Set2")+facet_wrap(chromosome~.,ncol = 3,nrow = 3)+theme_bw()

#female population

#male population
awk '{print $0"\tGene"}' female.autosome.gene.fst.value > male.autosome.gene.fst.1.value
awk '{print $0"\tTE"}'  female.autosome.TE.fst.value > male.autosome.TE.fst.1.value
awk '{print $0"\tTE"}'  female.sexchrome.TE.fst.value > male.sexchrome.TE.fst.1.value
awk '{print $0"\tGene"}'  female.sexchrome.gene.fst.value > male.sexchrome.gene.fst.1.value

cat male.autosome.gene.fst.1.value male.autosome.TE.fst.1.value male.sexchrome.gene.fst.1.value male.sexchrome.TE.fst.1.value |egrep 'chromo|chr1\b|chr2\b|chr3\b|chr4\b|chr5\b|chr6\b|chr7\b|chr8\b|chrU' > all.gene.TE.fst.1.value

setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/UV-specific/autosome/female-population")
female_pop=read.table("all.gene.TE.fst.1.value",header = T)
compaired <- list(c("Gene", "TE"))
ggplot(female_pop, aes(Type,Fst,fill=Type),guide="none")+geom_violin(draw_quantiles = c(0.5))+labs(y = "Fst") +geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test)+scale_fill_brewer(palette = "Set2")+facet_wrap(chromosome~.,ncol = 3,nrow = 3)+theme_bw()

#another figure with gene and TE split
# ggplot(female_pop, aes(chromosome,Fst,fill=Type),guide="none")+geom_violin(draw_quantiles = c(0.5))+labs(y = "Fst")+scale_fill_brewer(palette = "Set2")+facet_wrap(Type~.)+theme_bw()
# ggplot(male_pop, aes(chromosome,Fst,fill=Type),guide="none")+geom_violin(draw_quantiles = c(0.5))+labs(y = "Fst")+scale_fill_brewer(palette = "Set2")+facet_wrap(Type~.)+theme_bw()
# 
# ggplot(male_pop, aes(chromosome,Fst,fill=Type),guide="none")+geom_boxplot()+labs(y = "Fst")+scale_fill_brewer(palette = "Set2")+facet_wrap(Type~.)+theme_bw()
ggplot(male_pop, aes(chromosome,Fst,fill=Type),guide="none")+geom_boxplot(position = position_dodge(width = 0.9))+scale_fill_brewer(palette = "Set2")+stat_summary(fun.y = median,geom = 'line',aes(group = Type),position = position_dodge(width = 0.9))+facet_wrap(Type~.)+theme_bw()
ggplot(female_pop, aes(chromosome,Fst,fill=Type),guide="none")+geom_boxplot(position = position_dodge(width = 0.9))+scale_fill_brewer(palette = "Set2")+stat_summary(fun.y = median,geom = 'line',aes(group = Type),position = position_dodge(width = 0.9))+facet_wrap(Type~.)+theme_bw()

#final figure facet
library(PupillometryR)
ggplot(male_pop, aes(chromosome,Fst,fill=Type),guide="none")+geom_flat_violin(position = position_nudge(x = .2), alpha = .4)+geom_point(aes(color = Type),position = position_jitter(w = .15), size = 1,alpha = 0.4,show.legend = F) +geom_boxplot(width = .25,outlier.shape = NA,alpha = 0.5)+scale_fill_brewer(palette = "Set2")+scale_color_brewer(palette = "Set2")+stat_summary(fun.y = median,geom = 'line',aes(group = Type),position = position_dodge(width = 0.9))+facet_wrap(Type~.)+theme_bw()

ggplot(female_pop, aes(chromosome,Fst,fill=Type),guide="none")+geom_flat_violin(position = position_nudge(x = .2), alpha = .4)+geom_point(aes(color = Type),position = position_jitter(w = .15), size = 1,alpha = 0.4,show.legend = F) +geom_boxplot(width = .25,outlier.shape = NA,alpha = 0.5)+scale_fill_brewer(palette = "Set2")+scale_color_brewer(palette = "Set2")+stat_summary(fun.y = median,geom = 'line',aes(group = Type),position = position_dodge(width = 0.9))+facet_wrap(Type~.)+theme_bw()


## 1 autosome pattern with sex chromosome pattern
less -S all.gene.TE.fst.1.value |awk '{if($1~/chrV/){print $0"\tV"}else{print $0"\tAutosome"}}' > all.gene.TE.fst.2.value
New_male=read.table("all.gene.TE.fst.2.value",header = T)
compaired <- list(c("Autosome", "V"))

#median value
ggplot(New_male, aes(Category,Fst,fill=Category),guide="none")+geom_boxplot(position = position_dodge(width = 0.9))+scale_fill_npg()+geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test)+theme_bw()

#distribution pattern different 
ggplot(New_male, aes(Fst,Category,fill=Category), guide="none") +geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +theme_bw()+theme(axis.text.y=element_blank())+scale_fill_npg() +labs(x="Fst", y="Density")+xlim(0,0.5)

##redo whole genome level fst distribution##
nohup pixy --stats pi fst dxy --vcf test_filtered1.vcf.gz --populations U.sample.info --window_size 20000  --n_cores 10 --output_prefix female_pixy_20K  --bypass_invariant_check yes &
nohup pixy --stats pi fst dxy --vcf test_filtered1.vcf.gz --populations V.sample.info --window_size 20000  --n_cores 10 --output_prefix male_pixy_20K  --bypass_invariant_check yes &
  
less -S female_pixy_20K_fst.txt|egrep -v 'chrV'|awk '{OFS="\t"}{print $3,$6}'|awk '{OFS="\t"}{if($2<0){print $1,"0"}else{print $0}}'|grep -v NA|awk '{if($1~/chrU/){print $0"\tU"}else{print $0"\tAutosome"}}' > female_genome_level_autosome_U.fst.table

less -S male_pixy_20K_fst.txt|egrep -v 'chrU'|awk '{OFS="\t"}{print $3,$6}'|awk '{OFS="\t"}{if($2<0){print $1,"0"}else{print $0}}'|grep -v NA|awk '{if($1~/chrV/){print $0"\tV"}else{print $0"\tAutosome"}}' > male_genome_level_autosome_V.fst.table

setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/UV-specific/autosome")
Male_autso_sex=read.table("male_genome_level_autosome_V.fst.table",head=T)
library(ggsignif)
library("PupillometryR")
library(ggsci)
compaired <- list(c("Autosome", "V"))
library(ggridges)
#All autosomes ridege plot
ggplot(Male_autso_sex, aes(Fst,Category,fill=Category), guide="none") +geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5) +theme_bw()+theme(axis.text.y=element_blank())+scale_fill_npg() +labs(x="Fst", y="Density")+theme_bw()


#All autosomes boxplot with sign test
ggplot(Male_autso_sex, aes(Category,Fst,fill=Category),guide="none")+geom_boxplot(position = position_dodge(width = 0.9))+scale_fill_npg()+geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test)+theme_bw()

#All autosome boxplot violin and point#
ggplot(Male_autso_sex) +
  aes(x =Category , 
      y = Fst, 
      fill = Category) +
  geom_flat_violin(position = position_nudge(x = .2), 
                   alpha = .4) +
  geom_point(aes(color = Category), 
             position = position_jitter(w = .15), 
             size = 1,
             alpha = 0.4,
             show.legend = F) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5)+geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test)+scale_fill_npg()+scale_color_npg()+theme_bw()


#split autosome to different chromosome

#ridge
ggplot(Male_autso_sex, aes(Fst,chromosome,fill=chromosome), guide="none") +geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5)+scale_fill_npg()+theme_bw()

#ggplot(Male_autso_sex, aes(Fst,Category,fill=Category), guide="none") +geom_density_ridges(quantile_lines=TRUE, quantile_fun=function(x,...)median(x), size=0.5)+theme(axis.text.y=element_blank())+scale_fill_npg() +geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test)+theme_bw()
#box cloud point
ggplot(Male_autso_sex) +
  aes(x =chromosome , 
      y = Fst, 
      fill = chromosome) +
  geom_flat_violin(position = position_nudge(x = .2), 
                   alpha = .4) +
  geom_point(aes(color = Category),  
             position = position_jitter(w = .15), 
             size = 1,
             alpha = 0.4,
             show.legend = F) +
  geom_boxplot( width = .25, 
               outlier.shape = NA,
               alpha = 0.5)+scale_fill_npg()+theme_bw()
#for each autosome
ggplot(Male_autso_sex) +
  aes(x =chromosome , 
      y = Fst, 
      fill = chromosome) +
  geom_flat_violin(position = position_nudge(x = .2), 
                   alpha = .4) +
  geom_point(aes(color = chromosome), 
             position = position_jitter(w = .15), 
             size = 1,
             alpha = 0.4,
             show.legend = F) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5)+scale_fill_npg()+scale_color_npg()+stat_summary(fun.y = median,geom = 'line',aes(group = NA),position = position_dodge(width = 0.9))+theme_bw()



#female 
Female_autso_sex=read.table("female_genome_level_autosome_U.fst.table",head=T)
compaired <- list(c("Autosome", "U"))
ggplot(Female_autso_sex) +
  aes(x =Category , 
      y = Fst, 
      fill = Category) +
  geom_flat_violin(position = position_nudge(x = .2), 
                   alpha = .4) +
  geom_point(aes(color = Category), 
             position = position_jitter(w = .15), 
             size = 1,
             alpha = 0.4,
             show.legend = F) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5)+geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test)+scale_fill_npg()+theme_bw()

#Test same number female or male
setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/UV-specific/autosome/Test-same-number")
Female_autso_sex=read.table("female_genome_level_autosome_U.fst.table",head=T)
compaired <- list(c("Autosome", "U"))
p8=ggplot(Female_autso_sex) +
  aes(x =Category , 
      y = Fst, 
      fill = Category) +
  geom_flat_violin(position = position_nudge(x = .2), 
                   alpha = .4) +
  geom_point(aes(color = Category), 
             position = position_jitter(w = .15), 
             size = 1,
             alpha = 0.4,
             show.legend = F) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5)+geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test)+scale_fill_npg()+theme_bw()

male_autso_sex=read.table("male_genome_level_autosome_V.fst.table",head=T)
compaired <- list(c("Autosome", "V"))
p7=ggplot(male_autso_sex) +
  aes(x =Category , 
      y = Fst, 
      fill = Category) +
  geom_flat_violin(position = position_nudge(x = .2), 
                   alpha = .4) +
  geom_point(aes(color = Category), 
             position = position_jitter(w = .15), 
             size = 1,
             alpha = 0.4,
             show.legend = F) +
  geom_boxplot(width = .25, 
               outlier.shape = NA,
               alpha = 0.5)+geom_signif(comparisons = compaired,step_increase = 0.3,map_signif_level = F,test = wilcox.test)+scale_fill_npg()+theme_bw()





setwd("/groups/dolan/user/shuangyang.wu/CBE/Backup-data-9/rudealis/All-data/new/gvcf-209/78-no-dup-fq/old-assembly-fullUV/3-Pi-Fst-1/Euro-Japan/UV-specific/autosome")
corrFemale=read.table("corrTEGene.female",header = T)
head(corrFemale,8)
cor_value <- round(cor(corrFemale$Gene, corrFemale$TE), 2)
p1=ggplot(corrFemale, aes(x = Gene, y = TE)) +
  geom_point() +
  stat_smooth(method = "lm", level = 0.95, fill = "lightblue", alpha = 0.3)+geom_text(aes(label = "r=0.85"), x = 0.1, y = 0.125)+xlab("Autosome+U Gene Fst")+ylab("Autosome+U TE Fst")+theme_bw()

#without sex chr
withoutU=head(corrFemale,8)
cor_value <- round(cor(withoutU$Gene, withoutU$TE), 2)
p2=ggplot(withoutU, aes(x = Gene, y = TE)) +
  geom_point() +
  stat_smooth(method = "lm", level = 0.95, fill = "lightblue", alpha = 0.3)+geom_text(aes(label = "r=0.86"), x = 0.1, y = 0.125)+xlab("Autosome Gene Fst")+ylab("Autosome TE Fst")+theme_bw()


corrMale=read.table("corrTEGene.male",header = T)
cor_value <- round(cor(corrMale$Gene, corrMale$TE), 2)
p3=ggplot(corrMale, aes(x = Gene, y = TE)) +
  geom_point() +
  stat_smooth(method = "lm", level = 0.95, fill = "lightblue", alpha = 0.3)+geom_text(aes(label = "r=0.22"), x = 0.06, y = 0.06)+xlab("Autosome+V Gene Fst")+ylab("Autosome+V TE Fst")+theme_bw()

#without sex chr
withoutV=head(corrMale,8)
cor_value <- round(cor(withoutV$Gene, withoutV$TE), 2)
p4=ggplot(withoutV, aes(x = Gene, y = TE)) +
  geom_point() +
  stat_smooth(method = "lm", level = 0.95, fill = "lightblue", alpha = 0.3)+geom_text(aes(label = "r=0.85"), x = 0.06, y = 0.06)+xlab("Autosome Gene Fst")+ylab("Autosome TE Fst")+theme_bw()






