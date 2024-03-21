library(ggplot2)
library(dplyr)
library(reshape2)

#load data
mash_weights_all <- read.csv("./mash_weights_all_noReplace.csv")
mash_weights_0.05 <- read.csv("./mash_weights_p0.05_Replace.csv")


# create column of difference of male and female mash weights
mash_weights <- mash_weights_all %>%
  mutate(diff_fm = sum_weight_f - sum_weight_m) %>%     # amplification difference
  select(c(1,7))

#plot
sum_w_v = ggplot(data=mash_weights_all, aes(x=sum_weight_f, y= sum_weight_m)) +
geom_point(aes(col=`cancer`),size=2) +
geom_abline(intercept = 0, slope = 1,color="red",linetype=2)+
labs(x='Female-biased amplification', y = "Male-biased amplification") + 
		scale_fill_manual(values=c(BLCA= "#AABB66",COAD ="#EEEE00",HNSC="#FFD700",LIHC= "#CC9955",KIRC= "#FF6600",KIRP= "#FFAA00",LUAD ="#006600",LUSC= "#FF00BB",THCA ="#552200")) +  #black NAT; red:LUSCs
		scale_colour_manual(values=c(BLCA= "#AABB66",COAD ="#EEEE00",HNSC="#FFD700",LIHC= "#CC9955",KIRC= "#FF6600",KIRP= "#FFAA00",LUAD ="#006600",LUSC= "#FF00BB",THCA ="#552200") ) + 
		theme_Publication()+ xlim(0,100)+ylim(0,100)+
		theme(legend.position = "none")+
		 geom_text_repel(aes(label=cancer),  max.overlaps = Inf, box.padding=0.7, segment.color="grey")



###proportion of sex-effect
mash_weights_all$sum1= rowSums(mash_weights_all[,2:4])
dat <- melt(mash_weights_all[,c("cancer","sum_weight_f","sum_weight_m","sum_weight_e","sum1")],id.var=c("cancer","sum1"))
dat$weight = dat$value/dat$sum1*100
dat$variable=as.character(dat$variable)
dat$variable[dat$variable=="sum_weight_f"]="female > male"
dat$variable[dat$variable=="sum_weight_e"]="female = male"
dat$variable[dat$variable=="sum_weight_m"]="female < male"
dat$variable=factor(dat$variable,levels=c("female > male","female < male","female = male"))

mag_plot <- 
  ggplot(dat, aes(x=cancer, y=weight, fill=variable)) +
  geom_bar(position="stack", stat="identity") +
  theme_Publication() +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL, y="Weight", fill="Magnitude") + 
  scale_fill_manual(values = c("female > male"="#AF382F","female < male"="#0D71A8","female = male"="#33CCCC")) +
  theme( axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5),legend.position="right")
print(mag_plot)
dev.off()


###proportion of correlation
mash_weights_all$sum2= rowSums(mash_weights_all[,5:8])
dat <- melt(mash_weights_all[,c("cancer","perf_corr","neg_corr","partial_pos","uncor","sum2" )],id.var=c("cancer","sum2"))
dat$weight = dat$value/dat$sum2*100
dat$variable=as.character(dat$variable)
dat$variable[dat$variable=="perf_corr"]="perfect"
dat$variable[dat$variable=="neg_corr"]="negative"
dat$variable[dat$variable=="partial_pos"]="partial"
dat$variable[dat$variable=="uncor"]="uncorrelated"
dat$variable=factor(dat$variable,levels=c("perfect","partial","uncorrelated","negative"))

colors <- c("#d1b724", "#563f61", "#b0464f", "#2b62d9")

cor_plot <- 
  ggplot(dat, aes(x=cancer, y=weight, fill=variable)) +
  geom_bar(position="stack", stat="identity") +
  theme_Publication() +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=NULL, y="Weight", fill="Correlation") + 
  scale_fill_manual(values = colors) +
  theme( axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5),legend.position="right")
print(cor_plot)


####for non-sex effect
THCA <- read.csv("THCA_mixprop_91_all_noReplace.txt", sep="\t")
BLCA <- read.csv("BLCA_mixprop_77_all_noReplace.txt", sep="\t")
LIHC<-read.csv("LIHC_mixprop_92_all_noReplace.txt", sep="\t")
LUAD <- read.csv("LUAD_mixprop_35_all_noReplace.txt", sep="\t")
LUSC <- read.csv("LUSC_mixprop_19_all_noReplace.txt", sep="\t")
COAD <- read.csv("COAD_mixprop_28_all_noReplace.txt",sep="\t")
KIRC <- read.csv("KIRC_mixprop_25_all_noReplace.txt", sep="\t")
KIRP <-read.csv("KIRP_mixprop_8_all_noReplace.txt", sep="\t")
HNSC <-read.csv("HNSC_mixprop_13_all_noReplace.txt", sep="\t")

prepare_df <- function(df) {
  df$magnitude <- factor(df$magnitude, levels = c('f1','f3', 'f2', 'f1.5','equal1','m1.5','m2','m3','m1'))
  df <- df %>% mutate_at(1, as.numeric) %>%
    arrange(correlation, magnitude)
  return(df)
}

df_null_cp <- data.frame()

for(ca in cancer_list){
	df <- get(ca)
	df_values <- data.frame( Name = df$X, Mean = rowMeans(df[2:ncol(df)]), 
                         SE = rowSds(as.matrix(df[2:ncol(df)])) / sqrt(length(colnames(df)) - 1) )
	# split matrice names
	df_values <- df_values %>%
  separate(Name, c("sex","correlation","magnitude"), sep="[_]", fill="right") %>%
  mutate(magnitude = paste0(sex, magnitude))


# split between null and values
	df_ave <- prepare_df(df_values[2:nrow(df_values),c(2,3,4,5)])
	df_null <- prepare_df(df_values[1,c(2,3,4,5)])
	df_null <- df_null %>% 
 	 mutate(Mean = as.numeric(Mean)) %>%
	  mutate(mean_lab = ifelse(Mean < 0.0005, "0%", sprintf("%.1f%%", round(Mean*100,1)) )) %>%
 	 mutate(se_lab = ifelse(SE < 0.0005, "0%", sprintf("%.1f%%", round(SE*100,1)) ))
	df_null$cancer <- ca
	df_null_cp <- rbind(df_null_cp, df_null)

}

df_null_cp$p_thresh <- "1"

THCA0.05 <- read.csv("./mash_p0.05/THCA_mixprop_100_all.txt",sep="\t")
BLCA0.05 <- read.csv("./mash_p0.05/BLCA_mixprop_100_all.txt",sep="\t")
LIHC0.05 <- read.csv("./mash_p0.05/LIHC_mixprop_100_all.txt",sep="\t")
LUAD0.05 <- read.csv("./mash_p0.05/LUAD_mixprop_100_all.txt",sep="\t")
LUSC0.05 <- read.csv("./mash_p0.05/LUSC_mixprop_100_all.txt",sep="\t")
COAD0.05 <- read.csv("./mash_p0.05/COAD_mixprop_100_all.txt",sep="\t")
KIRC0.05 <- read.csv("./mash_p0.05/KIRC_mixprop_100_all.txt",sep="\t")
KIRP0.05 <- read.csv("./mash_p0.05/KIRP_mixprop_100_all.txt",sep="\t")
HNSC0.05 <- read.csv("./mash_p0.05/HNSC_mixprop_100_all.txt",sep="\t")

for(ca in cancer_list){
	df <- get(paste0(ca,"0.05"))
	df_values <- data.frame( Name = df[,1], Mean = rowMeans(df[2:ncol(df)]), 
                         SE = rowSds(as.matrix(df[2:ncol(df)])) / sqrt(length(colnames(df)) - 1) )
	# split matrice names
	df_values <- df_values %>%
  separate(Name, c("sex","correlation","magnitude"), sep="[_]", fill="right") %>%
  mutate(magnitude = paste0(sex, magnitude))


# split between null and values
	df_ave <- prepare_df(df_values[2:nrow(df_values),c(2,3,4,5)])
	df_null <- prepare_df(df_values[1,c(2,3,4,5)])
	df_null <- df_null %>% 
 	 mutate(Mean = as.numeric(Mean)) %>%
	  mutate(mean_lab = ifelse(Mean < 0.0005, "0%", sprintf("%.1f%%", round(Mean*100,1)) )) %>%
 	 mutate(se_lab = ifelse(SE < 0.0005, "0%", sprintf("%.1f%%", round(SE*100,1)) ))
	df_null$cancer <- ca
	df_null$p_thresh <-"0.05"
	df_null_cp <- rbind(df_null_cp, df_null)

}


write.csv(df_null_cp, file="./df_null_cp_allCpGs_nominal_weights.csv",row.names=F)

#plot
null_plot <- ggplot(df_null_cp, aes(x=p_thresh, y=Mean, group=cancer, color=cancer)) +
  geom_point() +
  geom_line() +
  labs(x="P-value Threshold", y="Weights", title="Weights on \nno-sex effects matrix") +
  scale_colour_manual(values=c("BLCA"="#AABB66" , "COAD"="#EEEE00", "HNSC"= "#FFD700", "LIHC"= "#CC9955", "KIRC"="#FF6600", "KIRP"="#FFAA00", "LUAD"="#006600", "LUSC"="#FF00BB", "THCA"="#552200"))+
  theme_Publication() +
  theme(legend.position="none", title=element_text(size=11))
  
ggsave("./weights_nonSexEffect.pdf",h=3.3,w=2.5)

