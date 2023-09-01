#/lustre/home/llwang/TCR/TRUST/Muti_region_new/Merge_all_airr/all_sample_new
#https://shazam.readthedocs.io/en/stable/topics/findThreshold/

library(data.table)
library(dplyr)
library(shazam)
library(ggplot2)

dir="/lustre/home/llwang/TCR/TRUST/Muti_region_new/Merge_all_airr/"
a=fread(paste0(dir,"merged.alldata.tsv")) %>% as.data.frame()
b=fread(paste0(dir,"merged.alldata_igblast_db-pass_germ-pass.tsv")) %>% as.data.frame()

b=b[,c(grep("sequence_id",colnames(b)),which(!colnames(b) %in% colnames(a)))]
a[1:5,1:5]
b[1:5,1:5]

x=merge(a,b,by="sequence_id")

write.table(x,file=paste0(dir,"merged.alldata_igblast_db-pass_germ-pass.add_c_call.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
x=fread("merged.alldata_igblast_db-pass.add_c_call.tsv")
dist_ham <- distToNearest(x, sequenceColumn="junction", 
                          vCallColumn="v_call", jCallColumn="j_call",
                          model="ham", normalize="len", nproc=6)
write.csv(dist_ham,file=paste0(dir,"dist_ham.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
# Find threshold using density method
pdf(paste0(dir,"dist_density_threshold.new.pdf"),4,3)

p1 <- ggplot(subset(dist_ham, !is.na(dist_nearest)),
             aes(x=dist_nearest)) + 
    theme_bw() + 
    xlab("Hamming distance") + 
    ylab("Count") +
    scale_x_continuous(breaks=seq(0, 1, 0.1)) +
    geom_histogram(color="white", binwidth=0.02) +
    geom_vline(xintercept=0.12, color="firebrick", linetype=2)
plot(p1)

output <- findThreshold(dist_ham$dist_nearest, method="density")
threshold <- output@threshold

print(paste0("threshold is",threshold))
# Plot distance histogram, density estimate and optimum threshold
plot(output, title="Density Method")
dev.off()
