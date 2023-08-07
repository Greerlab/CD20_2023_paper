#install.packages("ggplot2")
library(ggplot2)
#########################################################################################
setwd("/Users/plg/Desktop/HaoChing") #set directory to the place you put the 
X=read.csv("6C-0728-1.csv", header = T, sep=",") #read the file you want to analyze
Projectname="6C-0728-1"  # This name will appear in the title and the file name
#########################################################################################
X$frames=NULL # remove the "frame" column
F0=apply(X[450:500,],2,mean) # generate F0
delF=sweep(X,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

Celltrace=function(G){
Figure=ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y=paste(G))) +
           labs(x="frame", y=expression(paste(Delta,"F/F0")))+
           scale_x_continuous(breaks=seq(0,1250,250))+
           geom_vline(xintercept=500, linetype="dashed", color = "red")+
           geom_line(color="blue")+
           ggtitle(paste(Projectname,G,sep = "_"))+
           theme_classic()+
           theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste(Projectname,"_",G,".png",sep = ""))
assign(paste(G), Figure, envir = globalenv())
}# function to create plots and save it 

lapply(cells,Celltrace) # apply the function Celltrace to all the target in the list





