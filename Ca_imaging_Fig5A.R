library(openxlsx)
library(ggplot2)
## 2,5-DMP

File = "ca imaging_isoamyl acetate/20190403.xlsx"
raw1 = read.xlsx(File, sheet = "20190403-CD20-2,5-DMP-1")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell40"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/2,5-DMP.pdf", width = 5, height = 5)

#17 ms4a screening___20211119_Ms4a2-sheet-20211119_well2_mCD20_cell12.pdf
File = "17 ms4a screening/20211119_Ms4a2.xlsx"
raw1 = read.xlsx(File, sheet = "20211119_well2_mCD20")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell12"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/2,5-DMP_1.pdf", width = 5, height = 5)

#20211001, 1007, 1012_Human MS4A1 (WT)_2,5-DMP-sheet-20211001_mCD20_2,5-DMP_well5_cell7
File = "20211001, 1007, 1012_Human MS4A1 (WT)_2,5-DMP.xlsx"
raw1 = read.xlsx(File, sheet = "20211001_mCD20_2,5-DMP_well5")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell7"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/2,5-DMP_2.pdf", width = 5, height = 5)

#20211001, 1007, 1012_Human MS4A1 (WT)_2,5-DMP-sheet-20211007_mCD20_2,5-DMP_well5_cell48
File = "20211001, 1007, 1012_Human MS4A1 (WT)_2,5-DMP.xlsx"
raw1 = read.xlsx(File, sheet = "20211007_mCD20_2,5-DMP_well5")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell48"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/2,5-DMP_3.pdf", width = 5, height = 5)

#ca imaging_2,5-DMP___20181016_2,5-dmp-sheet-1016_CD20_cell2
File = "ca imaging_2,5-DMP/20181016_2,5-dmp.xlsx"
raw1 = read.xlsx(File, sheet = "1016_CD20")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell2"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/2,5-DMP_4.pdf", width = 5, height = 5)

#ca imaging_2,5-DMP___20181016-sheet-1016_CD20_cell2.pdf
File = "ca imaging_2,5-DMP/20181016.xlsx"
raw1 = read.xlsx(File, sheet = "1016_CD20")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell2"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/2,5-DMP_5.pdf", width = 5, height = 5)

#ca imaging_2,5-DMP___20181107-HEK_2,5-dmp-sheet-1107-CD20-1_cell12
File = "ca imaging_2,5-DMP/20181107-HEK_2,5-dmp.xlsx"
raw1 = read.xlsx(File, sheet = "1107-CD20-1")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell12"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/2,5-DMP_5.pdf", width = 5, height = 5)

#ca imaging_acetophenone___0313-sheet-20190313-CD20-2,5-dmp-1_cell97.pdf
File = "ca imaging_acetophenone/0313.xlsx"
raw1 = read.xlsx(File, sheet = "20190313-CD20-2,5-dmp-1")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell97"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/2,5-DMP_5.pdf", width = 5, height = 5)


p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell68"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_classic()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/scale.pdf", width = 5, height = 5)


## 2,3-DMP
File = "ca imaging_vanillin/20190228.xlsx"
raw1 = read.xlsx(File, sheet = "20190228-CD20-2,3-DMP-1")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell58"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/2,3-DMP.pdf", width = 5, height = 5)


## 2,6-DMP
File = "20210113_2,6-DMP.xlsx"
raw1 = read.xlsx(File, sheet = "20210113_well3_CD20_2,6-DMP")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell45"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/2,6-DMP.pdf", width = 5, height = 5)



#pyridine
#ca imaging_N cyclics___20190612-sheet-20190612_pyridine_CD20-2_cell31
File = "ca imaging_N cyclics/20190612.xlsx"
raw1 = read.xlsx(File, sheet = "20190612_pyridine_CD20-2")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell31"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/pyridine.pdf", width = 5, height = 5)

#quinoline
#ca imaging_N cyclics___20190612-sheet-20190612_quinoline_CD20-2_cell8.pdf
File = "ca imaging_N cyclics/20190612.xlsx"
raw1 = read.xlsx(File, sheet = "20190612_quinoline_CD20-2")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell8"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/quinoline.pdf", width = 5, height = 5)

#RS
#20210116_2,6-DMP-sheet-20210116_well7_CD20_RS_cell73
File = "20210116_2,6-DMP.xlsx"
raw1 = read.xlsx(File, sheet = "20210116_well7_CD20_RS")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell73"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/RS.pdf", width = 5, height = 5)


#ca imaging_vanillin/0308.xlsx/20190308-CD20-vanillin-1_cell21
File = "ca imaging_vanillin/0308.xlsx"
raw1 = read.xlsx(File, sheet = "20190308-CD20-vanillin-1")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell21"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/vanillin.pdf", width = 5, height = 5)


#ca imaging_vanillin/0308.xlsx/20190312-CD20-isoamylacetate-1_cell18
File = "ca imaging_isoamyl acetate/0312.xlsx"
raw1 = read.xlsx(File, sheet = "20190312-CD20-isoamylacetate-1")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell18"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/IAA.pdf", width = 5, height = 5)

#ca imaging_N cyclics___20190607-sheet-20190607_CD20_pyrrolidine-1_cell8
File = "ca imaging_N cyclics/20190607.xlsx"
raw1 = read.xlsx(File, sheet = "20190607_CD20_pyrrolidine-1")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell18"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 5A (Ca imaging)/pyrrolidine.pdf", width = 5, height = 5)
