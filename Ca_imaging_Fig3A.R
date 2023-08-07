library(openxlsx)
library(ggplot2)
#Figure 3A
## MS4A1 2,5-DMP
File = "17 ms4a screening/20211203_Ms4a2,6D.xlsx"
raw1 = read.xlsx(File, sheet = "20211203_well6_cd20")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell68"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 3A (Ca imaging)/Ms4a1_2,5-DMP.pdf", width = 5, height = 5)

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell68"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_classic()
ggsave(plot = p, "../CD20 paper figures/Fig 3A (Ca imaging)/scale.pdf", width = 5, height = 5)


## MS4A1 RS
File = "20210113_2,6-DMP.xlsx"
raw1 = read.xlsx(File, sheet = "20210113_well7_CD20_RS")
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
ggsave(plot = p, "../CD20 paper figures/Fig 3A (Ca imaging)/Ms4a1_RS.pdf", width = 5, height = 5)



## MS4A6C 2,5-DMP
File = "ca imaging_2,5-DMP/20181210.xlsx"
raw1 = read.xlsx(File, sheet = "20181210_HEK_2,5-DMP_6C_1")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell4"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 3A (Ca imaging)/Ms4a6c_2,5-DMP.pdf", width = 5, height = 5)

## MS4A6C RS
File = "20200622.xlsx"
raw1 = read.xlsx(File, sheet = "well1_6C_NI")
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
ggsave(plot = p, "../CD20 paper figures/Fig 3A (Ca imaging)/Ms4a6c_RS.pdf", width = 5, height = 5)

# 4B 2,5-DMP
File = "ca imaging_2,5-DMP/20181107-HEK_2,5-dmp.xlsx"
raw1 = read.xlsx(File, sheet = "1107-4B-1")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell1"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
ggsave(plot = p, "../CD20 paper figures/Fig 3A (Ca imaging)/Ms4a4b_2,5-DMP.pdf", width = 5, height = 5)


# 4B RS
File = "20200701.xlsx"
raw1 = read.xlsx(File, sheet = "well9_4B_RS")
raw1 = raw1[,grep("cell",colnames(raw1)), drop = F]
raw1 = raw1[-(1:329),]
F0=apply(raw1[1:30,],2,mean) # generate F0
delF=sweep(raw1,2,F0,"-") #generate delta F
F_F0=sweep(delF,2,F0,"/") # delta F / F0
cells=colnames(F_F0) # list of input, which are the column names of F_F0: cell1, cell2, cell3...

p = ggplot(data=F_F0, aes_string(x=c(1:nrow(F_F0)), y="cell4"))+
  scale_x_continuous(breaks=seq(0,260,25))+xlim(0, 200)+ylim(-0.5, 4.25)+
  geom_vline(xintercept=c(30,40), color = c("red", "red"))+
  geom_vline(xintercept= 160, color = "green")+
  geom_line(color="blue",size=2)+theme_void()
p
ggsave(plot = p, "../CD20 paper figures/Fig 3A (Ca imaging)/Ms4a4b_RS.pdf", width = 5, height = 5)


