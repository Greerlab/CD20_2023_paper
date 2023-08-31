Calcium image
================

## Installation

``` r
#install.packages("openxlsx","ggplot2")
```

## Read input file and plot

``` r
suppressMessages(library(openxlsx))
suppressMessages(library(ggplot2))
File = "20211203_Ms4a2,6D.xlsx"
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
p
```

![](Calcium_imaging_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
