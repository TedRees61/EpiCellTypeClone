#         read results summaries
#
res1<-read.csv("CellTypingResults1.csv")
res2<-read.csv("CellTypingResults2.csv")
res3<-read.csv("CellTypingResults3.csv")
#
#         simple plot - replace with violin
#
p1<-boxplot(res1$auPR)
p2<-boxplot(res2$auPR)
p3<-boxplot(res3$auPR)
p1+p2+p3
#
#       Simple scatter plot for overall results
#
#

RPR<-read.csv("ResultperRef.csv")
plot(x=log(as.number(RPR$Ref.Cells)),y=RPR$auPR)
x=RPR$Ref.Cells
z=log(x)
y=RPR$auPR
plot(x,y)
plot(z,y)
RPR$GT<-RPR$TP+RPR$FN
RPR$GT
#
#     Apply cell count filter > 8
#
sRPR<-RPR[RPR$GT>8,]
x=sRPR$Ref.Cells
z=log(x)
y=sRPR$auPR
plot(z,y)

#
#   GG plot section - use these plots
#

library(ggplot2)

df <- data.frame(
  x = RPR$auPR,
  y = log(RPR$TP+RPR$FN),
  group = RPR$Ref,each=10
)
ggplot(df, aes(x = x, y = y, color = group)) + 
  geom_violin()
  
df2<-data.frame(
    x=RPR$Ref,
    y=RPR$auPR
)
ggplot(df2, aes(x = x, y = y)) + 
  geom_violin()



comp_mathys<-data.frame(x=log(RPR$Ref.Cells),y=RPR$auPR,CellType=RPR$LegendCellList)
#
#   basic plot with smoothing line (delete later)
#
plot_mathys <-   ggplot(data=comp_mathys,aes(x = x, y = y),colour=RPR$LegendCellList) +
    geom_smooth(method = "lm", se = T,formula='y ~ x',              
                  color = 'grey20', alpha = 0.4,fill="white") +
  geom_point(size=3) +  annotate('text', x = -Inf, y = Inf, hjust = -1, vjust = 1.2,           
            label="Correlation of Precision Recall to Reference Cell Count") + 
  labs(y= "Area under Precision/Recall Curve", x = "Log (Ref Cell Counts)",colour="blue") +
  ggtitle("Mathys et al. -\npseudoreplication")+  theme_cowplot()+  theme(axis.text = element_text(size=9),        
            plot.title = element_text(size = 13, face = "plain"),        
            axis.title = element_text(size=11),        
            legend.text = element_text(size=11))
+  
  scale_colour_manual(values=c(0,1,2,3,4,5,6,7,8,9)) #just use a colour palette you like and replace pal here 
plot_mathys
#install.packages('devtools')
# devtools::install_github('smin95/smplot2', force = TRUE) # requires VPN if you are in China
#
#   Nice plot - use this one
#

library(cowplot)
library(smplot2)    # not sure we need this now
correl<-cor.test(comp_mathys$x,comp_mathys$y)
Rcorrel<-paste("Correlation of AUPR to Reference Cell Count\nP:",round(as.numeric(correl[3]),3)," Pearson:",round(as.numeric(correl[4]),3))
plot_mathys <-   ggplot(data=comp_mathys,aes(x = x, y = y),colour=RPR$LegendCellList)
plot_mathys <-   ggplot(data=comp_mathys,aes(x=x,y=y)) + 
  geom_point(size=2,mapping=aes(x=x,y=y,color=CellType)) + 
  geom_smooth(method = "lm", se = T,formula='y ~ x',              
              color = 'grey20', alpha = 0.4,fill="white") +
  annotate('text', x = -Inf, y = Inf, hjust = 0, vjust = 1, size=5,          
           label=Rcorrel) + 
  labs(y= "Area under Precision/Recall Curve", x = "Log (Ref Cell Counts)") 
plot_mathys
#
#   New plot - tidy up and get rid of Grey
#

correl<-cor.test(comp_mathys$x,comp_mathys$y)
Rcorrel<-paste(" P:",round(as.numeric(correl[3]),3)," Pearson:",round(as.numeric(correl[4]),3))
plot_mathys <-   ggplot(data=comp_mathys,aes(x = x, y = y),colour=RPR$LegendCellList)
plot_mathys <-   ggplot(data=comp_mathys,aes(x=x,y=y)) + 
  geom_point(size=2,mapping=aes(x=x,y=y,color=CellType)) + 
  geom_smooth(method = "lm", se = T,formula='y ~ x',              
              color = 'grey20', alpha = 0.4,fill="white") +
  annotate('text', x = -Inf, y = Inf, hjust = 0, vjust = 1, size=5,          
           label=Rcorrel) + 
  guides(color=guide_legend(title="Cell Type")) +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.3, vjust = 2, size = 20)) +
  labs(title="Correlation of AUPR to Reference Cell Count",y= "Area under Precision/Recall Curve", x = "Log (Ref Cell Counts)") 
plot_mathys

#
#   Second curve showing average prediction score
#
nRPR<-res1
nRPR$ref<-1
Nres2<-res2
Nres2$ref<-2
nRPR<-rbind(nRPR,Nres2)
Nres3<-res3
Nres3$ref<-3
nRPR<-rbind(nRPR,Nres3)
nRPR
BinRef1<-nRPR$ref==1
nRPR$ref_name[BinRef1] <- "1. Cerebellar Vermis"
BinRef2<-nRPR$ref==2
nRPR$ref_name[BinRef2] <- "2. Cerebellar Deep Nuclei"
BinRef3<-nRPR$ref==3
nRPR$ref_name[BinRef3] <- "3. Lateral Hemisphere"
comp_mathys<-NULL
comp_mathys<-data.frame(x=nRPR$MeanPredScore,
                        y=nRPR$auPR,
                        CellType=nRPR$LegendCellList,
                        Ref=as.character(nRPR$ref),
                        RefName=nRPR$ref_name)
correl<-cor.test(comp_mathys$x,comp_mathys$y)
Rcorrel<-paste("P:",round(as.numeric(correl[3]),3)," Pearson:",round(as.numeric(correl[4]),3))
#plot_mathys <-   ggplot(data=comp_mathys,aes(x = x, y = y),colour=nRPR$LegendCellList)
plot_mathys <-   ggplot(data=comp_mathys,aes(x=x,y=y)) + 
  geom_point(size=2,mapping=aes(x=x,y=y,color=CellType)) + 
  geom_smooth(method = "lm", se = T,formula='y ~ x',              
              color = 'grey20', alpha = 0.4,fill="white") +
  annotate('text', x = -Inf, y = Inf, hjust = -0.1, vjust = 1, size=5, label=Rcorrel) + 
  labs(title=" Correlation of AUPR to Mean Prediction Score", 
       y= "Area under Precision/Recall Curve", 
       x = "Mean Prediction Score") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, size = 20)) +
  guides(color=guide_legend(title="Cell Type"))
plot_mathys
#
# same but coloured by reference
plot_mathys <-   ggplot(data=comp_mathys,aes(x=x,y=y)) + 
  geom_point(size=2,mapping=aes(x=x,y=y,color=RefName)) + 
  geom_smooth(method = "lm", se = T,formula='y ~ x',              
              color = 'grey20', alpha = 0.4,fill="white") +
#              color = comp_mathys$RefName, alpha = 0.4,fill="white") +
  annotate('text', x = -Inf, y = Inf, hjust = -0.1, vjust = 1, size=5, label=Rcorrel) + 
  labs(title=" Correlation of AUPR to Mean Prediction Score", 
       y= "Area under Precision/Recall Curve", 
       x = "Mean Prediction Score") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, size = 20)) +
  guides(color=guide_legend(title="Reference"))
plot_mathys
# Nice violin plot by reference type
#
nRPR$ref_name[BinRef1] <- "1. Cerebellar Vermis"
nRPR$ref_name[BinRef2] <- "2. Cerebellar Deep Nuclei"
nRPR$ref_name[BinRef3] <- "3. Cerebellum Lat. Hemisphere"
comp_mathys<-data.frame(x=nRPR$MeanPredScore,
                        y=nRPR$auPR,
                        CellType=nRPR$LegendCellList,
                        Ref=as.character(nRPR$ref),
                        RefName=nRPR$ref_name)
ggplot(comp_mathys,aes(y=y,x=RefName))+
  geom_violin()+geom_boxplot(width=0.15)+  geom_point() +
  labs(title=" Violin Plot per Reference", 
       y= "Area under Precision/Recall Curve", 
       x = "Reference Brain Region") 


# Combined plot with AUPR
combined_plot <- ggplot(data = comp_mathys, aes(x = x, y = y, color = CellType, shape = RefName)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, formula = 'y ~ x', color = 'grey20', alpha = 0.4, fill = "white", aes(group = 1)) + # Added aes(group = 1) for single line
  annotate('text', x = -Inf, y = Inf, hjust = -0.1, vjust = 1, size = 5, label = Rcorrel) +
  labs(title = "Correlation of AUPR to Mean Prediction Score", 
       y = "Area under Precision/Recall Curve", 
       x = "Mean Prediction Score") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, size = 20)) +
  guides(color = guide_legend(title = "Cell Type"), shape = guide_legend(title = "Reference"))

# Show the combined plot
combined_plot
#
# plot data for correlation with mcc
#
comp_mathys<-data.frame(x=nRPR$MeanPredScore,
                        y=nRPR$MCC,
                        CellType=nRPR$LegendCellList,
                        Ref=as.character(nRPR$ref),
                        RefName=nRPR$ref_name)
correl<-cor.test(comp_mathys$x,comp_mathys$y)
Rcorrel<-paste("P:",round(as.numeric(correl[3]),5)," Pearson:",round(as.numeric(correl[4]),3))
#
# Violin Plot not used
#

ggplot(comp_mathys,aes(y=y,x=RefName))+
  geom_violin()+geom_boxplot(width=0.15)+  geom_point() +
  labs(title=" Violin Plot per Reference", 
       y= "Matthews Correlation Coefficient", 
       x = "Reference Brain Region") 

#
# Combined plot prediction score vs MCC
#
combined_plot <- ggplot(data = comp_mathys, aes(x = x, y = y, color = CellType, shape = RefName)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, formula = 'y ~ x', color = 'grey20', alpha = 0.4, fill = "white", aes(group = 1)) + # Added aes(group = 1) for single line
  annotate('text', x = -Inf, y = Inf, hjust = -0.1, vjust = 1, size = 5, label = Rcorrel) +
  labs(title = "Correlation of MCC to Mean Prediction Score", 
       y = "Matthews Correlation Coefficient", 
       x = "Mean Prediction Score") +
  theme_cowplot() +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2, size = 20)) +
  guides(color = guide_legend(title = "Cell Type"), shape = guide_legend(title = "Reference"))

# Show the combined plot
combined_plot

