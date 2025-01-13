
library(ggplot2)


pop<-read.csv("demography/Hcan[br]_GMRFSkiride_df.csv", sep="\t", header = T)

 g <- ggplot(pop, aes(x=Time))+
      geom_ribbon(aes(ymin=Lower, ymax=Upper), fill = "grey50")+            
      geom_line(aes(y=Median), size = 1)+
      #geom_line(aes(y=Lower), colour = "gray60", size = 1)+
      #geom_line(aes(y=Upper), colour = "gray60", size = 1)+
      geom_vline(xintercept = 0.0042, linetype = "dashed", colour = "gray75", size=.65)+
      geom_vline(xintercept = 0.0129, linetype = "dashed", colour = "gray75", size=.65)+
      coord_cartesian(xlim=c(0, 0.02))+
      scale_x_continuous(name = "Time (Kybp)", labels = c(0, 5, 10, 15, 20))+
      ylab("Effective population size")+
      theme_classic()+
      theme(text = element_text(family = "serif", size = 32),
            plot.background = element_rect(linetype = "blank", fill = NA),
            panel.background = element_rect(linetype = "blank", fill = NA),
            panel.grid = element_blank())

g

png(filename = "Ne.png",
    width = 20, height = 20, units = "cm", bg = "transparent", res = 600)
g
dev.off()
