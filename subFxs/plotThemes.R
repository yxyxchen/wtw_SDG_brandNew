library('ggplot2')
# plot theme 
saveTheme = theme(panel.background = element_rect(fill = "white",colour = "white"),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank()) + 
  theme(title =element_text(size = 20, face='bold'), 
        plot.title = element_text(hjust = 0.5)) + 
  theme(axis.title=element_text(size= 20), axis.text = element_text(size=20), axis.line= element_line(color="black", size = 1)) +
  theme(strip.text = element_text(face="bold", size=15)) + 
  theme(legend.text=element_text(size= 15))

# plot theme 
displayTheme = theme(panel.background = element_rect(fill = "white",colour = "white"),
                  panel.grid.minor = element_blank(),
                  panel.grid.major = element_blank()) +
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+ 
  theme(axis.title=element_text(size=15), title =element_text(size=15, face='bold'), 
        plot.title = element_text(hjust = 0.5))+
  theme(axis.text = element_text(size=15), axis.line= element_line(color="black", size = 0.5)) +
  theme(strip.text = element_text(face="bold", size=15)) + 
  theme(legend.text=element_text(size= 15))