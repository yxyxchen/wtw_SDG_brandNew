library('dplyr')
library("tidyr")
library("ggplot2")
library("ggpubr")
source("subFxs/plotThemes.R")

# output dir
dir.create('figures/MFplot')

# load experiment parameters
load("expParas.RData")

# plot WTW timecourses in two environments
MFResults = MFAnalysis(isTrct = F)
sumStats = MFResults[['sumStats']]
timeWTW_ = MFResults[['timeWTW_']]
nSub = nrow(sumStats)
data.frame(wtw = unlist(timeWTW_),
           time = rep(tGrid, nSub),
           condition = rep(sumStats$condition, each = length(tGrid))) %>%
  group_by(condition, time) %>%
  dplyr::summarise(mu = mean(wtw, na.rm = T), se = sd(wtw, na.rm = T) / sqrt(sum(!is.na(wtw))),
                   min = mu- se, max = mu + se) %>%
  ggplot(aes(time, mu, linetype = condition)) +
  geom_ribbon(aes(ymin=min, ymax=max), fill = 'pink') +
  geom_line(color = themeColor, size = 1) +
  xlab("Elapsed time (s)") + ylab("WTW (s)") + 
  myTheme
ggsave("figures/MFPlot/wtw_timecourse.eps", width = 5, height = 4) 

# plot average WTWs in two environments
MFResults = MFAnalysis(isTrct = T)
sumStats = MFResults[['sumStats']]
sumStats %>% ggplot(aes(condition, muWTW)) + geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center',
               color = themeColor, fill = "pink",
               binwidth = 1.5, stroke = 6) + 
  stat_compare_means(comparisons = list(c("HP", "LP")),
                     aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                     bracket.size = 1, size = 6, label.y = 22) +
  xlab("") + ylab("AUC (s)") + ylim(c(0, 25)) + 
  myTheme 
ggsave("figures/MFPlot/muWTW_comparison.eps", width = 6, height = 6)


