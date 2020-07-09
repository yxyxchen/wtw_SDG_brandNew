library('dplyr')
library("tidyr")
library("ggplot2")
library("ggpubr")
library("lme4")
source("subFxs/plotThemes.R")
source("MFAnalysis.R")
library("lmerTest")
options(contrasts = c("contr.treatment", "contr.poly"))

# output dir
dir.create('figures/MFplot')

# load experiment parameters
load("expParas.RData")

# plot WTW timecourses in two environments
MFResults = MFAnalysis(isTrct = F)
sumStats = MFResults[['sumStats']]
timeWTW_ = MFResults[['timeWTW_']]
nSub = nrow(sumStats)
## use background color to distinguish used and excluded data 
yellowData = data.frame(
  xmin = rep(blockSec * (0 : 2), 2), xmax = rep(blockSec * (0 : 2), 2) + blockSec - max(delayMaxs),
  condition = rep(c("HP", "LP"), each = nBlock)
)
greyData = data.frame(
    xmin = rep(blockSec * (1 : 3), 2) - max(delayMaxs), xmax = rep(blockSec * (1 : 3), 2),
    condition = rep(c("HP", "LP"), each = nBlock)
  )
plotData = data.frame(wtw = unlist(timeWTW_),
                      time = rep(tGrid, nSub),
                      condition = rep(sumStats$condition, each = length(tGrid)))
isSig = lapply(1 : length(tGrid), function(i) {
  t = tGrid[i]
  HP = plotData$wtw[plotData$time == t & plotData$condition == "HP"]
  LP = plotData$wtw[plotData$time == t & plotData$condition == "LP"]
  tempt = wilcox.test(HP, LP)
  tempt$p.value
})
plotData %>%
  group_by(condition, time) %>%
  dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))),
                   min = mu- se, max = mu + se) %>%
  ggplot(aes(time, mu, color = condition)) +
  geom_rect(data = yellowData,
            aes(xmin=xmin, xmax=xmax, ymin=0, ymax=16), fill = "white",inherit.aes = F) +
  geom_rect(data = greyData,
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 16), fill = "#d9d9d9", inherit.aes = F) +
  geom_ribbon(aes(ymin=min, ymax=max, fill = condition, color = NA), alpha = 0.5) +
  geom_line(aes(color = condition), size = 1) +
  geom_point(data = data.frame(t = tGrid, isSig = ifelse(unlist(isSig) < 0.05, 16, NA)),
             aes(t, isSig), inherit.aes = F, shape = 4) +
  xlab("Task time (min)") + ylab("WTW (s)") + 
  myTheme + ylim(0, 16)  +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) +
  scale_x_continuous(breaks = 0:3 * 420, labels = 0:3 * 7) + 
  theme(legend.position = "none") +
  scale_fill_manual(values = conditionColors) +
  scale_color_manual(values = conditionColors)
ggsave("figures/MFPlot/wtw_timecourse.eps", width = 5, height = 3) 
ggsave("figures/MFPlot/wtw_timecourse.png", width = 5, height = 3) 


# plot average WTWs in two environments
MFResults = MFAnalysis(isTrct = T)
sumStats = MFResults[['sumStats']]
aovRes = aov(muWTW ~ condition, data = sumStats)
blockStats = MFResults[['blockStats']]
blockStats$blockNum = as.numeric(blockStats$blockNum)

sumStats %>% group_by(condition) %>% summarise(median(muWTW))
wilcox.test(sumStats$muWTW[sumStats$condition == "HP"],
            sumStats$muWTW[sumStats$condition == "LP"])
sumStats %>% ggplot(aes(condition, muWTW))  +
  geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) + 
  stat_compare_means(comparisons = list(c("HP", "LP")),
                     aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                     bracket.size = 1, size = 6, label.y = 22) + ylab("AUC (s)") + 
  myTheme  + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) + 
  scale_y_continuous(breaks = c(0, 12, 24), limits = c(0, 26)) +
  scale_fill_manual(values = conditionColors) +
  theme(legend.position = "none") + xlab("")
ggsave("figures/MFPlot/muWTW_comparison.eps", width = 4, height = 4)
ggsave("figures/MFPlot/muWTW_comparison.png", width = 4, height = 4)

# plot CIP
sumStats %>% group_by(condition) %>% summarise(median(stdWTW))
aovRes = aov(stdWTW ~ condition, data = sumStats)
wilcox.test(sumStats$stdWTW[sumStats$condition == "HP"],
            sumStats$stdWTW[sumStats$condition == "LP"])
library("latex2exp")
sumStats %>% ggplot(aes(condition, stdWTW))  +
  geom_dotplot(binaxis='y', stackdir='center', aes(fill = condition)) + 
  stat_compare_means(comparisons = list(c("HP", "LP")),
                     aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                     bracket.size = 1, size = 6, label.y = 10) + ylab(TeX("CIP ($s^2$)")) + 
  myTheme  + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) + 
  scale_y_continuous(breaks = c(0, 5, 10), limits = c(0, 11)) +
  scale_fill_manual(values = conditionColors) +
  theme(legend.position = "none") + xlab("")
ggsave("figures/MFPlot/stdWTW_comparison.eps", width = 4, height = 4)
ggsave("figures/MFPlot/stdWTW_comparison.png", width = 4, height = 4)

# correlation between 
cor.test(sumStats$muWTW, sumStats$stdWTW, method = "spearman")
cor.test(sumStats$muWTW[sumStats$condition == "HP"], sumStats$stdWTW[sumStats$condition == "HP"], method = "spearman")
cor.test(sumStats$muWTW[sumStats$condition == "LP"], sumStats$stdWTW[sumStats$condition == "LP"], method = "spearman")
sumStats %>% ggplot(aes(muWTW, stdWTW)) +
  geom_point(aes(color = condition), size = 3) +
  facet_grid(~condition) + myTheme +
  xlab("AUC (s)") + ylab(TeX("CIP ($s^2$)")) +
  scale_color_manual(values = conditionColors) +
  theme(legend.position = "none")
ggsave("figures/MFPlot/stdWTW_muWTW.eps", width = 8, height = 4)
ggsave("figures/MFPlot/stdWTW_muWTW.png", width = 8, height = 4)
  
  

# survival curve
## data for the optimal lines 
optim = data.frame(
  t = rep(kmGrid,  2),
  surv = rep(1, length(kmGrid) * 2),
  condition = rep(conditions, each = length(kmGrid)),
  select = rep(1:2, each = length(kmGrid))
) 
optim$surv[optim$condition == "LP" & kmGrid> optimWaitThresholds$LP] = 0 # quit after 2.2 s
optim$surv[optim$condition == "HP" & kmGrid> optimWaitThresholds$LP] = NA # don't plot after 2.2 s
optim$select[optim$condition == "HP" & kmGrid <= optimWaitThresholds$LP] = rep(1:2, each = 3) # plot interleaving colors 
optim$select[optim$condition == "LP" & kmGrid <= optimWaitThresholds$LP] = rep(1:2, each = 3) # plot interleaving colors
## stats test
survCurve_ = MFResults$survCurve_
plotData = data.frame(survCurve = unlist(survCurve_),
                      time = rep(kmGrid, nSub),
                      condition = rep(sumStats$condition, each = length(kmGrid)))
isSig = lapply(1 : length(kmGrid) , function(i)
{
  t = kmGrid[i]
  HP = plotData$survCurve[plotData$condition == "HP" & plotData$time == t]
  LP = plotData$survCurve[plotData$condition == "LP" & plotData$time == t]
  tempt = wilcox.test(HP, LP)
  tempt$p.value
}
  )
sigData = data.frame(
  t = kmGrid,
  isSig = ifelse(unlist(isSig) < 0.05, 1.02, NA)
)


## plot
plotData %>%
  group_by(condition, time) %>%
  dplyr::summarise(mu = mean(survCurve, na.rm = F), se = sd(survCurve, na.rm = F) / sqrt(sum(!is.na(survCurve))),
                   min = mu- se, max = mu + se) %>%
  ggplot(aes(time, mu, color = condition, fill = condition)) + geom_line() +
  geom_ribbon(aes(time, ymin = min, ymax = max), alpha = 0.5, color = NA) +
  geom_line(data = optim, aes(t, surv, color = condition, linetype = condition, alpha = condition), size = 1.2) +
  geom_line(data = data.frame(t = kmGrid[kmGrid > 2],surv = 1),
            aes(t, surv), color = conditionColors[1], size = 1.2, inherit.aes = F, alpha = 0.8) + 
  geom_point(data = sigData, aes(t, isSig), inherit.aes = F, color = "black", shape = 4, size = 0.8) + 
  scale_fill_manual(values = conditionColors) +
  scale_color_manual(values = conditionColors) +
  scale_linetype_manual(values = c("solid", "dotted")) +
  scale_alpha_manual(values = c(0.8, 1))+
  xlab("Elapsed time (s)") + ylab("Survival rate") + myTheme +
    theme(legend.position = "none") 
ggsave("figures/MFPlot/survival_curve.eps", width = 4, height = 4)
ggsave("figures/MFPlot/survival_curve.png", width = 4, height = 4) 

# mixed effect anova, I should be more interested in non-stressed. 
fit = lmer(muWTW ~  condition * blockNum + (1 | id), data)
summary(fit)
wilcox.test(blockStats[blockStats$manipulation == 1, "muWTW"],
            blockStats[blockStats$manipulation == 2, "muWTW"], paired = F)
