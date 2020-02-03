library('dplyr')
library("tidyr")
library("ggplot2")
library("ggpubr")
library("lme4")
source("subFxs/plotThemes.R")
source("MFAnalysis.R")

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
  xmin = rep(blockSec * (0 : 2), 2), xmax = rep(blockSec * (0 : 2), 2) + blockSec - max(tMaxs),
  condition = rep(c("HP", "LP"), each = nBlock)
)
greyData = data.frame(
    xmin = rep(blockSec * (1 : 3), 2) - max(tMaxs), xmax = rep(blockSec * (1 : 3), 2),
    condition = rep(c("HP", "LP"), each = nBlock)
  )
data.frame(wtw = unlist(timeWTW_),
           time = rep(tGrid, nSub),
           condition = rep(sumStats$condition, each = length(tGrid))) %>%
  group_by(condition, time) %>%
  dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))),
                   min = mu- se, max = mu + se) %>%
  ggplot(aes(time, mu, linetype = condition)) +
  geom_rect(data = yellowData,
            aes(xmin=xmin, xmax=xmax, ymin=0, ymax=16), fill = "#ffffcc",inherit.aes = F) +
  geom_rect(data = greyData,
            aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = 16), fill = "#d9d9d9", inherit.aes = F) +
  geom_ribbon(aes(ymin=min, ymax=max), fill = 'pink') +
  geom_line(color = themeColor, size = 1) +
  xlab("Task time (s)") + ylab("WTW (s)") + 
  myTheme + ylim(0, 16) + ggtitle(expName) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) +
  scale_x_continuous(breaks = 0:3 * 420, labels = 0:3 * 700) + 
  theme(legend.position = "none")
ggsave("figures/MFPlot/wtw_timecourse.eps", width = 5, height = 3) 
ggsave("figures/MFPlot/wtw_timecourse.png", width = 5, height = 3) 


# plot the truncated wtw timecourse
blockTime = rep(seq(0, blockSec-1, by = 1), nBlock)
select = blockTime < (blockSec - max(tMaxs))
df = data.frame(wtw = unlist(timeWTW_),
           time = rep(tGrid, nSub),
           condition = rep(sumStats$condition, each = length(tGrid))) %>%
  group_by(condition, time) %>%
  dplyr::summarise(mu = mean(wtw, na.rm = F), se = sd(wtw, na.rm = F) / sqrt(sum(!is.na(wtw))),
                   min = mu- se, max = mu + se) 
df$mu[!select] = NA; df$max[!select] = NA; df$min[!select] = NA
df %>% ggplot(aes(time, mu, linetype = condition))  +
  geom_ribbon(aes(ymin=min, ymax=max), fill = 'pink') +
  geom_line(color = themeColor, size = 1) +
  xlab("Task time (s)") + ylab("WTW (s)") + 
  myTheme + ylim(0, 16) + ggtitle(expName) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) +
  scale_x_continuous(breaks = 0:3 * 420, labels = 0:3 * 700) + 
  theme(legend.position = "none")
ggsave("figures/MFPlot/wtw_timecourse_truc.eps", width = 5, height = 3) 
ggsave("figures/MFPlot/wtw_timecourse_truc.png", width = 5, height = 3) 


# plot average WTWs in two environments
MFResults = MFAnalysis(isTrct = T)
# wilcox.test(sumStats[sumStats$condition == "HP", "muWTW"],
#             sumStats[sumStats$condition == "LP", "muWTW"],paired = F)
sumStats = MFResults[['sumStats']]
blockStats = MFResults[['blockStats']]
sumStats %>% ggplot(aes(condition, muWTW)) + geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center',
               color = themeColor, fill = "pink",
               binwidth = 1.5, stroke = 6) + 
  stat_compare_means(comparisons = list(c("HP", "LP")),
                     aes(label = ..p.signif..), label.x = 1.5, symnum.args= symnum.args,
                     bracket.size = 1, size = 6, label.y = 22) +
  xlab("") + ylab("AUC (s)") + 
  myTheme + ggtitle(expName) + 
  theme(plot.title = element_text(face = "bold", hjust = 0.5, color = themeColor)) + 
  scale_y_continuous(breaks = c(0, 12, 24), limits = c(0, 26))
ggsave("figures/MFPlot/muWTW_comparison.eps", width = 4, height = 4)
ggsave("figures/MFPlot/muWTW_comparison.png", width = 4, height = 4)


# mixed effect anova
data = blockStats
data$condition = ifelse()
#data$condition = ifelse(blockStats$condition == "HP", 0, 1)
fit = lm(muWTW ~  condition , blockStats)
summary(fit)

data %>% group_by(manipulation, condition) %>%
  summarise(AUC = mean(muWTW))

