library("logitnorm")
library("ggplot2")

ggplot(expPara[passCheck,], aes(phi)) + geom_histogram()
ggplot(expPara[passCheck,], aes(log(phi))) + geom_histogram()
ggplot(expPara[passCheck,], aes(tau)) + geom_histogram()
ggplot(expPara[passCheck,], aes(prior)) + geom_histogram()


mx = seq(0, 0.1, length.out = 1000)
median = logit(0.0052 + 0.1)
y = dlogitnorm(x, median, sigma = 0.8)
plot(x, y)
