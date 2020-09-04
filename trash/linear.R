library("lme4")
library("tidyverse")

df <- read.csv("../data/dataset.csv", na.strings=".") %>% select(ID, TIME, DaT.scan, G2019S) %>% na.omit

plot(DaT.scan~TIME, data=df)

model <- lmer(DaT.scan~TIME+(TIME|G2019S), data=df)

summary(model)
plot(model)
abline(model)

model <- lmer(data=df, DaT.scan~TIME+(1|ID))
