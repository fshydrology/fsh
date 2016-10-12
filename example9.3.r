mydata <- read.csv(url('https://github.com/fshydrology/fsh/raw/master/data9.3.csv'),row.names=1)

cor(mydata) # correlation matrix

step1 <- lm( Runoff ~ A, data=mydata)
summary(step1)
anova(step1)

step2 <- lm( Runoff ~ A + S, data=mydata)
summary(step2)
anova(step2)

step3 <- lm( Runoff ~ A + S + Pr, data=mydata)
summary(step3)
anova(step3)