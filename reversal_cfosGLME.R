## GLME on binned response data from Manning et al (2018) 
#prepare data as follows:
#timestamp data as excel spreadsheet (e.g. "SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx") with three columns, A) subject ID, B) Correct lever press timestamp and C) Incorrect lever press timestamp. where there are more values in B than C just have zeros
#Blind details as excel spreadsheet (e.g. "Reversal cFos cohort blind.xlsx") with two columns, A) subject ID ("Phsyical Tag") and B)Genotype
#cFos cell counts as csv file (e.g. "cfos_040618_finalmice.xlsx") with columns including subject ID and cFos measurments arranged/labelled by region of interest


library(readr)
library(lme4)
library(ggplot2)
#library(dplyr)
library(tidyr)
library(psych)
library(gdata)
library(R.matlab)
library(xtable)
library(Hmisc)
library(nnet)
library(reshape2)
library(ggbiplot)
library(corrplot)
library(lsmeans)
library(factoextra)
library(ggfortify)
library(readxl)
library(linbin)
library(MASS)
library(lattice)
library(fastICA)
library(plotly)
library(stargazer)

setwd("/Users/Lizzie/Documents/GitHub/SAPAP3_reversal/")
d <- read_excel("/Users/Lizzie/Documents/GitHub/SAPAP3_reversal/SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx")
g <- read_excel("/Users/Lizzie/Documents/GitHub/SAPAP3_reversal/Reversal cFos cohort blind.xlsx")
cfos <- read_csv("/Users/Lizzie/Documents/GitHub/SAPAP3_reversal/cFos_040618_finalmice.csv")


# read in the data
setwd("~/code/SAPAP3_reversal/")
d <- read_excel("~/code/SAPAP3_reversal/SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx")
View(d)


# read in genotype
g <- read_excel("~/code/SAPAP3_reversal/Reversal cFos cohort blind.xlsx")
g$id <- g$`Physical Tag`
g <- g[,2:3]


# read in cfos cell counts

cfos <- read_csv("~/code/SAPAP3_reversal/cFos_040618_finalmice.csv")
names(cfos)[names(cfos)=="ID"] <- "id"
names(cfos)[names(cfos)=="correct"] <- "tot_correct"
names(cfos)[names(cfos)=="incorrrect"] <- "tot_incorrect"
names(cfos)[names(cfos)=="Genotype"] <- "genotype01"
names(cfos)[names(cfos)=="NAc S"] <- "NAcS"
names(cfos)[names(cfos)=="NAc C"] <- "NAcC"

View(cfos)

# check PCA on cfos

just_rois <- cfos[,7:16]
#head(just_rois)
feedROIcor <- cor(just_rois)
pdf("cfos correlations by regions.pdf", width=12, height=12)
corrplot(feedROIcor, cl.lim=c(0,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="orange", addCoefasPercent = FALSE,
         p.mat = 1-feedROIcor, sig.level=0.5, insig = "blank")
dev.off()

# look
summary(cfos.pca)


# make tall cfos
cshort <- cfos[,c(1:2,7:16)]
cshort <- as.data.frame(cshort)
t_cfos <- melt(cshort,id.vars = c("id","genotype01"))
names(t_cfos)[names(t_cfos)=="value"] <- "cfos"
names(t_cfos)[names(t_cfos)=="variable"] <- "region"

#test and plot cFos density by genoytpe across regions
summary(cm1 <- lm(cfos ~ region*genotype01, data = t_cfos))
anova(cm1)

pdf("cfos by region and genotype.pdf", width=12, height=6)
boxplot(cfos ~ genotype01 + region, data = t_cfos, main = "cfos activity by SAPAP3 genotype (blue = WT, purple = KO) and region",
        xlab = "Region", ylab = "cfos", varwidth = TRUE, col =  cm.colors(2))
dev.off()

# generate binned data. binsize is in 0.1 sec
d$id <- d$`subject ID`
d$corr <- d$`Correct lever press`
d$inc <- d$`Incorrect lever press`
d <- d[,4:6]
d[d==0] <- NA
View(d)
ids <- unique(d$id)
describe(d$corr)
describe(d$inc)
start <- 0
binsize <- 40
maxbin <- 16800

#censor data for mice that completed training in <30 minutes
c <- data.frame()
for (id in ids)
## the if condition that sets a different end for some ids, or just set end to equal the time of last response
  {end <- max(c(max(d$inc[d$id==id]),max(na.omit(d$corr[d$id==id]))))
  # make sure we are not removing idle time at the end of sessions for mice with full session durations
  if (end>maxbin) {end <- 18000}
  bins <- seq(start,end,binsize)
  t=cut(d$corr[d$id==id],c(bins), labels = FALSE)
b <- table( factor(t, levels = 1:length(bins)-1))
b <- as.data.frame(b)
b$id <- id
names(b)[names(b)=="Var1"] <- "time"
names(b)[names(b)=="Freq"] <- "corr"
tinc=cut(d$inc[d$id==id],c(bins), labels = FALSE)
i <- table( factor(tinc, levels = 1:length(bins)-1))
i <- as.data.frame(i)
b$inc <- i$Freq
c <- rbind(c,b)
}

c <- c[c$time!="0",]

c$t.ord <- as.ordered(c$time)
c$t.num <- as.numeric(c$time)

# plot(c$time,c$inc)
# plot(c$time,c$corr)
# 
# hist(c$inc)
# hist(c$corr)
length(c$inc)

#plot responses over time
ggplot(aes(x=t.num, y=corr),data = c) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Correct responses") +
  stat_smooth(data=c, aes(x=t.num, y=corr), size=0.4, alpha=0.1, se=FALSE, method="loess")

ggplot( aes(x=t.num, y=inc),data = c) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Incorrect responses") +
  stat_smooth(data=c, aes(x=t.num, y=inc), size=0.4, alpha=0.1, se=FALSE, method="loess")

# set up negative binomial models
# estimate overdispersion parameter theta
theta.corr <- theta.ml(c$corr, mean(c$corr), 513, limit = 50, eps = .Machine$double.eps^.25, trace = FALSE)
theta.inc <- theta.ml(c$inc, mean(c$inc), 513, limit = 50, eps = .Machine$double.eps^.25, trace = FALSE)



# for correct presses
#summary(mcorr1 <- glm(corr ~ t.ord + (1:id), family = negative.binomial(theta = theta.corr), data = c))
summary(mcorr2 <- glm(corr ~ t.num + (1:id), family = negative.binomial(theta = theta.corr), data = c))

#anova(mcorr1,mcorr2)

# for incorrect presses
summary(minc1 <- glm(inc ~ t.num + (1:id), family = negative.binomial(theta = theta.corr), data = c))

# combine all responses into one variable with correct/incorrect as factor
c.all <- melt(c,id.vars = c("id","time","t.num","t.ord"),value.name = "response", variable.name = "type")
names(c.all)[names(c.all)=="value"] <- "response"
names(c.all)[names(c.all)=="variable"] <- "type"

View(c.all)

nbfit <- suppressWarnings(fitdistr(c.all$response, "negative binomial"))
nbfit$estimate
logLik(nbfit)

c.all$response <- as.numeric(c.all$response)

#Plot and save response (log) curves
pdf("Response curves gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = type),data = c.all) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=c.all, aes(x=t.num, y=log(response + .1), color = type), size=2, alpha=0.2, se=TRUE, method="loess") + geom_count(position = "jitter")
dev.off()
ggplot(aes(x=t.num, y=response, color = type),data = c.all) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) +
  theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses") +
  geom_count(position = "jitter")

# merge with genotype
# corr/inc separated
hdf <- merge(c,g) # horizontal data frame
#View(hdf)

# all together
bdf <- merge(c.all,g)
#View(bdf)

names(cfos)[names(cfos)=="grooming time"] <- "grooming.time"

# merge with cell counts
hdfc <- merge(hdf,cfos)
bdfc <- merge(bdf,cfos)


pdf("Log response curves by genotype gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = type),data = bdf) + theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Responses (log)") +
  stat_smooth(data=bdf, aes(x=t.num, y=log(response + .1), color = type), size=2, alpha=0.2, se=TRUE, method="loess") + geom_count(position = "jitter") +
  facet_grid(Genotype ~ .)
dev.off()

pdf("Linear response curves by genotype gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=response, color = type),data = bdf) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) +
  theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Responses") +
  geom_count(position = "jitter") +
  facet_grid(Genotype ~ .)
dev.off()


#  quick look at genotype
# correct
summary(mcorrg1 <- glm(corr ~ t.num*Genotype +  (1:id), family = negative.binomial(theta = theta.corr), data = hdf))
#incorrect
summary(mincg1 <- glm(inc ~ t.num*Genotype +  (1:id), family = negative.binomial(theta = theta.inc), data = hdf))

summary(lmc1 <- lm(cfos$tot_correct ~ cfos$genotype01))
summary(lmi1 <- lm(cfos$tot_incorrect ~ cfos$genotype01))
summary(lmg1 <- lm(log(cfos$grooming.time) ~ cfos$genotype01))

# calculate theta for both responses
theta.resp <- theta.ml(na.omit(c.all$response), mean(na.omit(c.all$response)), 972, limit = 50, eps = .Machine$double.eps^.25, trace = FALSE)
summary(mrespg1 <- glm(response ~ 1 + t.num*type + type*Genotype +  (1:id), family = negative.binomial(theta = theta.resp), data = bdf))
car::Anova(mrespg1)

# plot over time by genotype
#lsmip(mrespg1, Genotype ~ 1 + t.num | type, at = list(t.num = c(bins[2],bins[length(bins)]/2, bins[length(bins)])), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )
ls.respg1 <- lsmeans(mrespg1,"type", by = "Genotype")
plot(ls.respg1, type ~ response, horiz=F,ylab = "Response levels", xlab = "type")


#test cFos ROI (alone) model e.g. PrL
summary(mrespg2 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + PrL*t.num*Genotype + PrL*type*Genotype + PrL*type*t.num +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg2)

#quick graphs to summarize interactions
#for "list" use 3 values of cFos measurement ~mean, mean + 1 SD, mean - 1 SD

lsmip(mrespg2, PrL ~ type | Genotype, at = list(PrL = c(10,100,200)), ylab = "log(response rate)", xlab = "type ", type = "predicted" )
lsmip(mrespg2, PrL ~ Genotype | type, at = list(PrL = c(10,100,200)), ylab = "log(response rate)", xlab = "type ", type = "predicted" )
lsmip(mrespg2, PrL ~ t.num | type, at = list(PrL = c(10,100,200),t.num  = c(1,length(bins)/2, length(bins))), ylab = "log(response rate)", xlab = "time ", type = "predicted" )




