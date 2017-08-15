## LME on binned response data from Lizzie Manning and Suzanne Ahmari's free operant deterministic reversal experiment

library(readr)
library(lme4)
library(ggplot2)
library(dplyr)
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


# read in the data
setwd("~/code/SAPAP3_reversal/")
d <- read_excel("~/code/SAPAP3_reversal/SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx")
#View(d)


# read in genotype
g <- read_excel("~/code/SAPAP3_reversal/Reversal cFos cohort blind.xlsx")
g$id <- g$`Physical Tag`
g <- g[,2:3]


# read in cfos cell counts
cfos <- read_csv("~/code/SAPAP3_reversal/SAPAP3 reversal cFos density_Final mice.csv")
names(cfos)[names(cfos)=="ID"] <- "id"
names(cfos)[names(cfos)=="correct"] <- "tot_correct"
names(cfos)[names(cfos)=="incorrrect"] <- "tot_incorrect"
names(cfos)[names(cfos)=="Genotype"] <- "genotype01"
#View(cfos)
# check PCA on cfos

just_rois <- cfos[,7:18]
#head(just_rois)
feedROIcor <- cor(just_rois)
pdf("cfos correlations by regions.pdf", width=12, height=12)
corrplot(feedROIcor, cl.lim=c(0,1),
         method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         order = "hclust", diag = FALSE,
         addCoef.col="orange", addCoefasPercent = FALSE,
         p.mat = 1-feedROIcor, sig.level=0.5, insig = "blank")
dev.off()

cfos.pca = prcomp((just_rois),scale = TRUE, center = TRUE)

# run PCA and write component scores
cfos_pcas <- get_pca_ind(cfos.pca)
# hdf$val1 <- cfos_pcas$coord[,1]
# hdf$val2 <- cfos_pcas$coord[,2]
# hdf$val3 <- cfos_pcas$coord[,3]
# hdf$val4 <- cfos_pcas$coord[,4]

# look
summary(cfos.pca)

# check the variance explained by various factors
plot(cfos.pca,type = 'l')


# only striatal
# str_rois <- cfos[,13:18]
# head(str_rois)
# feedROIcor <- cor(str_rois)
# corrplot(feedROIcor, cl.lim=c(0,1),
#          method = "circle", tl.cex = 1, tl.col = 'black',
#          order = "hclust", diag = FALSE,
#          addCoef.col="black", addCoefasPercent = FALSE,
#          p.mat = 1-feedROIcor, sig.level=0.6, insig = "blank")





d$id <- d$`subject ID`
d$corr <- d$`Correct lever press`
d$inc <- d$`Incorrect lever press`
d <- d[,4:6]
d[d==0] <- NA
View(d)
ids <- unique(d$id)
describe(d$corr)
describe(d$inc)
end <- 18000
start <- 0
binsize <- 1000
bins <- seq(start,end,binsize)

c <- data.frame()
for (id in ids)
{t=cut(d$corr[d$id==id],c(bins), labels = FALSE)
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

plot(c$time,c$inc)
plot(c$time,c$corr)

hist(c$inc)
hist(c$corr)
length(c$inc)

ggplot(aes(x=t.num, y=corr),data = c) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Time, 100s bins") + ylab("Correct responses") +
  stat_smooth(data=c, aes(x=t.num, y=corr), size=0.4, alpha=0.1, se=FALSE, method="loess")

ggplot( aes(x=t.num, y=inc),data = c) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Time, 100s bins") + ylab("Incorrect responses") +
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
pdf("Response curves gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = type),data = c.all) + theme_bw(base_size=20) + xlab("Time, 100s bins") + ylab("Responses (log)") +
  stat_smooth(data=c.all, aes(x=t.num, y=log(response + .1), color = type), size=2, alpha=0.2, se=TRUE, method="loess") + geom_count(position = "jitter")
dev.off()
ggplot(aes(x=t.num, y=response, color = type),data = c.all) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) +
  theme_bw(base_size=20) + xlab("Time, 100s bins") + ylab("Responses") +
  geom_count(position = "jitter")

# merge with genotype
# corr/inc separated
hdf <- merge(c,g) # horizontal data frame
#View(hdf)

ggplot(aes(x=t.num, y=corr, color = Genotype),data = hdf) + stat_smooth(size=2, alpha=0.1, method="loess", size = 1) + theme_bw(base_size=20) +
  xlab("Time, 100s bins") + ylab("Correct responses") +
  geom_count(position = "jitter")

ggplot( aes(x=t.num, y=inc, color = Genotype),data = hdf) + stat_smooth(size=2, alpha=0.1, method="loess",  size = 1) + theme_bw(base_size=20) +
  xlab("Time, 100s bins") + ylab("Incorrect responses") +
  geom_count(position = "jitter")


# all together
bdf <- merge(c.all,g)
#View(bdf)

names(cfos)[names(cfos)=="grooming time"] <- "grooming.time"

# merge with cell counts
hdfc <- merge(hdf,cfos)
bdfc <- merge(bdf,cfos)


pdf("Log esponse curves by genotype gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = type),data = bdf) + theme_bw(base_size=20) + xlab("Time, 100s bins") + ylab("Responses (log)") +
  stat_smooth(data=bdf, aes(x=t.num, y=log(response + .1), color = type), size=2, alpha=0.2, se=TRUE, method="loess") + geom_count(position = "jitter") +
  facet_grid(. ~ Genotype)
dev.off()

pdf("Linear esponse curves by genotype gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=response, color = type),data = bdf) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) +
  theme_bw(base_size=20) + xlab("Time, 100s bins") + ylab("Responses") +
  geom_count(position = "jitter") +
  facet_grid(. ~ Genotype)
dev.off()


#  quick look at genotype
# correct
summary(mcorrg1 <- glm(corr ~ t.num*Genotype + I(t.num^2) + I(t.num^3) + (1:id), family = negative.binomial(theta = theta.corr), data = hdf))
#incorrect
summary(mincg1 <- glm(inc ~ t.num*Genotype + I(t.num^2) + I(t.num^3) + (1:id), family = negative.binomial(theta = theta.inc), data = hdf))

summary(lmc1 <- lm(cfos$tot_correct ~ cfos$genotype01))
summary(lmi1 <- lm(cfos$tot_incorrect ~ cfos$genotype01))
summary(lmg1 <- lm(cfos$grooming.time ~ cfos$genotype01))

# both
# calculate theta for both responses
theta.resp <- theta.ml(c.all$response, mean(c.all$response), 1026, limit = 50, eps = .Machine$double.eps^.25, trace = FALSE)
summary(mrespg1 <- glm(response ~ 1 + t.num*type+ type*Genotype + t.num*Genotype +  (1:id), family = negative.binomial(theta = theta.resp), data = bdf))
car::Anova(mrespg1)

# plot over time by genotype
#lsmip(mrespg1, Genotype ~ 1 + t.num | type, at = list(t.num = c(1000, 9000, 18000)), ylab = "log(response rate)", xlab = "Time, s ", type = "predicted" )
ls.respg1 <- lsmeans(mrespg1,"type", by = "Genotype")
plot(ls.respg1, type ~ response, horiz=F,ylab = "Response levels", xlab = "type")



# grooming time -- I understand this to be the index of OCD-like behavior
summary(mrespg2 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + grooming.time*t.num*type +  (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg2)

summary(mrespg3 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + DLS*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg3)


summary(mrespg4 <- glm(response ~ t.num*type+ t.num*Genotype + type*Genotype + bdfc$`NAc S`*t.num*type*Genotype +   (1:id), family = negative.binomial(theta = theta.resp), data = bdfc))
car::Anova(mrespg4)

# conclusion: strong effects of cfos, weak effect of genotype, + interactions