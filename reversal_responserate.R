## LME on binned response data from Lizzie Manning and Suzanne Ahmari's free operant deterministic reversal experiment

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
library(rio)
# read in the data
setwd("/Users/Lizzie/Documents/GitHub/SAPAP3_reversal/")
d <- read_excel("/Users/Lizzie/Documents/GitHub/SAPAP3_reversal/SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx")
#setwd("~/code/SAPAP3_reversal/")
#d <- read_excel("~/code/SAPAP3_reversal/SAPAP3 cFos cohort lever press timestamp reversal day 1.xlsx")
View(d)


# read in genotype
g <- read_excel("/Users/Lizzie/Documents/GitHub/SAPAP3_reversal/Reversal cFos cohort blind.xlsx")
#g <- read_excel("~/code/SAPAP3_reversal/Reversal cFos cohort blind.xlsx")
g$id <- g$`Physical Tag`
g <- g[,2:3]


# read in cfos cell counts
#removed M2/M1 January 2018
#cfos <- read_csv("/Users/Lizzie/Documents/GitHub/SAPAP3_reversal/cFos_040618_finalmice.csv")
#cfos <- read_csv("~/code/SAPAP3_reversal/SAPAP3 reversal cFos density_Final mice.csv")
#names(cfos)[names(cfos)=="ID"] <- "id"
#names(cfos)[names(cfos)=="correct"] <- "tot_correct"
#names(cfos)[names(cfos)=="incorrrect"] <- "tot_incorrect"
#names(cfos)[names(cfos)=="Genotype"] <- "genotype01"
#names(cfos)[names(cfos)=="NAcS"] <- "NAcS"
#names(cfos)[names(cfos)=="NAcC"] <- "NAcC"


#View(cfos)
# check PCA on cfos

#just_rois <- cfos[,7:16]
#head(just_rois)
#feedROIcor <- cor(just_rois)
#pdf("cfos correlations by regions_no.pdf", width=12, height=12)
#corrplot(feedROIcor, cl.lim=c(0,1),
         #method = "circle", tl.cex = 1.5, type = "upper", tl.col = 'black',
         #order = "hclust", diag = FALSE,
         #addCoef.col="orange", addCoefasPercent = FALSE,
         #p.mat = 1-feedROIcor, sig.level=0.5, insig = "blank")
#dev.off()


#cfos.pca = prcomp((just_rois),scale = TRUE, center = TRUE)




# run PCA and write component scores
#cfos_pcas <- get_pca_ind(cfos.pca)
#cfos$val1 <- cfos_pcas$coord[,1]
# hdf$val2 <- cfos_pcas$coord[,2]
# hdf$val3 <- cfos_pcas$coord[,3]
# hdf$val4 <- cfos_pcas$coord[,4]

# look
#summary(cfos.pca)

# check the variance explained by various factors
#plot(cfos.pca,type = 'l')




c#fos$genotype01 <- as.factor(cfos$genotype01)

#ggplot(aes(x=genotype01, y=IL),data = cfos) +
  #geom_count()
#ggplot(aes(x=genotype01, y=DMS),data = cfos) +
  #geom_count()
#ggplot(aes(x=genotype01, y=NAcS),data = cfos) +
  g#eom_count()



# make tall cfos
#cshort <- cfos[,c(1:2,7:16)]
#cshort <- as.data.frame(cshort)
#t_cfos <- melt(cshort,id.vars = c("id","genotype01"))
#View(t_cfos)
#names(t_cfos)[names(t_cfos)=="value"] <- "cfos"
n#ames(t_cfos)[names(t_cfos)=="variable"] <- "region"

#summary(cm1 <- lm(cfos ~ region*genotype01, data = t_cfos))
#anova(cm1)

#pdf("cfos by region and genotype.pdf", width=12, height=6)
#boxplot(cfos ~ genotype01 + region, data = t_cfos, main = "cfos activity by SAPAP3 genotype (blue = WT, purple = KO) and region",
       # xlab = "Region", ylab = "cfos", varwidth = TRUE, col =  cm.colors(2))
#dev.off()


d$id <- d$`subject ID`
d$corr <- d$`Correct lever press`
d$inc <- d$`Incorrect lever press`
d <- d[,4:6]
d[d==0] <- NA
View(d)
ids <- unique(d$id)
describe(d$corr)
describe(d$inc)
# end <- 18000
start <- 0
# binsize <- 1000
# try smaller bins:
binsize <- 40
maxbin <- 16800


c <- data.frame()
for (id in ids)
## Lizzie, here you can insert the if condition that sets a different end for some ids, or just set end to equal the time of last response
  {end <- max(c(max(d$inc[d$id==id]),max(na.omit(d$corr[d$id==id]))))
  # make sure we are not removing idle time at the end of sessions for those 5-6 mice
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

plot(c$time,c$inc)
plot(c$time,c$corr)

hist(c$inc)
hist(c$corr)
length(c$inc)

ggplot(aes(x=t.num, y=corr),data = c) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Correct responses") +
  stat_smooth(data=c, aes(x=t.num, y=corr), size=0.4, alpha=0.1, se=FALSE, method="loess")

ggplot( aes(x=t.num, y=inc),data = c) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) + theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Incorrect responses") +
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

# AD: we should not need manual censoring

#remove bins with no data
# censored <- c(54:55,124:126,156:162,230:234,249:252,335:342,359:360,371:378,411:414,449:450,540:541,610:612,642:648,716:720,735:738,821:828,845:846,857:864,897:900,935:936)
# censoredtemp <- censored*1000/binsize;
# c.all[c(censoredtemp),c('response')] <- "NA"
# c.all[c(54,124:126,156:162,230:234,249:252,335:342,359:360,371:378,411:414,449:450,540,610:612,642:648,716:720,735:738,821:828,845:846,857:864,897:900,935,936),c('response')] <- "NA"
c.all$response <- as.numeric(c.all$response)

pdf("Response curves gam.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = type),data = c.all) + theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Responses (log)") +
  stat_smooth(data=c.all, aes(x=t.num, y=log(response + .1), color = type), size=2, alpha=0.2, se=TRUE, method="loess") + geom_count(position = "jitter")
dev.off()
ggplot(aes(x=t.num, y=response, color = type),data = c.all) + stat_smooth(size=2, alpha=0.2, method="gam", formula = y ~ s(x), size = 1) +
  theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Responses") +
  geom_count(position = "jitter")

# merge with genotype
# corr/inc separated
hdf <- merge(c,g) # horizontal data frame
#View(hdf)

ggplot(aes(x=t.num, y=corr, color = Genotype),data = hdf) + stat_smooth(size=2, alpha=0.1, method="loess", size = 1) + theme_bw(base_size=20) +
  xlab("Time, 200s bins") + ylab("Correct responses") +
  geom_count(position = "jitter")

ggplot( aes(x=t.num, y=inc, color = Genotype),data = hdf) + stat_smooth(size=2, alpha=0.1, method="loess",  size = 1) + theme_bw(base_size=20) +
  xlab("Time, 200s bins") + ylab("Incorrect responses") +
  geom_count(position = "jitter")


# all together
bdf <- merge(c.all,g)
#View(bdf)

names(cfos)[names(cfos)=="grooming time"] <- "grooming.time"

# merge with cell counts
hdfc <- merge(hdf,cfos) #might get angry here becuase I removed all reference to cFos data
bdfc <- merge(bdf,cfos) #might get angry here becuase I removed all reference to cFos data


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

pdf("Log response curves by genotype gam combine responses1.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1)),data = bdf) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=bdf, aes(x=t.num, y=log(response + .1)), size=1, alpha=0.2, se=TRUE, method="loess")
dev.off()

pdf("Log response curves by genotype gam combine responses1a.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1)),data = bdf) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=bdf, aes(x=t.num, y=log(response + .1)), size=1, alpha=0.2, se=TRUE, method="loess")+ facet_grid(Genotype ~ .)
dev.off()

pdf("Log response curves by genotype gam combine responses1B.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = Genotype),data = bdf) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=bdf, aes(x=t.num, y=log(response + .1), color = Genotype), size=1, alpha=0.2, se=TRUE, method="loess")
dev.off()

pdf("Log response curves by genotype gam combine responses1inc.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(inc + .1), color = Genotype),data = hdf) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=hdf, aes(x=t.num, y=log(inc + .1), color = Genotype), size=1, alpha=0.2, se=TRUE, method="loess")
dev.off()


pdf("Log response curves by genotype gam combine responses1corr.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(corr + .1), color = Genotype),data = hdf) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=hdf, aes(x=t.num, y=log(corr + .1), color = Genotype), size=1, alpha=0.2, se=TRUE, method="loess")
dev.off()

pdf("Log response curves by genotype gam combine responses1C.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1), color = Genotype + type),data = bdf) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=bdf, aes(x=t.num, y=log(response + .1), color = type), size=1, alpha=0.2, se=TRUE, method="loess")
dev.off()

pdf("Log response curves by genotype gam combine responses2.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1)),data = bdf) + theme_bw(base_size=20) + xlab("Time, 4s bins") + ylab("Responses (log)") +
  stat_smooth(data=bdf, aes(x=t.num, y=log(response + .1)), size=1, alpha=0.2, se=TRUE, method="loess") +
  facet_grid(Genotype ~ .)
dev.off()

pdf("Log response curves by genotype gam3.pdf", width=8, height=6)
ggplot(aes(x=t.num, y=log(response + .1)),data = bdf) + theme_bw(base_size=20) + xlab("Time, 200s bins") + ylab("Responses (log)") +
  stat_smooth(data=bdf, aes(x=t.num, y=log(response + .1)), size=2, alpha=0.2, se=TRUE, method="loess") + geom_count(position = "jitter") +
  facet_grid(Genotype ~ .)
dev.off()

