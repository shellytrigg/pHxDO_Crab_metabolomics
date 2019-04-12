
#load libraries
library("reshape2", lib.loc="~/Library/R/3.3/library")
library("stringr", lib.loc="~/Library/R/3.3/library")
library("ggplot2", lib.loc="~/Library/R/3.3/library")
library("plyr", lib.loc="~/Library/R/3.3/library")
library("psych", lib.loc="~/Library/R/3.3/library")
library("tidyr", lib.loc="~/Library/R/3.3/library")
library("lme4", lib.loc="~/Library/R/3.3/library")
library("tidyr", lib.loc="~/Library/R/3.3/library")
library("survminer", lib.loc="~/Library/R/3.3/library")
library("survival", lib.loc="~/Library/R/3.3/library")
library("lmerTest", lib.loc="~/Library/R/3.3/library")

#read in stage data file
setwd("/Users/paul.mcelhany/Documents/Ocean Acidification/Experiments/Crab/Crab megalopae 2016/Data")
d <- read.csv("Updated Megalopae-Juvenile Data Sheet.csv", 
              stringsAsFactors=FALSE, skip = 1, header = TRUE)
str(d)
View(d)
#last day of observation
lastDayObs <- "9.10.16"
#remove M13
d <- subset(d, MOATS != "M13")

#read MOATS treatment table
dTreat <- read.csv("Treatments_2016_08_15.csv")
View(dTreat)

##### clean up input data
#remove blank rows
d <- subset(d, MOATS != "")
#convert to long skinny data
d <- gather(d, date, status, X6.8.16:X9.10.16)
# add treatments for each MOATS
d <- merge(d, dTreat, "MOATS")
#creat moats_jar ID
d$MOATSjar <- paste(d$MOATS, "_", d$JAR, sep = "")
#format date
d$dateString <- substring(as.character(d$date),2)
d$date <- as.POSIXlt(strptime(d$dateString, "%m.%d.%y"))
#remove respirometry crabs, not Dungness, and empty jars
d <-subset(d, substring(status,1,1) != "R" & status != "Not Dung." & status != "" & status != "---")
length(d$status)
levels(factor(d$status))

#create seperate rows if crabs removed and added on the same day
for(i in 1:length(d$status)){
  if(d$status[i] == "FREEZE/M1"){
    d$status[i] <- "FREEZE"
    newRow <- d[i,]
    newRow$status <- "M1"
    d <- rbind(d, newRow)
  }
  if(d$status[i] == "dead/M1"){
    d$status[i] <- "dead"
    newRow <- d[i,]
    newRow$status <- "M1"
    d <- rbind(d, newRow)
  }
}
d$status <- factor(d$status)
length(d$status)
levels(factor(d$status))

#sort by MOATSjar, date, and status (stats sort required for cases where crabs removed and started on same day)
d <- d[order(d$MOATSjar ,d$date, d$status), ]

#assign crabID 
d$crabID <- NA
for(i in 1:length(d$status)){
  if(d$status[i] == "M1"){
    cID <-  paste(d$MOATSjar[i], d$dateString[i], sep = "_")
  }
  d$crabID[i] <- cID
  if(d$status[i] == "dead" || d$status[i] == "FREEZE" || d$status[i] == "1J-dead"){
    cID <- NA
  }
}
#View(subset(d, is.na(crabID)))

#format as factors
d$treatment <- paste(d$pH_treatment, "_", d$DO_Treatment, sep = "")
d$MOATSjar <- as.factor(d$MOATSjar)
d$status <- as.factor(d$status)
as.data.frame(table(d[ , c("status")]))
d$crabID <- factor(d$crabID)
levels(d$status)

#make event columns
d$event <- ""

#assign events (assumes sorted by crabID, date and status)
for(i in 2:length(d$MOATS)){
  if(d$status[i-1] == "M1"){
    d$event[i-1] <- "Start"
  } 
  if(d$status[i-1] == "M1" && d$status[i] == "dead"){
    d$event[i] <- "mToDead"
  }
  if(d$status[i-1] == "M1" && d$status[i] == "J1-dead"){
    d$event[i] <- "mToJ1Dead"
  }
  if(d$status[i-1] == "m" && d$status[i] == "dead"){
    d$event[i] <- "mToDead"
  }
  if(d$status[i-1] == "m" && d$status[i] == "1J-dead"){
    d$event[i] <- "mToJ1Dead"
  }
  if(d$status[i-1] == "M1" && d$status[i] == "1J"){
    d$event[i] <- "mToJ1"
  }
  if(d$status[i-1] == "m" && d$status[i] == "1J"){
    d$event[i] <- "mToJ1"
  }
  if(d$status[i-1] == "1J" && d$status[i] == "dead"){
    d$event[i] <- "J1ToDead"
  }
  if(d$status[i-1] == "1J" && d$status[i] == "2J"){
    d$event[i] <- "J1ToJ2"
  }
  if(d$status[i-1] == "2J" && d$status[i] == "dead"){
    d$event[i] <- "J2ToDead"
  }
  if(d$status[i-1] == "2J" && d$status[i] == "FREEZE"){
    d$event[i] <- "J2ToFreeze"
  }
  if(d$event[i] == "" && d$dateString[i] == lastDayObs)
    d$event[i] <- "aliveNoMoltLastObs"
}
View(d)
#frequency of each event
as.data.frame(table(d[ , c("event")]))

#assign durations
d$duration <- NA
durCounter <- 0
for(i in 1:length(d$MOATS)){
  if(d$event[i] == "Start"){
    durCounter <- 0
  }
  d$duration[i] <- durCounter
  durCounter <- durCounter + 1
}


#########
#plot start dates
startDays <- levels(factor(as.character(subset(d, event == "Start")$date)))
startDays <- strptime(startDays, "%Y-%m-%d")
startDays <- as.POSIXct(startDays)
class(startDays)
aes(fill = factor(treatment))
levels(factor(subset(d, event == "Start")$treatment))
ggplot(subset(d, event == "Start"), aes(x=date)) + geom_histogram(aes(fill = treatment)) +
  theme_bw(base_size = 16) + 
  scale_x_datetime(breaks = startDays, date_labels = "%b %d") +
  xlab("Date in 2016") + ylab("Number of Megalopae Starting Experiment") + 
  theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=16))


#########
#Plot of time to transition and death by treatment
d$event2 <- ""
d$event2[d$event %in% c("J1ToDead", "J2ToDead", "mToDead", "mToJ1Dead")] <- "dead"
d$event2[d$event == "mToJ1"] <- "mToJ1"
d$event2[d$event == "J1ToJ2"] <- "J1ToJ2"

#Anova on time to J1ToJ2: no random effects
fit <- aov( duration~ pH_treatment*DO_Treatment, data = subset(d, event2 == "J1ToJ2"))
summary(fit)

#Anova on time to J1ToJ2: with random MOATS effects
fitME <- lmer(duration~ pH_treatment*DO_Treatment + (1|MOATS), data = subset(d, event2 == "J1ToJ2"))
summary(fitME)

#box plot
ggplot(subset(d, event2 != ""), aes(treatment, duration)) + 
  geom_boxplot(data = subset(d, event2 == "mToJ1")) +
  geom_boxplot(data = subset(d, event2 == "J1ToJ2")) + 
  geom_jitter(aes(color = event2), height = 0, width =  0.25, alpha = 0.5) + 
  geom_hline(aes(yintercept = 0)) +
  theme_bw(base_size = 18) + 
  scale_y_continuous(breaks = seq(0,65,5)) + 
  scale_x_discrete(labels = c("High pH\nHigh DO", "High pH\nLow DO", "Low pH\nHigh DO", "Low pH\nLow DO")) + 
  scale_color_manual(values=c("red", "blue", "green"), 
                     name="Molting Events",
                     breaks=c("mToJ1", "J1ToJ2", "dead"),
                     labels=c("Megalops to Juv1", "Juv1 to Juv2", "Died")) +
  ylab("Days Exposure") + xlab("Treatment") + 
  coord_flip()


#Plot with standard errors
#Run summarySE function at bottom of this page first
(meanHighPH <- mean(subset(d,pH_treatment == "High" & event2 == "J1ToJ2")$duration, na.rm = TRUE))
(meanLowPH <- mean(subset(d,pH_treatment == "Low"& event2 == "J1ToJ2")$duration, na.rm = TRUE))
meanLowPH - meanHighPH

summaryForErrorBars <- summarySE(subset(d, event2 == "J1ToJ2"), measurevar="duration", groupvars=c("treatment"), na.rm = TRUE)
View(summaryForErrorBars)
summaryForErrorBars$treatment <- factor(summaryForErrorBars$treatment)
summaryForErrorBars$treatment <- factor(summaryForErrorBars$treatment, levels= rev(levels(summaryForErrorBars$treatment)))
levels(summaryForErrorBars$treatment)
d$treatment <- as.character(d$treatment)
ggplot(summaryForErrorBars, aes(treatment, duration)) + 
  geom_errorbar(aes(ymin = duration - se, ymax = duration + se))+
  geom_point(size = 5) +
  geom_jitter(data = subset(d, event2 != ""), aes(x = treatment, y = duration, color = event2), 
              height = 0, width =  0.25, alpha = 0.5) + 
  theme_bw(base_size = 18) + 
  scale_y_continuous(breaks = seq(0,65,by=5)) +
  scale_x_discrete(labels = c("Low pH\nLow DO", "Low pH\nHigh DO", "High pH\nLow DO", "High pH\nHigh DO", "High pH\nLow DO")) +
  scale_color_manual(values=c("red", "blue", "green"),
                     name="Molting Events",
                     breaks=c("mToJ1", "J1ToJ2", "dead"),
                     labels=c("Megalops to Juv1", "Juv1 to Juv2", "Died")) +
  ylab("Days Exposure") + xlab("Treatment") + 
  coord_flip()


#############
# duration of Juvenile_1 stage

dEvents <- subset(d, event != ""& MOATS != "M13")
dCrab <- dcast(dEvents, crabID + treatment + pH_treatment + DO_Treatment + MOATS~ event, value.var = "duration")
dCrab$J1dur <- dCrab$J1ToJ2 - dCrab$mToJ1
#View(dCrab)

#Anova on Juv1 stage duration: no random effects
fit <- aov(J1dur ~ pH_treatment*DO_Treatment, data = dCrab)
summary(fit)

#Anova on Juv1 stage duration: with random MOATS effects
fitME <- lmer(J1dur ~ pH_treatment*DO_Treatment + (1|MOATS), data = dCrab)
summary(fitME)

#boxplot
ggplot(dCrab, aes(treatment, J1dur)) + 
  geom_boxplot() +
  geom_jitter(colour = "purple", height = 0, width =  0.25, alpha = 0.5) + 
  theme_bw(base_size = 18) + 
  scale_y_continuous(breaks = seq(15,60,5)) + 
  scale_x_discrete(labels = c("High pH\nHigh DO", "High pH\nLow DO", "Low pH\nHigh DO", "Low pH\nLow DO")) + 
  ylab("Duration of Juvenile Stage 1 (Days)") + xlab("Treatment") + 
  coord_flip()

#Plot with standard errors
#Run summarySE function at bottom of this page first
summaryForErrorBars <- summarySE(dCrab, measurevar="J1dur", groupvars=c("treatment"), na.rm = TRUE)
subset(dCrab,pH_treatment == "High")$J1dur
meanHighPH <- mean(subset(dCrab,pH_treatment == "High")$J1dur, na.rm = TRUE)
meanHighPH
meanLowPH <- mean(subset(dCrab,pH_treatment == "Low")$J1dur, na.rm = TRUE)
meanLowPH - meanHighPH
View(summaryForErrorBars)
summaryForErrorBars$treatment <- factor(summaryForErrorBars$treatment)
summaryForErrorBars$treatment <- factor(summaryForErrorBars$treatment, levels= rev(levels(summaryForErrorBars$treatment)))
levels(summaryForErrorBars$treatment)
dCrab$treatment <- as.character(dCrab$treatment)
ggplot(summaryForErrorBars, aes(treatment, J1dur)) + 
  geom_errorbar(aes(ymin = J1dur - se, ymax = J1dur + se))+
  geom_point(size = 5) +
  geom_jitter(data = dCrab, aes(treatment, J1dur), colour = "purple", height = 0, width =  0.25, alpha = 0.5) + 
  #geom_hline(aes(yintercept= meanHighPH)) +
  #geom_hline(aes(yintercept= meanLowPH)) +
  theme_bw(base_size = 18) + 
  scale_y_continuous(breaks = seq(15,60,by=5), limits = c(15,45)) + 
  scale_x_discrete(labels = c("Low pH\nLow DO", "Low pH\nHigh DO", "High pH\nLow DO", "High pH\nHigh DO", "High pH\nLow DO")) + 
  ylab("Duration of Juvenile Stage 1 (Days)") + xlab("Treatment") + 
  coord_flip()


################
#Survival plot
d$isDead <- NA
d$isDead[d$event %in% c("J1ToDead", "J2ToDead","mToDead", "mToJ1Dead")] <- 1
d$isDead[d$event == "J2ToFreeze"] <- 0
dSub <- subset(d, MOATS != "M13" & !is.na(isDead) & duration <= 44)
surv <- Surv(time = dSub$duration, event = dSub$isDead, type = "right")
pSurv <- ggsurvplot(survfit(surv ~ treatment, dSub), risk.table = FALSE, pval = FALSE, conf.int = FALSE,
                  font.main = 16, font.x =  16, font.y = 16, font.tickslab = 16, font.legend = 16,
                  break.time.by = 3, legend = c(0.2, 0.2), legend.title = "pHxDO Treatment",
                  legend.labs = levels(dSub$treatment))
pSurv + labs(title = "pHxDO Treatment") + xlab("Time (Days)") + ylim(0.9, 1)








########################
# extra stuff



ggplot(subset(d, event != "" & event != "Start"), aes(x=duration)) + geom_histogram() + facet_wrap(~ event)
ggplot(subset(d, status == "FREEZE"), aes(treatment, duration)) + geom_boxplot()
ggplot(subset(d, event == "J1ToJ2"), aes(treatment, duration)) + geom_boxplot()
ggplot(subset(d, event != "" & event != "Start"), aes(treatment, duration)) + geom_boxplot()+ facet_wrap(~ event)
ggplot(d, aes(treatment, duration)) + geom_boxplot()


#frequency of each event
as.data.frame(table(d[ , c("event")]))
dNo13 <- subset(d, MOATS != "M13")
moatsCount <- table(dNo13$status, dNo13$MOATS)
write.csv(file = "moatsCount.csv", moatsCount)
treatCount <- table(dNo13$status, dNo13$treatment)
write.csv(file = "treatCount.csv", treatCount)

ggplot(subset(d, event == "Start"), aes(x=date)) + geom_bar()

dStatusByTreat <- as.data.frame(table(dNo13$status, dNo13$treatment))
View(dStatusByTreat)

fracDead <- subset(dStatusByTreat, Var1 == "dead")$Freq / subset(dStatusByTreat, Var1 == "M1")$Freq
fracAlive <- subset(dStatusByTreat, Var1 == "FREEZE")$Freq / subset(dStatusByTreat, Var1 == "M1")$Freq
View(fracAlive)
treat <- subset(dStatusByTreat, Var1 == "M1")$Var2
dFracDead <- data.frame(treat, fracDead)
View(dFracDead)



dEvents <- subset(d, event != "")
dCrab <- dcast(dEvents, crabID + treatment + pH_treatment + DO_Treatment + MOATS~ event, value.var = "duration")
dCrab$J1dur <- dCrab$J1ToJ2 - dCrab$mToJ1
View(dCrab)


ggplot(subset(dCrab, MOATS != "M13"), aes(treatment, J1dur)) + geom_boxplot()
ggplot(subset(dCrab, MOATS != "M13"), aes(pH_treatment, J1dur)) + geom_boxplot()
ggplot(subset(dCrab, MOATS != "M13"), aes(DO_Treatment, J1dur)) + geom_boxplot()

table(dCrab$MOATS, dCrab$J2ToFreeze)
View(d)
View(dCrab)

fit <- aov(J1dur ~ pH_treatment*DO_Treatment, data = dCrab)
summary(fit)

fitME <- lmer(J1dur ~ pH_treatment*DO_Treatment + (MOATS), data = dCrab)
summary(fitME)



## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
} 
  
  
  
  
  