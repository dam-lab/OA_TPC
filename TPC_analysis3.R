setwd("C:/users/james/Documents/Grad_school/OA_Project/Cost_experiments/TPC/")

library(data.table)

library(dplyr)

library(ggplot2)

library(tidyr)

library(broom)

library(lme4)
library(car)
library(pROC)
library(emmeans)
library(sjPlot)

TPC_data <- fread(file = "TPC_Data_complete.txt")

## count the numbers of points for each temp
TPC.mean <- TPC_data %>%
  group_by(temp) %>%
  summarise(mean = mean(survivorship, na.rm = TRUE),
            sd = sd(survivorship, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)

xlab <- "Temperature (°C)"

## create unique dummy variables based on treatment, dev.treatment, and generation
TPC_data <- unite(TPC_data,#the data frame
                  unique, #the name of the new column
                  c(treatment, dev.treatment, generation), #the existing columns you want to combine
                  remove = FALSE) 

TPC.list <- split(TPC_data, f = TPC_data$unique)


#TPC_data$environment <- as.character(NA)

#TPC_data[1:896,8] <- as.character("Home")
#TPC_data[897:nrow(TPC_data),8] <- as.character("Transplant")

### complete graph comparing against AA and HH





TPC_plot_complete <- ggplot(data = TPC_data, aes(temp, #the x-axis value
                                                 survivorship, #the y-axis value
                                                 color = unique))+ # change the lines to be different based on sex
  geom_point(position = position_jitter(height = 0.025, width = 0.5))+ # make a scatter plot with the points slightly offset from each other
  
  geom_smooth(method = "glm", # make the plotted model match the model we use to test
              method.args = list(family = "binomial"),# make sure we know that there is only 2 possible data values
              se = TRUE)+ # do not plot standard error
  
  scale_y_continuous(breaks = c(0,0.5,1))+ # 3 values on the y-axis displayed
  scale_colour_manual(name = "Experimental\nCondition", #need to use scale_colour_manual for line graphs
                      values = c("Blue",
                                 "Brown",
                                 "Black",
                                 "Gold",
                                 "Red"),
                      labels = c("Ambient",
                                 "Ambient in GW",
                                 "GW in Ambient F1",
                                 "GW in Ambient F9",
                                 "GW"))+
  theme_light()+
  xlab(xlab)+
  ylab("Survivorship")+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.title.align = 0.5)
TPC_plot_complete

TPC_data_2 <- filter(TPC_data, !grepl("HH_AA_9", unique))

TPC_plot <- ggplot(data = TPC_data_2, aes(temp, #the x-axis value
                                        survivorship, #the y-axis value
                                        color = treatment))+ # change the lines to be different based on sex
  geom_point(alpha = 0.5,
             size = 4,
             position = position_jitter(height = 0.025, width = 0.5))+ # make a scatter plot with the points slightly offset from each other
  
  geom_smooth(method = "glm", # make the plotted model match the model we use to test
              method.args = list(family = "binomial"),# make sure we know that there is only 2 possible data values
              se = FALSE,
              size = 2)+ # do not plot standard error
  
  scale_y_continuous(breaks = c(0,0.5,1))+ # 3 values on the y-axis displayed
  scale_colour_manual(name = "Experimental\nCondition", #need to use scale_colour_manual for line graphs
                      values = c("Blue",
                                 "Purple",
                                 "Gold",
                                 "Red"),
                      labels = c(expression(AM[AM]),
                                 expression(AM[GH]),
                                 expression(GH[AM]),
                                 expression(GH[GH])))+
  theme_light()+
  xlab(xlab)+
  ylab("Survivorship")+
  theme(axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 40, colour = "black"),
        axis.text.y = element_text(size = 40, colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 40),
        legend.title = element_blank(),
        legend.title.align = 0.5)
TPC_plot

## try with different linetypes
TPC_plot2 <- ggplot(data = TPC_data_2, aes(temp, #the x-axis value
                                          survivorship, #the y-axis value
                                          color = treatment,
                                          linetype = dev.treatment))+ # change the lines to be different based on sex
  geom_point(alpha = 0.5,
             size = 4,
             position = position_jitter(height = 0.025, width = 0.5))+ # make a scatter plot with the points slightly offset from each other
  
  geom_smooth(method = "glm", # make the plotted model match the model we use to test
              method.args = list(family = "binomial"),# make sure we know that there is only 2 possible data values
              se = FALSE,
              size = 2)+ # do not plot standard error
  
  scale_y_continuous(breaks = c(0,0.5,1))+ # 3 values on the y-axis displayed
  
  scale_colour_manual(name = "Lineage", #need to use scale_colour_manual for line graphs
                      values = c("Blue",
                                 "Red"),
                      labels = c("AM", "GH"))+
  
  scale_linetype_discrete(name = "Dev. Env.",
                          labels = c("AM", "GH"))+
  
  theme_light()+
  xlab(xlab)+
  ylab("Survivorship")+
  theme(axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 40, colour = "black"),
        axis.text.y = element_text(size = 40, colour = "black"),
        legend.position = "right",
        legend.text = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.title.align = 0.5)
TPC_plot2
TPC_plot+facet_wrap(.~environment, nrow = 2) # create two graphs with one on top of the other (i.e. 2 rows instead of 2 columns)


## Compare HH to AA transplants between gens 1 and 9
TPC_data.HHAA <- TPC_data[(TPC_data$treatment == "HH") & (TPC_data$dev.treatment == "AA"),]


TPC_data.HHAA$generation <- as.factor(TPC_data.HHAA$generation)



TPC_plot_HHAA <- ggplot(data = TPC_data.HHAA, aes(temp, #the x-axis value
                                                  survivorship, #the y-axis value
                                                  color = generation))+ # change the lines to be different based on sex
  geom_point(alpha = 0.5, 
             size = 4,
             position = position_jitter(height = 0.025, width = 0.5))+ # make a scatter plot with the points slightly offset from each other
  
  geom_smooth(method = "glm", # make the plotted model match the model we use to test
              method.args = list(family = "binomial"),# make sure we know that there is only 2 possible data values
              se = FALSE,
              size = 2)+ # do not plot standard error
  
  #scale_y_continuous(breaks = c(0,0.5,1))+ # 3 values on the y-axis displayed
  scale_colour_manual(name = "Lineage", #need to use scale_colour_manual for line graphs
                      values = c(#"Blue",
                        #"Purple",
                        "Gold",
                        "Gray"),
                      labels = c(#"Ambient",
                        # "Ambient in GW",
                        "F1",
                        "F9"))+
  scale_y_continuous(breaks = c(0,0.5,1))+
  theme_light()+
  xlab(xlab)+
  ylab("Survivorship")+
  theme(axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 40, colour = "black"),
        axis.text.y = element_text(size = 40, colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 40),
        legend.title = element_text(size= 40),
        legend.title.align = 0.5)
TPC_plot_HHAA

##### AM AMGH and GH TPC plot #####

TPC_data.AAHH.AA.HH <- rbind(TPC_data.AAdev, TPC_data.AAHH)



TPC_plot_AAHH.AA.HH <- ggplot(data = TPC_data.AAHH.AA.HH, aes(temp, #the x-axis value
                                                              survivorship, #the y-axis value
                                                              color = unique))+ # change the lines to be different based on sex
  geom_point(alpha = 0.5,
             size = 4,
             position = position_jitter(height = 0.025, width = 0.5))+ # make a scatter plot with the points slightly offset from each other
  
  geom_smooth(method = "glm", # make the plotted model match the model we use to test
              method.args = list(family = "binomial"),# make sure we know that there is only 2 possible data values
              se = FALSE,
              size =4)+ # do not plot standard error
  
  #scale_y_continuous(breaks = c(0,0.5,1))+ # 3 values on the y-axis displayed
  scale_colour_manual(name = "Generation", #need to use scale_colour_manual for line graphs
                      values = c("Blue",
                                 "Pink",
                                 "Red"),
                      #"Gold",
                      #"Red"),
                      labels = c(expression(AM[AM]),
                                 expression(AM[GH]),
                                 expression(GH[GH])))+
  #"F1",
  #"F9"))+
  theme_light()+
  xlab(xlab)+
  ylab("Survivorship")+
  theme(axis.title.x = element_text(size = 60),
        axis.title.y = element_text(size = 60),
        axis.text.x = element_text(size = 60, colour = "black"),
        axis.text.y = element_text(size = 60, colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 60),
        legend.title = element_blank(),
        legend.title.align = 0.5)
TPC_plot_AAHH.AA.HH


### Blue and red lines only

TPC_data_blue_red <- filter(TPC_data, unique == c("AA_AA", "HH_HH"))

TPC_plot_blue_red <- ggplot(data = TPC_data_blue_red, aes(temp, #the x-axis value
                                                          survivorship, #the y-axis value
                                                          color = unique))+ # change the lines to be different based on sex
  geom_point(position = position_jitter(height = 0.025, width = 0.5))+ # make a scatter plot with the points slightly offset from each other
  
  geom_smooth(method = "glm", # make the plotted model match the model we use to test
              method.args = list(family = "binomial"),# make sure we know that there is only 2 possible data values
              se = FALSE)+ # do not plot standard error
  
  scale_y_continuous(breaks = c(0,0.5,1))+ # 3 values on the y-axis displayed
  scale_colour_manual(name = "Experimental\nCondition", #need to use scale_colour_manual for line graphs
                      values = c("Blue",
                                 #"Purple",
                                 #"Gold",
                                 "Red"),
                      labels = c("Ambient",
                                 # "Ambient in GW",
                                 # "GW in Ambient",
                                 "GW"))+
  theme_light()+
  xlab(xlab)+
  ylab("Survivorship")+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.title.align = 0.5)
TPC_plot_blue_red

### Red and gold lines

TPC_data_red_yellow <- filter(TPC_data, unique == c("HH_AA", "HH_HH"))

TPC_plot_red_yellow <- ggplot(data = TPC_data_red_yellow, aes(temp, #the x-axis value
                                                              survivorship, #the y-axis value
                                                              color = unique))+ # change the lines to be different based on sex
  geom_point(position = position_jitter(height = 0.025, width = 0.5))+ # make a scatter plot with the points slightly offset from each other
  
  geom_smooth(method = "glm", # make the plotted model match the model we use to test
              method.args = list(family = "binomial"),# make sure we know that there is only 2 possible data values
              se = FALSE)+ # do not plot standard error
  
  scale_y_continuous(breaks = c(0,0.5,1))+ # 3 values on the y-axis displayed
  scale_colour_manual(name = "Experimental\nCondition", #need to use scale_colour_manual for line graphs
                      values = c(#"Blue",
                        #"Purple",
                        "Gold",
                        "Red"),
                      labels = c(#"Ambient",
                        # "Ambient in GW",
                        "GW in Ambient",
                        "GW"))+
  theme_light()+
  xlab(xlab)+
  ylab("Survivorship")+
  theme(axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20, colour = "black"),
        axis.text.y = element_text(size = 20, colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.title.align = 0.5)
TPC_plot_red_yellow

TPC_glm <- glm(survivorship~treatment*dev.treatment*sex, data = TPC_data)
summary(TPC_glm)















############ Statistics ##############

setwd("C:/users/james/Documents/Grad_school/OA_Project/Cost_experiments/TPC/Statistics")
#catherine notes

#str(TPC_data) # structure of the data frame

## make all the variables as "factors"... This will change the results of the ANOVA
TPC_data$replicate<-as.factor(TPC_data$replicate)
#TPC_data$temp <- as.factor(TPC_data$temp)
TPC_data$treatment <- as.factor(TPC_data$treatment)
TPC_data$dev.treatment <- as.factor(TPC_data$dev.treatment)


TPC_data$sex <- as.factor(TPC_data$sex)

# construct the model 
# linear mixed model with replicates designated as random effects
m1<-glmer(survivorship~treatment*dev.treatment*temp*generation*sex + # do not include sex since it has no effect in a sex*temp anova
            (1|replicate), #include the replicate as a random effect
          family=binomial, #indicate binomial data
          data = TPC_data)


plot_model(m1)
tab_model(m1)


TPC.Anova <- Anova(m1, type = 2) # use type 3 because interactions ARE significant and the data is unbalanced. Type 3 therefore gives a more robust test.
TPC.Anova
summary(TPC.Anova)
TPC.Anova$factors <- rownames(TPC.Anova)
fwrite(TPC.Anova, file = "TPC_Anova_total_type2.txt", sep = "\t")


m1a <- glmer(survivorship~treatment*dev.treatment*temp + # do not include sex since it has no effect in a sex*temp anova
               (1|replicate), #include the replicate as a random effect
             family=binomial, #indicate binomial data
             data = TPC_data)

TPC.Anova_no_gen_no_sex <- Anova(m1a, type = 2) # use type 3 because interactions ARE significant and the data is unbalanced. Type 3 therefore gives a more robust test.
TPC.Anova_no_gen_no_sex
summary(TPC.Anova_no_gen_no_sex)
TPC.Anova_no_gen_no_sex$factors <- rownames(TPC.Anova_no_gen_no_sex)
fwrite(TPC.Anova_no_gen_no_sex, file = "TPC_Anova_total_no_gen_no_sex_type2.txt", sep = "\t")

#TPC.emm <- emmeans(m1, pairwise ~ treatment | dev.treatment)


#pairs(TPC.emm)


#library(multcomp)
#m1.posthoc <- glht(m1)
#summary(m1.posthoc)

#TukeyHSD(m1)



############################################################################################################################################################################################
###### post-hoc models for each of the orignial treatments with Type II ANOVAs


TPC_data.HH <- filter(TPC_data, TPC_data$treatment == "HH")


m2<-glmer(survivorship~generation*temp + 
            (1|replicate), #include the replicate as a random effect
          family=binomial, #indicate binomial data
          data = TPC_data.HH)

TPC.Anova2 <- Anova(m2, type = 2) # use type 2 because interactions are NOT significant. Type 2 therefore gives a more robust test.
TPC.Anova2$factors <- rownames(TPC.Anova2)

fwrite(TPC.Anova2, file = "TPC_gen_anova_type2.txt", sep = "\t")



TPC_data.AA <- filter(TPC_data, TPC_data$treatment == "AA")
m3<-glmer(survivorship~dev.treatment*temp + 
            (1|replicate), #include the replicate as a random effect
          family=binomial, #indicate binomial data
          data = TPC_data.AA)
TPC.Anova3 <- Anova(m3, type = 2) # use type 2 because interactions are NOT significant. Type 2 therefore gives a more robust test.
TPC.Anova3$factors <- rownames(TPC.Anova3)

fwrite(TPC.Anova3, file = "TPC_AA_anova_type2.txt", sep = "\t")

TPC_data.orig <- filter(TPC_data, TPC_data$treatment == TPC_data$dev.treatment)
m4<-glmer(survivorship~treatment*temp + 
            (1|replicate), #include the replicate as a random effect
          family=binomial, #indicate binomial data
          data = TPC_data.orig)

TPC.Anova4 <- Anova(m4, type = 2) # use type 2 because interactions are NOT significant. Type 2 therefore gives a more robust test.
TPC.Anova4$factors <- rownames(TPC.Anova4)

fwrite(TPC.Anova4, file = "TPC_orig_anova_type2.txt", sep = "\t")

TPC_data.HHF1 <- filter(TPC_data, TPC_data$treatment == "HH" & TPC_data$generation == 1)
m5<-glmer(survivorship~dev.treatment*temp + 
            (1|replicate), #include the replicate as a random effect
          family=binomial, #indicate binomial data
          data = TPC_data.HHF1)

TPC.Anova5 <- Anova(m5, type = 2) # use type 2 because interactions are NOT significant. Type 2 therefore gives a more robust test.
TPC.Anova5$factors <- rownames(TPC.Anova5)

fwrite(TPC.Anova5, file = "TPC_HHF1_anova_type2.txt", sep = "\t")

TPC_data.rescue <- filter(TPC_data, TPC_data$treatment == "AA" & TPC_data$dev.treatment == "AA")
TPC_data.rescue2 <- filter(TPC_data, TPC_data$generation == 9)
TPC_data.rescue3 <- rbind(TPC_data.rescue, TPC_data.rescue2)


m6<-glmer(survivorship~treatment*temp + 
            (1|replicate), #include the replicate as a random effect
          family=binomial, #indicate binomial data
          data = TPC_data.rescue3)
TPC.Anova6 <- Anova(m6, type = 2) # use type 2 because interactions are NOT significant. Type 2 therefore gives a more robust test.
TPC.Anova6$factors <- rownames(TPC.Anova6)

fwrite(TPC.Anova6, file = "TPC_rescue_anova_type2.txt", sep = "\t")



TPC_data.AAHH <- filter(TPC_data, TPC_data$dev.treatment == "HH")
m7<-glmer(survivorship~treatment*temp + 
            (1|replicate), #include the replicate as a random effect
          family=binomial, #indicate binomial data
          data = TPC_data.AAHH)
TPC.Anova7 <- Anova(m7, type = 2) # use type 2 because interactions are NOT significant. Type 2 therefore gives a more robust test.
TPC.Anova7$factors <- rownames(TPC.Anova7)

fwrite(TPC.Anova7, file = "TPC_HHdev_anova_type2.txt", sep = "\t")

TPC_data$transplant <- if_else(TPC_data$treatment == TPC_data$dev.treatment, 0, 1)

TPC_data.AAdev <- filter(TPC_data, TPC_data$treatment == "AA", TPC_data$dev.treatment == "AA")
TPC_data.AAdev2 <- filter(TPC_data, TPC_data$transplant == 1 & TPC_data$generation == 1)
TPC_data.AAdev3 <- rbind(TPC_data.AAdev, TPC_data.AAdev2)

m8<-glmer(survivorship~transplant*temp + 
            (1|replicate), #include the replicate as a random effect
          family=binomial, #indicate binomial data
          data = TPC_data.AAdev3)
TPC.Anova8 <- Anova(m8, type = 2) # use type 2 because interactions are NOT significant. Type 2 therefore gives a more robust test.
TPC.Anova8$factors <- rownames(TPC.Anova8)

fwrite(TPC.Anova8, file = "TPC_AAtrans_anova_type2.txt", sep = "\t")

m9<-glmer(survivorship~treatment*temp + 
            (1|replicate), #include the replicate as a random effect
          family=binomial, #indicate binomial data
          data = TPC_data.AAdev2)
TPC.Anova9 <- Anova(m9, type = 2) # use type 2 because interactions are NOT significant. Type 2 therefore gives a more robust test.
TPC.Anova9$factors <- rownames(TPC.Anova9)

fwrite(TPC.Anova9, file = "TPC_transplants_anova_type2.txt", sep = "\t")

TPC_data.Gen1 <- filter(TPC_data, TPC_data$generation == 1)

m10<-glmer(survivorship~treatment*dev.treatment*temp*sex + 
             (1|replicate), #include the replicate as a random effect
           family=binomial, #indicate binomial data
           data = TPC_data.Gen1)
TPC.Anova10 <- Anova(m10, type = 2) # use type 2 because interactions are NOT significant. Type 2 therefore gives a more robust test.
TPC.Anova10$factors <- rownames(TPC.Anova10)

fwrite(TPC.Anova10, file = "TPC_no_gen_anova_type2.txt", sep = "\t")


# test to see if transplants at max performance are different from each other - post-hoc pairwise comparison
TPC.Tukey.9 <- TukeyHSD(aov(survivorship~treatment*as.factor(temp), data = TPC_data.AAdev2))
TPC.Tukey.9

TPC_data.AAdev2$treatment <- as.factor(TPC_data.AAdev2$treatment)
TPC_data.AAdev2$temp <- as.factor(TPC_data.AAdev2$temp)
TPC.pairwise.9 <- pairwise.t.test(TPC_data.AAdev2$survivorship, TPC_data.AAdev2$treatment:TPC_data.AAdev2$temp, p.adj = "none")
TPC.pairwise.9 <- tidy(TPC.pairwise.9)

TPC.Anova9$factors <- rownames(TPC.Anova9)

fwrite(TPC.Anova9, file = "TPC_trans_anova.txt", sep = "\t")
fwrite(TPC.pairwise.9, file = "TPC_trans_pairwise.txt", sep = "\t")


##### Try finding the LD50 values ######


library(MASS)


LD50 <- dose.p(m1, p = 0.5)



## Try using a list of models

TPC.list <- split(TPC_data, f=TPC_data$unique)

# apply the "glm" function (create models) to data in TPC.list
TPC.model.list <- lapply(TPC.list, function (x) glm(survivorship~temp, data = x, family = binomial)) 


# find the 50 percent probabilities for each member of the TPC.model.list
TPC.LD50 <- lapply(TPC.model.list, function (x) dose.p(x, p=0.5))




TPC.LD50$AA_AA_1
TPC.LD50$AA_HH_1
TPC.LD50$HH_AA_1
TPC.LD50$HH_AA_9
TPC.LD50$HH_HH_1

# create vectors for the data frame
treatment <- c("AM", "AM", "GH", "GH", "GH")
dev.treatment <- c( "AM", "GH", "AM", "AM", "GH")
generation <- c("1", "1", "1", "9", "1")
LD50 <- c(29.33632, 28.12996, 28.12321, 29.84187, 27.31829)
stand.err <- c(0.3496183, 0.5283111, 0.3112534, 0.6210961, 0.4229522)

# put the vectors in a data frame
LD50.df <- data.frame(treatment, dev.treatment, generation, LD50, stand.err)


Rxn.norm.df <- LD50.df[-4,]

library(ggpubr)

LD50.reaction.norm <- ggplot(Rxn.norm.df, aes(dev.treatment, LD50, color = factor(treatment), group = treatment))+
  geom_point(size = 5)+
  geom_line(size = 2)+
  geom_errorbar(aes(ymin = LD50-stand.err, ymax = LD50+stand.err, 
                   
                    #size = 10,
                    
                    width = 0.1), 
                #position = position_dodge(width = 0.5)
  )+
  theme_classic()+
  
  labs(y=expression(atop(LD[50]," stress"*" temperature"*" (°C)")),
       x="Developmental Environment")+
  scale_color_manual(name="Lineage",
                     labels=c("AM",
                              "GH"),
                     values = c("blue",
                                "red"))+
  theme(legend.title = element_text(colour = "black", size=40),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 40),
        axis.text = element_text(size = 40),
        legend.text = element_text(size = 40))


LD50.reaction.norm




xlab

# find the max survivorship for each model
TPC.18 <- filter(TPC_data, TPC_data$temp == 18)

TPC.max <- TPC.18 %>%
  group_by(unique) %>%
  summarise(mean = mean(survivorship, na.rm = TRUE),
            sd = sd(survivorship, na.rm = TRUE),
            n = n()) %>%
  mutate(se = sd/sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n-1)*se,
         upper.ci = mean + qt(1 - (0.05 / 2), n-1)*se)

TPC.max <- separate(TPC.max, "unique", into = c("treatment", "dev.treatment", "generation"))







TPC.list.18 <- split(TPC.18, f=TPC.18$unique)

# apply the "glm" function (create models) to data in TPC.list
TPC.model.list.18 <- lapply(TPC.list.18, function (x) glm(survivorship~temp, data = x, family = binomial)) 



TPC.max <- lapply(TPC.model.list.18, function (x) dose.p(x, p=max(x$fitted.values)))

TPC.max$AA_AA_1
TPC.max$AA_HH_1
TPC.max$HH_AA_1
TPC.max$HH_AA_9
TPC.max$HH_HH_1

## changed the numbers to reflect percentages (i.e. multiplied by 100)
maximum <- c(098.62849, 089.19029, 098.75173, 090.53305, 092.73588)
max.stand.err <- c(002.502318, 005.692939, 003.472660, 008.333333, 004.718416)


# put the vectors in a data frame
Max.df <- data.frame(treatment, dev.treatment, generation, maximum, max.stand.err)

Max.df.plot <- Max.df[-4,]

Max.reaction.norm <- ggplot(Max.df.plot, aes(dev.treatment, maximum, color = factor(treatment), group = treatment))+
  geom_point(size = 5)+
  geom_line(size = 2)+
  geom_errorbar(aes(ymin = maximum-max.stand.err, ymax = maximum+max.stand.err, width = 0.1), 
                #position = position_dodge(width = 0.5)
  )+
  theme_classic()+
  
  labs(y=expression(atop("Maximum", "Survivorship (%)")),
       x="Developmental Treatment")+
  scale_color_manual(name="Lineage",
                     labels=c("AM",
                              "GH"),
                     values = c("blue",
                                "red"))+
  theme(legend.title = element_text(colour = "black", size=28),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.title = element_text(size = 40),
        axis.text = element_text(size = 40),
        legend.text = element_text(size = 40))

Max.reaction.norm


fwrite(Max.df, file = "Max_surv_temp.txt", sep = "\t")
