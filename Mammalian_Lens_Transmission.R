### Veilleux Lens Transmission & Spatial Resolution Phylogenetic Analyses Code ###

#packages
library(ape)
library(dplyr)
library(geiger)
library(phytools)
library(phylolm)
library(car)
library(ggplot2)
library(nationalparkcolors) ## for colors in these figures, any palette can be substituted

##### 1. Calculate %UVA and 50T nm ####
## %UVA transmitted and 50%T (nm) calculated from McDowell et al. 2024 lens spectra,
## computed using Douglas & Jeffery 2014 definitions:
##   %UVA = area under the transmission curve 315-400 nm, as a % of the
##          85 nm x 100% rectangle
##   50%T = wavelength where transmission first reaches 50%

dat <- read.csv("McDowell2024_UVA_data_subset.csv",stringsAsFactors=T,header=T)

uva <- dat$wl >= 315 & dat$wl <= 400. #set UVA
pct_uva <- colSums(dat[uva, -1]) / 85 #calculate percent UVA

t50 <- sapply(dat[, -1], function(T) dat$wl[which(T >= 50)[1]]) #identify 50.T.nm
t50[t50 == min(dat$wl)] <- NA
t50

res <- data.frame(species = names(pct_uva), pct_uva = round(pct_uva, 1), t50 = t50)

write.csv(res, "McDowell_UVA_T50.csv", row.names = FALSE)

## Calculate for hamster digitized data ##
## hamster data from WebPlotDigitizer
hamster<-read.csv("Hamster_digitized_data.csv",header=F,col.names=c("wl","T"))

hamster<-hamster[order(hamster$wl),]
wl<-315:700
Transmission<-approx(hamster$wl,hamster$T,xout=wl)$y

#PropUV
sum(Transmission[wl >=315 & wl <=400])/85

#50Tnm
wl[which(Transmission>=50)[1]]

##### 2. Compare with D&J's values for shared species #####
# Exclude column: New = new taxa not in D&J 2014; Y = exclude from correlation
comps<-read.csv("McDowell_UVA_T50_analysis.csv",header=T, stringsAsFactors = T)

D.J<-subset(comps,Exclude!="New")
summary(D.J)

no.bat.rat<-subset(D.J,Exclude!="Y")
cor.test(no.bat.rat$pct_uva,no.bat.rat$D.J_pctUVA)
cor.test(no.bat.rat$t50,no.bat.rat$D.J_t50)


### 3. Prepare Dataset and Trees ####
mam.dat<-read.csv("Mammal_UVA_Dataset.csv",header=T,stringsAsFactors = T)

#order activity pattern to compare with different reference categories
mam.dat$AP.2[mam.dat$AP=="N"]<-"1.Nocturnal"
mam.dat$AP.2[mam.dat$AP=="C"]<-"2.Cathemeral"
mam.dat$AP.2[mam.dat$AP=="D"]<-"3.Diurnal"

mam.dat$AP.3[mam.dat$AP=="N"]<-"3.Nocturnal"
mam.dat$AP.3[mam.dat$AP=="C"]<-"2.Cathemeral"
mam.dat$AP.3[mam.dat$AP=="D"]<-"1.Diurnal"

#convert to data frame
data.frame(mam.dat[2:14]) -> df.mam
rownames(df.mam) <- mam.dat[,1]

#Transform variables: log transform visual acuity (VA) and logit-transform proportion lens transmission
df.mam <- df.mam %>%
  mutate(log_VA  = log(VA),
         logit_Lens = logit(PropUVA))

#set levels of Acuity Type
df.mam$Acuity.Type<-factor(df.mam$Acuity.Type,levels = c("A", "B"))

summary(df.mam)

df.mam$AP.4[df.mam$AP=="N"]<-"Non-diurnal"
df.mam$AP.4[df.mam$AP=="C"]<-"Non-diurnal"
df.mam$AP.4[df.mam$AP=="D"]<-"Diurnal"
df.mam$AP.4<-as.factor(df.mam$AP.4)
df.mam$AP.3<-as.factor(df.mam$AP.3)
df.mam$AP.2<-as.factor(df.mam$AP.2)

###### Define Data Subsets #####
nohap <- subset(df.mam, Haplorhine != "Y")   # non-haplorhines
dark  <- subset(df.mam, AP != "D")           # non-diurnal
day   <- subset(df.mam, AP == "D")           # diurnal
dayNH <- subset(day, Haplorhine != "Y")      # diurnal non-haplorhines (post hoc)


###### Load Phylogenetic Trees###### 
trees_1000_trimmed<-read.nexus("mammal_trees_1000_Nov2026.nex")

#check overlap that it really trim to the 47 species on a random tree
tree489<-trees_1000_trimmed[[489]]
overlap <- name.check(tree489, df.mam)
overlap #should say OK


### 4. Perform Phylogenetic Regressions ####

###### Function for running phylogenetic analysis across trees ######
run_phylogenetic_analysis <- function(
    trees, formula, data, model_type,
    lower_bound = 1e-15, upper_bound = 100
) {
  results <- lapply(seq_along(trees), function(i) {
    tr <- trees[[i]]
    
    # Diagnostic wrapper for each tree
    fit <- try(phylolm(formula, data = data, phy = tr, model = model_type,
                       lower.bound = lower_bound, upper.bound = upper_bound),
               silent = TRUE)
    
    if (inherits(fit, "try-error")) {
      cat("Tree", i, "failed:", conditionMessage(attr(fit, "condition")), "\n")
      return(NULL)
    }
    
    if (is.null(fit)) return(NULL)
    
    # Extract coefficients and stats
    summ <- summary(fit)
    ci_lower <- summ$coefficients[,1] - 1.96 * summ$coefficients[,2]
    ci_upper <- summ$coefficients[,1] + 1.96 * summ$coefficients[,2]
    
    # Extract AIC values if they exist
    model.aic <- AIC(fit)
    
    # Extract R2 values if they exist
    r2 <- ifelse(!is.null(summ$r.squared), summ$r.squared, NA)
    r2_adj <- ifelse(!is.null(summ$adj.r.squared), summ$adj.r.squared, NA)
    
    # Extract model-specific parameter
    if (model_type == "BM") {
      model_param <- fit$sigma2
      param_name <- "sigma2"
    } else if (model_type == "lambda") {
      model_param <- fit$optpar
      param_name <- "lambda"
    } else if (model_type == "OUrandomRoot") {
      model_param <- fit$optpar
      param_name <- "alpha"
    } else {
      model_param <- fit$optpar
      param_name <- "optpar"
    }
    
    data.frame(
      tree = i,
      term = rownames(summ$coefficients),
      estimate = summ$coefficients[,1],
      se = summ$coefficients[,2],
      t = summ$coefficients[,3],
      p = summ$coefficients[,4],
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      model_param = model_param,
      param_name = param_name,
      R2 = r2,
      R2_adj = r2_adj,
      model.aic = model.aic,
      stringsAsFactors = FALSE
    )
  })
  
  results_clean <- results[!sapply(results, is.null)]
  total_trees <- length(trees)
  successful_trees <- length(results_clean)
  cat("Summary:", successful_trees, "of", total_trees, "trees fit successfully.\n")
  
  if (successful_trees == 0) return(NULL)
  do.call(rbind, results_clean)
}

###### Function for summarizing phylogenetic results ######
summarize_phylogenetic_results <- function(results_df) {
  if (is.null(results_df)) return(NULL)
  
  coef_summary <- results_df %>%
    group_by(term) %>%
    summarize(
      mean_estimate = mean(estimate, na.rm = TRUE),
      median_estimate = median(estimate, na.rm = TRUE),
      sd_estimate = sd(estimate, na.rm = TRUE),
      mean_p = mean(p, na.rm = TRUE),
      median_p = median(p, na.rm = TRUE),
      prop_sig = mean(p < 0.05, na.rm = TRUE),
      mean_ci_lower = mean(ci_lower, na.rm = TRUE),
      mean_ci_upper = mean(ci_upper, na.rm = TRUE),
      median_ci_lower = median(ci_lower, na.rm = TRUE),
      median_ci_upper = median(ci_upper, na.rm = TRUE),
      n_trees = n()
    )
  
  param_summary <- results_df %>%
    group_by(tree) %>%
    slice(1) %>%
    ungroup() %>%
    summarize(
      mean_param = mean(model_param, na.rm = TRUE),
      median_param = median(model_param, na.rm = TRUE),
      sd_param = sd(model_param, na.rm = TRUE),
      param_name = first(param_name)
    )
  
  AIC_summary <- results_df %>%
    group_by(tree) %>%
    slice(1) %>%
    ungroup() %>%
    summarize(
      mean_AIC = mean(model.aic, na.rm = TRUE),
      median_AIC = median(model.aic, na.rm = TRUE),
      sd_AID = sd(model.aic, na.rm = TRUE)
    )
  
  r2_summary <- results_df %>%
    group_by(tree) %>%
    slice(1) %>%  # one row per tree
    ungroup() %>%
    summarize(
      mean_R2 = mean(R2, na.rm = TRUE),
      sd_R2 = sd(R2, na.rm = TRUE),
      mean_R2_adj = mean(R2_adj, na.rm = TRUE),
      median_R2_adj = median(R2_adj, na.rm = TRUE),
      sd_R2_adj = sd(R2_adj, na.rm = TRUE)
    )
  
  list(
    coefficients = coef_summary,
    model_params = param_summary,
    AIC_summary = AIC_summary,
    r2_summary = r2_summary,
    successful_trees = length(unique(results_df$tree))
  )
}

#combine the two functions for one
run_model <- function(formula, data) {
  vars <- all.vars(formula)
  d <- data[complete.cases(data[, vars]), ]
  res <- run_phylogenetic_analysis(trees_1000_trimmed, formula, d,
                                   model_type = "lambda",
                                   lower_bound = 0, upper_bound = 1)
  summarize_phylogenetic_results(res)
}


###### Analysis 1: Lens transmission metrics and activity pattern #######
#Table S4 - Activity pattern
tS4_prop <- run_model(logit_Lens ~ AP.3, df.mam)
tS4_propB <- run_model(logit_Lens ~ AP.2, df.mam) # change AP reference to Nocturnal

tS4_t50  <- run_model(X50.T.nm ~ AP.3, df.mam)
tS4_t50B  <- run_model(X50.T.nm ~ AP.2, df.mam) # change AP reference to Nocturnal

#Table S5 - Activity pattern + path length
tS5_prop <- run_model(logit_Lens ~ AP.3+Path.Length, df.mam) #note, warnings of dropping 4 taxa from tree due to missing data
tS5_propB <- run_model(logit_Lens ~ AP.2+Path.Length, df.mam) # change AP reference to Nocturnal

tS5_t50 <- run_model(X50.T.nm ~ AP.3+Path.Length, df.mam) #note, warnings of dropping 4 taxa from tree due to missing data
tS5_t50B <- run_model(X50.T.nm ~ AP.2+Path.Length, df.mam)  #change AP reference to Nocturnal


###### Analysis 2: Lens transmission metrics and visual acuity: Table 1 ###### 
# note, there will be warnings to let you know that the model is dropping taxa with missing data

t1_all_prop <- run_model(logit_Lens ~ log_VA + AP.4 + Acuity.Type, df.mam)
t1_nohap_prop <- run_model(logit_Lens ~ log_VA + AP.4 + Acuity.Type, nohap)
t1_dark_prop  <- run_model(logit_Lens ~ log_VA + Acuity.Type, dark)
t1_day_prop  <- run_model(logit_Lens ~ log_VA + Acuity.Type, day)

t1_all_t50 <- run_model(X50.T.nm ~ log_VA + AP.4 + Acuity.Type, df.mam)
t1_nohap_t50 <- run_model(X50.T.nm ~ log_VA + AP.4 + Acuity.Type, nohap)
t1_dark_t50  <- run_model(X50.T.nm ~ log_VA + Acuity.Type, dark)
t1_day_t50  <- run_model(X50.T.nm ~ log_VA + Acuity.Type, day)


###### Analysis 3: Relative cornea size and acuity: Table S6 ###### 

tS6 <- run_model(log_VA ~ Rel.Cornea + Acuity.Type, df.mam)


###### Analysis 4: Lens transmission metrics and relative cornea size: Table 2 ###### 

t2_all_prop <- run_model(logit_Lens ~ Rel.Cornea + AP.4, df.mam)
t2_nohap_prop <- run_model(logit_Lens ~ Rel.Cornea + AP.4, nohap)
t2_dark_prop  <- run_model(logit_Lens ~ Rel.Cornea, dark)
t2_day_prop  <- run_model(logit_Lens ~ Rel.Cornea, day)

t2_all_t50 <- run_model(X50.T.nm ~ Rel.Cornea + AP.4, df.mam)
t2_nohap_t50 <- run_model(X50.T.nm ~ Rel.Cornea + AP.4, nohap)
t2_dark_t50  <- run_model(X50.T.nm ~ Rel.Cornea, dark)
t2_day_t50  <- run_model(X50.T.nm ~ Rel.Cornea, day)


#### 5. Model Diagnostic Checks #####

###### Check for multicollinearity: run a regular LM to check VIF ###### 
# VIF > 5-10 indicates problematic multicollinearity
vif(lm(X50.T.nm ~ log_VA + Acuity.Type + AP.4, data = df.mam))  # Table 1 model
vif(lm(logit_Lens ~ AP.3 + Path.Length, data = df.mam))         # Table S5 model
vif(lm(log_VA ~ Rel.Cornea + Acuity.Type, data = df.mam))       # Table S6 model
vif(lm(log_VA ~ Rel.Cornea + AP.4, data = df.mam))


###### Check residuals for normality and homoscedasticity ###### 
#use a sample tree

check_resid <- function(formula, data, tree = trees_1000_trimmed[[420]]) {
  d   <- data[complete.cases(data[, all.vars(formula)]), ]
  fit <- phylolm(formula, data = d, phy = tree, model = "lambda")
  tr  <- transf.branch.lengths(tree, model = "lambda",
                               parameters = list(lambda = fit$optpar))$tree
  r   <- residuals(fit)
  V   <- ape::vcv(tr)[names(r), names(r)]
  r_c <- solve(t(chol(V)), r)
  r_c <- r_c / sd(r_c)
  par(mfrow = c(1, 2))
  plot(fit$fitted.values, r_c, xlab = "Fitted", ylab = "Phylo-corrected residuals",
       main = deparse(formula))
  qqnorm(r_c); qqline(r_c)
  c(max = names(r_c)[which.max(r_c)], min = names(r_c)[which.min(r_c)])
}

check_resid(logit_Lens ~ AP.3, df.mam)                          # Table S4
check_resid(X50.T.nm ~ AP.3, df.mam)                          # Table S4

check_resid(logit_Lens ~ AP.3 + Path.Length, df.mam)            # Table S5
check_resid(X50.T.nm ~ AP.3 + Path.Length, df.mam)            # Table S5

check_resid(logit_Lens ~ log_VA + Acuity.Type + AP.4, df.mam)   # Table 1, all taxa
check_resid(logit_Lens ~ log_VA + Acuity.Type + AP.4, nohap)   # Table 1, non-haplorhines
check_resid(logit_Lens ~ log_VA + Acuity.Type, dark)            # Table 1, non-diurnal
check_resid(logit_Lens ~ log_VA + Acuity.Type, day)            # Table 1, diurnal only

check_resid(X50.T.nm ~ log_VA + Acuity.Type + AP.4, df.mam)   # Table 1, all taxa
check_resid(X50.T.nm ~ log_VA + Acuity.Type + AP.4, nohap)   # Table 1, non-haplorhines
check_resid(X50.T.nm ~ log_VA + Acuity.Type, dark)            # Table 1, non-diurnal
check_resid(X50.T.nm ~ log_VA + Acuity.Type, day)            # Table 1, diurnal only

check_resid(logit_Lens ~ Rel.Cornea + AP.4, df.mam)   # Table 2, all taxa
check_resid(logit_Lens ~ Rel.Cornea + AP.4, nohap)   # Table 2, non-haplorhines
check_resid(logit_Lens ~ Rel.Cornea, dark)            # Table 2, non-diurnal
check_resid(logit_Lens ~ Rel.Cornea, day)            # Table 2, diurnal only

check_resid(X50.T.nm ~ Rel.Cornea + AP.4, df.mam)   # Table 2, all taxa
check_resid(X50.T.nm ~ Rel.Cornea + AP.4, nohap)   # Table 2, non-haplorhines
check_resid(X50.T.nm ~ Rel.Cornea, dark)            # Table 2, non-diurnal
check_resid(X50.T.nm ~ Rel.Cornea, day)            # Table 2, diurnal only

###### Testing Effects of Specific Taxa ###### 
#Does Bubalus outlier in Activity Pattern analyses influence results?
df.mam.noBd<-df.mam[rownames(df.mam) != "Bubalus_depressicornis",]
alt_t4_prop <- run_model(logit_Lens ~ AP.3, df.mam.noBd)
alt_t4_T50 <- run_model(X50.T.nm ~ AP.3, df.mam.noBd)

alt_t5_prop <- run_model(logit_Lens ~ AP.3+Path.Length, df.mam.noBd)
alt_t5_T50 <- run_model(X50.T.nm ~ AP.3+Path.Length, df.mam.noBd)

#Rhabodmys acuity measured using looming test
#does excluding it affect direction/significance of results?
df.mam.noRh<-df.mam[rownames(df.mam) != "Rhabdomys_pumilio",]
nohap.noRh<-nohap[rownames(nohap) != "Rhabdomys_pumilio",]
day.noRh<-day[rownames(day) != "Rhabdomys_pumilio",]

alt_t1_all_prop <- run_model(logit_Lens ~ log_VA + AP.4 + Acuity.Type, df.mam.noRh)
alt_t1_nohap_prop <- run_model(logit_Lens ~ log_VA + AP.4 + Acuity.Type, nohap.noRh)
alt_t1_day_prop  <- run_model(logit_Lens ~ log_VA + Acuity.Type, day.noRh)

alt_t1_all_t50 <- run_model(X50.T.nm ~ log_VA + AP.4 + Acuity.Type, df.mam.noRh)
alt_t1_nohap_t50 <- run_model(X50.T.nm ~ log_VA + AP.4 + Acuity.Type, nohap.noRh)
alt_t1_day_t50  <- run_model(X50.T.nm ~ log_VA + Acuity.Type, day.noRh)


#### 6. Making Figures #### 

###### Fig 1. Lens transmission, activity pattern, and spatial resolution ###### 
#Fig 1a: Logit transformed UV transmission by activity pattern
vars <- all.vars(logit_Lens~AP.2)
data.df <- df.mam[complete.cases(df.mam[, vars]),] 
p1<-ggplot(data=df.mam, aes(x=AP.2, y=logit_Lens,fill=AP.2)) + 
  geom_boxplot()+theme_classic(base_size = 14)
p1  +ylab("Logit-Transformed PropUVA")+xlab("Activity Pattern")+scale_fill_manual(values=park_palette("Everglades"))+
  scale_x_discrete(labels=c("Nocturnal","Cathemeral","Diurnal"))+xlab("")+theme(legend.position = "none")

#Fig 1d: λT0.5 and activity pattern
p1<-ggplot(data=df.mam, aes(x=AP.2, y=X50.T.nm,fill=AP.2)) + 
  geom_boxplot()+theme_classic(base_size = 12)
p1  +ylab("λT0.5")+xlab("Activity Pattern")+scale_fill_manual(values=park_palette("Everglades"))+
  scale_x_discrete(labels=c("Nocturnal","Cathemeral","Diurnal"))+xlab("")+theme(legend.position = "none")

#Fig. 1b: Log visual acuity by proportion UV lens transmission, with activity pattern and haplorhine primates coded
vars <- all.vars(logit_Lens~log_VA+AP.2+Acuity.Type)
data.df <- df.mam[complete.cases(df.mam[, vars]),] 
p1<-ggplot(data=data.df, aes(x=(log_VA), y=logit_Lens,col=AP.2,shape=Haplorhine)) + 
  geom_point(size=2.5)+theme_classic(base_size = 12)
p1  +ylab("Logit-Transformed PropUVA")+xlab("Log-Transformed Visual Acuity (cpd)")+
  scale_color_manual(values = park_palette("Everglades"),labels = c("Nocturnal", "Cathemeral","Diurnal"))+ 
  labs(color = "Activity Pattern")

#Fig. 1e: Log visual acuity and 50.T.nm, with activity pattern and haplorhine
vars <- all.vars(X50.T.nm~log_VA+AP.2+Acuity.Type)
data.df <- df.mam[complete.cases(df.mam[, vars]),] 
p1<-ggplot(data=data.df, aes(x=log_VA, y=X50.T.nm,col=AP.2,shape=Haplorhine)) + 
  geom_point(size=2)+theme_classic(base_size = 12)
p1  +ylab("λT0.5")+xlab("Log-Transformed Visual Acuity (cpd)")+ scale_color_manual(values = park_palette("Everglades"),labels = c("Nocturnal", "Cathemeral","Diurnal"))+ labs(color = "Activity Pattern")

#Fig. 1c: relative cornea size X proportion UV lens transmission, with activity pattern and haplorhine coded
vars <- all.vars(logit_Lens~Rel.Cornea+AP.2)
data.df <- df.mam[complete.cases(df.mam[, vars]),] 
p1<-ggplot(data=data.df, aes(x=(Rel.Cornea), y=logit_Lens,col=AP.2,shape=Haplorhine)) + 
  geom_point(size=2)+theme_classic(base_size = 12)
p1  +ylab("Logit-Transformed PropUVA")+xlab("Relative Cornea Size")+ scale_color_manual(values = park_palette("Everglades"),labels = c("Nocturnal", "Cathemeral","Diurnal"))+ labs(color = "Activity Pattern")

#Fig. 1f: Relative cornea size and 50.T.nm, with activity pattern and haplorhine
vars <- all.vars(X50.T.nm~Rel.Cornea+AP.2)
data.df <- df.mam[complete.cases(df.mam[, vars]),] 
p1<-ggplot(data=data.df, aes(x=(Rel.Cornea), y=X50.T.nm,col=AP.2,shape=Haplorhine)) + 
  geom_point(size=2)+theme_classic(base_size = 12)
p1  +ylab("λT0.5")+xlab("Relative Cornea Size")+ scale_color_manual(values = park_palette("Everglades"),labels = c("Nocturnal", "Cathemeral","Diurnal"))+ labs(color = "Activity Pattern")

###### Supplemental Fig S1. Published vs Calculate Transmission ###### 
comps<-read.csv("McDowell_UVA_T50_analysis.csv",header=T, stringsAsFactors = T)
D.J<-subset(comps,Exclude!="New")
p1<-ggplot(data=D.J, aes(x=(pct_uva/100), y=D.J_pctUVA/100)) + 
  geom_point(size=2.5)+theme_classic(base_size = 12)
p1  +ylab("Douglas & Jeffery PropUVA")+xlab("Estimated PropUVA from McDowell data")+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")

p1<-ggplot(data=D.J, aes(x=(t50), y=D.J_t50)) + 
  geom_point(size=2.5)+theme_classic(base_size = 12)
p1  +ylab("Douglas & Jeffery λT0.5")+xlab("Estimated λT0.5 from McDowell data")+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed")

###### Supplemental Fig S2: Lens transmission, activity pattern and path length ###### 
vars <- all.vars(logit_Lens~AP.2+Path.Length)
data.df <- df.mam[complete.cases(df.mam[, vars]),] 
p1<-ggplot(data=data.df, aes(x=(Path.Length), y=logit_Lens,col=AP.2)) + 
  geom_point(size=2)+theme_classic(base_size = 12)
p1  +ylab("Logit-Transformed PropUVA")+xlab("Path Length (mm)")+
  scale_color_manual(values = park_palette("Everglades"),labels = c("Nocturnal", "Cathemeral","Diurnal"))+ 
  labs(color = "Activity Pattern")

vars <- all.vars(X50.T.nm~AP.2+Path.Length)
data.df <- df.mam[complete.cases(df.mam[, vars]),] 
p1<-ggplot(data=data.df, aes(x=(Path.Length), y=X50.T.nm,col=AP.2)) + 
  geom_point(size=2)+theme_classic(base_size = 12)
p1  +ylab("λT0.5")+xlab("Path Length (mm)")+
  scale_color_manual(values = park_palette("Everglades"),labels = c("Nocturnal", "Cathemeral","Diurnal"))+ 
  labs(color = "Activity Pattern")

###### Supplemental Fig S3: Relative cornea size and visual acuity by acuity type ###### 
vars <- all.vars(log_VA~Rel.Cornea+Acuity.Type)
data.df <- df.mam[complete.cases(df.mam[, vars]),] 
p1<-ggplot(data=data.df, aes(x=(Rel.Cornea), y=log_VA,col=AP.2,shape=Acuity.Type)) + 
  geom_point(size=2)+theme_classic(base_size = 12)
pal <- park_palette("Everglades")
p1 + ylab("Log-Transformed Visual Acuity (cpd)") + xlab("Relative Cornea Size") +
  scale_color_manual(values = c("2.Cathemeral" = pal[2], "3.Diurnal" = pal[3]),
                     labels = c(C = "Cathemeral", D = "Diurnal")) +labs(color = "Activity Pattern")
