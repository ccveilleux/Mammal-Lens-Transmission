### Veilleux Lens Transmission - Spatial Resolution Phylogenetic Analyses Code ###

setwd("~/Dropbox/ResearchProjects/LensTransmission/Data")

#packages
library(ape)
library(dplyr)
library(geiger)
library(phytools)
library(phylolm)
library(car)
library(ggplot2)

### 1. Prepare Mammalian Dataset ####
mam.dat<-read.csv("Mammal_UV_Dataset.csv",header=T,stringsAsFactors = T)

#order activity pattern to compare with different reference categories
mam.dat$AP.2[mam.dat$AP=="N"]<-"1.Nocturnal"
mam.dat$AP.2[mam.dat$AP=="C"]<-"2.Cathemeral"
mam.dat$AP.2[mam.dat$AP=="D"]<-"3.Diurnal"

mam.dat$AP.3[mam.dat$AP=="N"]<-"3.Nocturnal"
mam.dat$AP.3[mam.dat$AP=="C"]<-"2.Cathemeral"
mam.dat$AP.3[mam.dat$AP=="D"]<-"1.Diurnal"

#Convert UV lens transmission percentage to a proportion
mam.dat$Prop.Lens.UVA<-mam.dat$Percent.UVA.Transmitted/100

#convert to data frame
data.frame(mam.dat[2:15]) -> df.mam
rownames(df.mam) <- mam.dat[,1]

#Transform variables: log transform visual acuity (VA) and logit-transform proportion lens transmission
df.mam <- df.mam %>%
  mutate(log_VA  = log(VA),
         logit_Lens = logit(Prop.Lens.UVA))

#set levels of Acuity Type
df.mam$Acuity.Type<-factor(df.mam$Acuity.Type,levels = c("A", "B"))

summary(df.mam)

df.mam$AP.4[df.mam$AP=="N"]<-"Non-diurnal"
df.mam$AP.4[df.mam$AP=="C"]<-"Non-diurnal"
df.mam$AP.4[df.mam$AP=="D"]<-"Diurnal"
df.mam$AP.4<-as.factor(df.mam$AP.4)

### 2. Prepare Trees ####
file <- "~/Dropbox/ResearchProjects/LensTransmission/Data/Data_S7_Mammalia_credibleTreeSets_tipDR/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_FBDasZhouEtAl_all10k_v2_nexus.trees"

#total number of trees in file
n_total <- 10000

# Read in nexus file
trees <- read.nexus(file)

# Check length = 10k
length(trees)

set.seed(247)  # reproducible
trees_1000.new <- sample(trees, size = 1000, replace = FALSE)
length(trees_1000.new)  # check it worked: should be 1000
class(trees_1000.new) # check it worked: should be multiPhylo
#fix species binomial headers (some have "X_Anolis", some have family and order)
# Extract binomial name from trees
extract_binomial <- function(tips) {
  sapply(strsplit(tips, "_"), function(x) paste(x[1:2], collapse="_"))
}

# Apply to all trees
trees_1000_all <- lapply(trees_1000.new, function(tr) {
  tr$tip.label <- extract_binomial(tr$tip.label)
  tr
})
#check that it worked
head(trees_1000_all[[1]]$tip.label)

#reset to multiPhylo
class(trees_1000_all) <- "multiPhylo"

#check it all worked:
length(trees_1000_all)
class(trees_1000_all)

# prune tree tips to all species in this dataset only (38 species)
species_keep <- rownames(df.mam)

length(intersect(species_keep, trees_1000_all[[1]]$tip.label))

trees_1000_trimmed <- lapply(trees_1000_all, function(tr) {
  drop.tip(tr, setdiff(tr$tip.label, species_keep))
})

class(trees_1000_trimmed) <- "multiPhylo"

## if reading in trees_1000_trimmed instead:
write.nexus(trees_1000_trimmed, file = "~/Dropbox/ResearchProjects/LensTransmission/Data/mammal_trees_1000.nex")
trees_1000_trimmed<-read.nexus("~/Dropbox/ResearchProjects/LensTransmission/Data/mammal_trees_1000.nex") 

#check overlap that it really trim to the 38 species on a random tree
tree489<-trees_1000_trimmed[[489]]
overlap <- name.check(tree489, df.mam)
overlap


### 3. Perform Phylogenetic Regressions ####

#### Function for running phylogenetic analysis across trees ####
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

#### Function for summarizing phylogenetic results ####
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

#pick a formula, then a data subset if needed

#### Analysis 1: Lens transmission metrics and activity pattern #####
my_formula <- logit_Lens ~ AP.3 #Diurnal as reference
my_formula <- logit_Lens ~ AP.2 #Nocturnal as reference
my_formula <- logit_Lens ~ AP.4 #diurnal vs. non-diurnal

my_formula <- X50.T.nm ~ AP.3 #Diurnal as reference
my_formula <- X50.T.nm ~ AP.2 #Nocturnal as reference
my_formula <- X50.T.nm ~ AP.4 #diurnal vs. non-diurnal

#Set variables and dataset
vars <- all.vars(my_formula)
data.df <- df.mam[complete.cases(df.mam[, vars]),] 

#Run model
lambda_results <- run_phylogenetic_analysis(
  trees = trees_1000_trimmed,
  formula = my_formula,
  data = data.df, 
  model_type = "lambda",
  lower_bound = 0,
  upper_bound = 1
)

lambda_summary <- summarize_phylogenetic_results(lambda_results)
print(lambda_summary)

#### Analysis 2: Lens transmission metrics and visual acuity #####
# note, there will be warnings to let you know that the model is dropping taxa with missing data

#pick a formula and data subset
my_formula <- logit_Lens ~ log_VA+AP.4+Acuity.Type #acuity, acuity type, diurnal v. non-diurnal
my_formula <- logit_Lens ~ log_VA+Acuity.Type 

my_formula <- X50.T.nm ~ log_VA+AP.4+Acuity.Type #acuity, acuity type, diurnal v. non-diurnal
my_formula <- X50.T.nm ~ log_VA+Acuity.Type #acuity, acuity type for non-diurnal species subset

#subset data to match formula
#all taxa
vars <- all.vars(my_formula)
data.df <- df.mam[complete.cases(df.mam[, vars]),] 

#excluding haplorhines
nohap<-subset(df.mam,Haplorhine!="Y")
data.df <- nohap[complete.cases(nohap[, vars]),] 

#excluding diurnal taxa
dark<-subset(df.mam,AP!="D")
data.df <- dark[complete.cases(dark[, vars]),] 

# Run model    
lambda_results <- run_phylogenetic_analysis(
  trees = trees_1000_trimmed,
  formula = my_formula,
  data = data.df, 
  model_type = "lambda",
  lower_bound = 0,
  upper_bound = 1
)

lambda_summary <- summarize_phylogenetic_results(lambda_results)
print(lambda_summary)

#### Analysis 3: Lens transmission metrics and relative cornea size #####

#pick a formula and data subset
my_formula <- log_VA ~ Rel.Cornea+Acuity.Type #relative cornea size and acuity analysis

my_formula <- logit_Lens ~ Rel.Cornea+AP.4 #relative cornea size, diurnal v. non-diurnal
my_formula <- logit_Lens ~ Rel.Cornea #relative cornea size for non-diurnal species subset

my_formula <- X50.T.nm ~ Rel.Cornea+AP.4 #relative cornea size, diurnal v. non-diurnal
my_formula <- X50.T.nm ~ Rel.Cornea#relative cornea size for non-diurnal species subset

#subset data to match formula
#all taxa
vars <- all.vars(my_formula)
data.df <- df.mam[complete.cases(df.mam[, vars]),] 

#excluding haplorhines
nohap<-subset(df.mam,Haplorhine!="Y")
data.df <- nohap[complete.cases(nohap[, vars]),] 

#excluding diurnal taxa
dark<-subset(df.mam,AP!="D")
data.df <- dark[complete.cases(dark[, vars]),] 

# Run model    
lambda_results <- run_phylogenetic_analysis(
  trees = trees_1000_trimmed,
  formula = my_formula,
  data = data.df, 
  model_type = "lambda",
  lower_bound = 0,
  upper_bound = 1
)

lambda_summary <- summarize_phylogenetic_results(lambda_results)
print(lambda_summary)

#### Check Collinearity ####
# Check for multicollinearity
library(car)

# Run regular lm to check VIF (phylogenetic VIF is complex)
regular_lm <- lm(logit_Lens ~ log_VA+Acuity.Type + AP.4, 
                 data = data.df)

vif_results <- vif(regular_lm)
print(vif_results)  # VIF > 5-10 indicates problematic multicollinearity

###### Check Sample Tree for Residuals ######
model_tree<- trees_1000_trimmed[[420]] #pick one

my_formula <- logit_Lens ~ AP.3
my_formula <- logit_Lens ~ log_VA+Acuity.Type+AP.4
my_formula <- X50.T.nm ~ AP.3
my_formula <- X50.T.nm ~ Rel.Cornea

vars <- all.vars(my_formula)
data.df <- df.mam[complete.cases(df.mam[, vars]),] 
model_main <- phylolm(my_formula, phy = model_tree, data = data.df,model="lambda")
summary(model_main)

# Check residuals for normality/homoscedasticity
par(mfrow = c(1,2))
plot(model_main$fitted.values, residuals(model_main),
     xlab = "Fitted", ylab = "Residuals")
qqnorm(residuals(model_main)); qqline(residuals(model_main))


### 4. Making Plots #####
library(nationalparkcolors)
library(ggplot2)

###### Scatterplots ######

#Fig. 1b: Log visual acuity by proportion UV lens transmission, with activity pattern and haplorhine primates coded
p1<-ggplot(data=df.mam, aes(x=(log_VA), y=logit_Lens,col=AP.2,shape=Haplorhine)) + 
  geom_point(size=2.5)+theme_classic(base_size = 12)
p1  +ylab("Logit-Transformed UV Transmission")+xlab("Log-Transformed Visual Acuity (cpd)")+
  scale_color_manual(values = park_palette("Everglades"),labels = c("Nocturnal", "Cathemeral","Diurnal"))+ 
  labs(color = "Activity Pattern")

#Fig. 1c: relative cornea size X proportion UV lens transmission, with activity pattern and haplorhine coded
p1<-ggplot(data=df.mam, aes(x=(Rel.Cornea), y=logit_Lens,col=AP.2,shape=Haplorhine)) + 
  geom_point(size=2)+theme_classic(base_size = 12)
p1  +ylab("Logit-Transformed UV Transmission")+xlab("Relative Cornea Size")+ scale_color_manual(values = park_palette("Everglades"),labels = c("Nocturnal", "Cathemeral","Diurnal"))+ labs(color = "Activity Pattern")

#Fig. 1e: Log visual acuity and 50.T.nm, with activity pattern and haplorhine
p1<-ggplot(data=df.mam, aes(x=log_VA, y=X50.T.nm,col=AP.2,shape=Haplorhine)) + 
  geom_point(size=2)+theme_classic(base_size = 12)
p1  +ylab("λT0.5")+xlab("Log-Transformed Visual Acuity (cpd)")+ scale_color_manual(values = park_palette("Everglades"),labels = c("Nocturnal", "Cathemeral","Diurnal"))+ labs(color = "Activity Pattern")

#Fig. 1f: Relative cornea size and 50.T.nm, with activity pattern and haplorhine
p1<-ggplot(data=df.mam, aes(x=(Rel.Cornea), y=X50.T.nm,col=AP.2,shape=Haplorhine)) + 
  geom_point(size=2)+theme_classic(base_size = 12)
p1  +ylab("λT0.5")+xlab("Relative Cornea Size")+ scale_color_manual(values = park_palette("Everglades"),labels = c("Nocturnal", "Cathemeral","Diurnal"))+ labs(color = "Activity Pattern")

#Supplemental Fig S1: Relative cornea size and visual acuity by acuity type
p1<-ggplot(data=df.mam, aes(x=(Rel.Cornea), y=log_VA,col=AP.2,shape=Acuity.Type)) + 
  geom_point(size=2)+theme_classic(base_size = 12)
p1  +ylab("Relative Cornea Size")+xlab("Log-Transformed Visual Acuity (cpd)")+
  scale_color_manual(values = park_palette("Everglades"),labels = c("Nocturnal", "Cathemeral","Diurnal"))+ 
  labs(color = "Activity Pattern")

##### Boxplots #####

#Fig 1a: Logit transformed UV transmission by activity pattern
p1<-ggplot(data=df.mam, aes(x=AP.2, y=logit_Lens,fill=AP.2)) + 
  geom_boxplot()+theme_classic(base_size = 12)
p1  +ylab("Logit-Transformed UV Transmission")+xlab("Activity Pattern")+scale_fill_manual(values=park_palette("Everglades"))+
  scale_x_discrete(labels=c("Nocturnal","Cathemeral","Diurnal"))+xlab("")+theme(legend.position = "none")

#Fig 1d: λT0.5 and activity pattern
p1<-ggplot(data=df.mam, aes(x=AP.2, y=X50.T.nm,fill=AP.2)) + 
  geom_boxplot()+theme_classic(base_size = 12)
p1  +ylab("λT0.5")+xlab("Activity Pattern")+scale_fill_manual(values=park_palette("Everglades"))+
  scale_x_discrete(labels=c("Nocturnal","Cathemeral","Diurnal"))+xlab("")+theme(legend.position = "none")

