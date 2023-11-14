#### 4 WAY BLOCK RELAXATION ALGORITHM

# load libraries
library(tidyverse)
library(BEDASSLE)
library(dplyr)
library(RColorBrewer)
library(summixdr)

set.seed(316)

dat = data.frame(read.delim("/Users/rileylamont/Desktop/Unknown Ancestry Simulations/Data Sets/chr22contacan.txt"))

###  sample 100k SNPs for simulations
###  set parameters - what your unknown group is, reference groups are, and their props

ssample = dat %>%
  sample_n(100000)

## simulate reference ancestry groups functions

simulate_reference <- function(anc_list, AncFrame){
  N_EUR = 500
  N_AFR = 500
  N_IAM = 500
  N_SAS = 500
  N_EAS = 500
  N_FIN = 500
  ancvec = AncFrame
  ancnum = length(ancvec)
  sample_N <- c(paste(rep("N_", each = length(anc_list)), anc_list, sep = ""))
  sample_counts <- as.data.frame(lapply(sample_N, get))

  refdat = ssample %>%
    select(CHROM, POS, REF, ALT, paste(rep("AF_", each = length(anc_list)), anc_list, sep=""))
  names(refdat) <- gsub("AF_", "", names(refdat))


  #loop
  refsims = as.data.frame(matrix(0, nrow = nrow(refdat), ncol = length(anc_list)))
  for (k in 1:ancnum){
    refsimcount = t(sapply(refdat[[as.character(ancvec[k])]], function(x){x2<-as.numeric(x); rmultinom(1, as.numeric(sample_counts[k]), prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
    refsim = (2 * refsimcount[,1] + refsimcount[,2]) / (2 * as.numeric(sample_counts[k]))
    refsims[,k] = refsim
  }

  refdat_sim <- as.data.frame(cbind(refdat[,1], refdat[,2], refdat[,3], refdat[,4], refsims))
  names(refdat_sim) <- c("CHROM", "POS", "REF", "ALT", anc_list)
  refdat_sim_flt <- refdat_sim %>% filter_at(colnames(refdat_sim), any_vars(. > .01 & .< .99))
  refdat_sim_flt <<- refdat_sim_flt
}


popu_name = "EAS"
popu_N = 30000
popu_prop = popu_N/100000

pop1_name = "SAS"
pop1_N = 10000
pop1_prop = pop1_N/100000


pop2_name = "AFR"
pop2_N = 20000
pop2_prop = pop2_N/100000

pop3_name = "IAM"
pop3_N = 40000
pop3_prop = pop3_N/100000

### simulate observed group (simulate each indivudally, then combine)

popu_geno = t(sapply(ssample[,paste("AF_", popu_name, sep="")], function(x){x2<-as.numeric(x); rmultinom(1, popu_N, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
popu_dat = data.frame(AC = 2*popu_geno[,1] + popu_geno[,2],
                      AN = 2*popu_N) %>%
  mutate(AF = AC / AN)

pop1_geno = t(sapply(ssample[,paste("AF_", pop1_name, sep="")], function(x){x2<-as.numeric(x); rmultinom(1, pop1_N, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
pop1_dat = data.frame(AC = 2*pop1_geno[,1] + pop1_geno[,2],
                      AN = 2*pop1_N) %>%
  mutate(AF = AC / AN)

pop2_geno = t(sapply(ssample[,paste("AF_", pop2_name, sep="")], function(x){x2<-as.numeric(x); rmultinom(1, pop2_N, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
pop2_dat = data.frame(AC = 2*pop2_geno[,1] + pop2_geno[,2],
                      AN = 2*pop2_N) %>%
  mutate(AF = AC / AN)

pop3_geno = t(sapply(ssample[,paste("AF_", pop3_name, sep="")], function(x){x2<-as.numeric(x); rmultinom(1, pop3_N, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
pop3_dat = data.frame(AC = 2*pop3_geno[,1] + pop3_geno[,2],
                      AN = 2*pop3_N) %>%
  mutate(AF = AC / AN)

combsubset = data.frame(AC = pop1_dat$AC + pop2_dat$AC + pop3_dat$AC + popu_dat$AC,
                        AN = pop1_dat$AN + pop2_dat$AN + pop3_dat$AN + popu_dat$AN) %>%
  mutate(AF = AC / AN)


### create df with reference groups and observed pop
sumframe = data.frame(POS = ssample$POS,
                      REF = ssample$REF,
                      ALT = ssample$ALT,
                      pop1_ref = ssample[,paste("AF_", pop1_name, sep="")],
                      pop2_ref = ssample[,paste("AF_", pop2_name, sep="")],
                      pop3_ref = ssample[,paste("AF_", pop3_name, sep="")],
                      obs = combsubset$AF)

## simulate reference for all continental reference groups (does not change w/ different observed pop)
anc_list <- c("AFR", "EUR", "IAM", "SAS", "EAS", "FIN")
AncFrame <- c('AFR', 'EUR', 'IAM', 'SAS', 'EAS', 'FIN')

simulate_reference(anc_list, AncFrame)

# Pull reference data
ref_dat = refdat_sim_flt %>%
  select(POS, REF, ALT, all_of(anc_list))

# Set observed data
obs_dat = sumframe %>%
  select(POS, REF, ALT, obs)

# Merge reference and observed
mergeframe = merge(ref_dat, obs_dat, by = c("POS","REF","ALT")) %>%
  select(POS, REF, ALT, all_of(anc_list), obs)
mergeframe[,4:ncol(mergeframe)] <- sapply(mergeframe[,4:ncol(mergeframe)],as.numeric)
mergeframe <- mergeframe %>% filter_at(colnames(mergeframe), any_vars(. > .01 & .< .99))

# Run initial Summix!
# r_sum = summix(mergeframe, reference = c(paste("", pop1_name, sep=""), paste("", pop2_name, sep=""), paste("", pop3_name, sep="")), observed = "obs")

## find iniital guess for pi_u based on loss
# blr_init_guess(r_sum)
#
# pi_1_init <- r_sum[paste("", pop1_name, sep="")]*(1 - pi_u_init)
# pi_2_init <- r_sum[paste("", pop2_name, sep="")]*(1 - pi_u_init)
# pi_3_init <- r_sum[paste("", pop3_name, sep="")]*(1 - pi_u_init)
#
# AFu_init = (mergeframe$obs
#             - pull(pi_1_init)*(mergeframe[paste("", pop1_name, sep="")])
#             - pull(pi_2_init)*mergeframe[paste("", pop2_name, sep="")]
#             - pull(pi_3_init)*mergeframe[paste("", pop3_name, sep="")])/pi_u_init
#
# AFu_init[AFu_init > 1] = 1
# AFu_init[AFu_init < 0] = 0
#
# mergeframe$unk <- AFu_init

mergeframe = mergeframe %>%
  mutate(unk = runif(nrow(mergeframe)))
mergeframe$unk[mergeframe$unk > 1] = 1
mergeframe$unk[mergeframe$unk < 0] = 0

### initialize other variables for BLRA
threshold = .001 #CHANGE
diff_unknown = 1 #initialize so loop will run
emiter = 0
max_iter = 500

columns = c(paste("", pop1_name, sep=""), paste("", pop2_name, sep=""), paste("", pop3_name, sep=""), "unk", "objective", "iter")
finalframe_small = data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(finalframe_small) = columns
finalframe_small[1,] = 1

while (diff_unknown >= threshold){
  # count iterations
  emiter = emiter + 1
  # summix
  summix_br = summix(mergeframe, reference = c(paste("", pop1_name, sep=""), paste("", pop2_name, sep=""), paste("", pop3_name, sep=""), "unk"), observed ="obs", override_removeSmallAnc = TRUE, scaled_objective = FALSE)

  if (pull(summix_br[8]) == 0){
    summix_br[8] <- .1
  }

  if (emiter > max_iter){
    break
  }
  # update AFu
  update_AFu = (mergeframe$obs - pull(summix_br[5])*(mergeframe[paste("", pop1_name, sep="")])
                - pull(summix_br[6])*(mergeframe[paste("", pop2_name, sep="")])
                - pull(summix_br[7])*(mergeframe[paste("", pop3_name, sep= "")]))/pull(summix_br[8])

  update_AFu[update_AFu > 1] = 1
  update_AFu[update_AFu < 0] = 0

  mergeframe$unk = update_AFu

  # save iteration output
  outframe = data.frame(pop1 = summix_br[5],
                        pop2 = summix_br[6],
                        pop3 = summix_br[7],
                        unk = summix_br[8],
                        objective = summix_br[1],
                        iter = emiter)
  colnames(outframe) = c(paste("", pop1_name, sep=""), paste("", pop2_name, sep=""), paste("", pop3_name, sep=""), "unk", "objective", "iter")
  # bind to output frame
  finalframe_small = rbind(finalframe_small, outframe)
  diff_unknown = abs(finalframe_small$unk[emiter+1] - finalframe_small$unk[emiter])
}

finalframe_small <- finalframe_small[-1,]
tail(finalframe_small)

pal = c('AFR' = "#FDE725FF", 'EAS' = "#5DC863FF", 'EUR' = "#21908CFF", 'IAM' = "#3B528BFF", 'SAS' = "#440154FF")

ggplot(finalframe_small) +
  geom_line(aes(iter, EAS), color =  "#5DC863FF")+
  geom_line(aes(iter, EUR), color = "#21908CFF")+
  geom_line(aes(iter, SAS), color ="#440154FF")+
  geom_line(aes(iter, unk), color = "#3B528BFF") +
  geom_hline(aes(yintercept = pop1_N/100000), color = "#5DC863FF", linetype = "dashed")+
  geom_hline(aes(yintercept = pop2_N/100000), color = "#21908CFF", linetype = "dashed")+
  geom_hline(aes(yintercept = pop3_N/100000), color = "#440154FF", linetype = "dashed") +
  geom_hline(aes(yintercept = popu_N/100000), color = "#3B528BFF", linetype = "dashed") +
  geom_hline(aes(yintercept = popu_N/100000 + 0.05), color = "#3B528BFF", alpha = .2)
  labs(x = "Iterations",
       y = "Ancestry Proportion") +
  theme_bw()



cumulative_dat <- data.frame(read.csv("/Users/rileylamont/Downloads/cumulative_dat.csv"))

cumulative_dat$Unknown <- as.factor(cumulative_dat$Unknown)

ggplot(cumulative_dat, aes(x=Unknown, y=Unk_Diff)) +
  geom_boxplot(fill = pal, alpha = 0.5) + theme_minimal() + labs(x = "Hidden Ancestry Group", y = "Difference Between Actual and Estimated Proportion")+
  geom_jitter(width = .1, size = 1, alpha = .5)


ggplot(cumulative_dat, aes(x=Unknown, y=Unk_Diff)) +
  geom_violin(
    aes(fill = Unknown, fill = after_scale(colorspace::lighten(fill, .5))),
    size = 1.2, bw = .2, color = NA
  ) +
  geom_boxplot(
    width = .1, size = 1.2, outlier.shape = NA
  ) +
  stat_summary(
    geom = "point",
    fun = median,
    color = "white",
    size = 5
  )
