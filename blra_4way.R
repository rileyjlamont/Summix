#### 4 WAY BLOCK RELAXATION ALGORITHM

# load libraries
library(tidyverse)
library(BEDASSLE)
library(dplyr)
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

simulate_observed <- function(obs_grps, obs_pops){
  popu_name <<- obs_grps[1]
  popu_N <<-obs_pops[1]
  
  pop1_name <<-obs_grps[2]
  pop1_N <<- obs_pops[2]
  
  pop2_name <<- obs_grps[3]
  pop2_N <<- obs_pops[3]
  
  pop3_name <<- obs_grps[4]
  pop3_N <<- obs_pops[4]
  
  ssample <<- dat %>%
    sample_n(100000)
  
  ### simulate observed genotype (simulate each individually, then combine)
  
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
  
  combsubset <<- data.frame(AC = pop1_dat$AC + pop2_dat$AC + pop3_dat$AC + popu_dat$AC,
                          AN = pop1_dat$AN + pop2_dat$AN + pop3_dat$AN + popu_dat$AN) %>%
    mutate(AF = AC / AN)
  
  
  ### create df with reference groups and observed pop
  sumframe <<- data.frame(POS = ssample$POS,
                        REF = ssample$REF,
                        ALT = ssample$ALT,
                        pop1_ref = ssample[,paste("AF_", pop1_name, sep="")],
                        pop2_ref = ssample[,paste("AF_", pop2_name, sep="")],
                        pop3_ref = ssample[,paste("AF_", pop3_name, sep="")],
                        obs = combsubset$AF)
}

## empty df for 100x iterations

save_frame = data.frame(matrix(nrow = 100, ncol = 6))

for (i in 1:100){
  ###  sample 100k SNPs for simulations
  obs_grps <- c("IAM", "EUR", "AFR", "EAS")
  obs_pops <- c(18000, 22000, 28000, 32000)
  
  simulate_observed(obs_grps, obs_pops)

  ## simulate reference for all continental reference groups (this does not change w/ different observed pop)
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

  ## method 1 to initiate block relaxation algorithm - commented out because we will use method 2, but kept in code

  ## Run initial Summix!
  ## r_sum = summix(mergeframe, reference = c(paste("", pop1_name, sep=""), paste("", pop2_name, sep=""), paste("", pop3_name, sep="")), observed = "obs")

  ## find initial guess for pi_u based on loss
  ## blr_init_guess(r_sum)
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

  ## method 2 for initializing - randomly simulating unknown AFs to run Summix

  mergeframe = mergeframe %>%
    mutate(unk = runif(nrow(mergeframe)))
  mergeframe$unk[mergeframe$unk > 1] = 1
  mergeframe$unk[mergeframe$unk < 0] = 0

  ### initialize other variables for BLRA
  threshold = .001  #cutoff
  diff_unknown = 1  #initialize so loop will run
  emiter = 0        #iterations
  max_iter = 500    #max iterations just in case

  ### save iteration output
  columns = c(paste("", pop1_name, sep=""), paste("", pop2_name, sep=""), paste("", pop3_name, sep=""), "unk", "objective", "iter")
  finalframe_small = data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(finalframe_small) = columns
  finalframe_small[1,] = 1

  while (diff_unknown >= threshold){
    # count iterations
    emiter = emiter + 1
    # summix
    summix_br = summix(mergeframe, reference = c(paste("", pop1_name, sep=""), paste("", pop2_name, sep=""), paste("", pop3_name, sep=""), "unk"), observed ="obs", override_removeSmallAnc = TRUE, scaled_objective = FALSE)

    # can't divide by 0 - just in case
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

  save_frame[i, ] = finalframe_small[emiter, ]
}

colnames(save_frame) = c(paste("", pop1_name, sep=""), paste("", pop2_name, sep=""), paste("", pop3_name, sep=""), "unk", "objective", "iter")
csv_name <-paste("2023127_PopU_", popu_name, "_", popu_N, "_Ref_", pop1_name, "_", pop1_N, "_",
                 pop2_name, "_", pop2_N, "_", pop3_name, "_", pop3_N, sep="")

write.csv(save_frame, file=paste("~/Desktop/BLRA_Sims/", csv_name, ".csv", sep="")) #save to desktop

# pal = c('AFR' = "#FDE725FF", 'EAS' = "#5DC863FF", 'EUR' = "#21908CFF", 'IAM' = "#3B528BFF", 'SAS' = "#440154FF")

ggplot(finalframe_small) +
  geom_line(aes(iter, EAS), color =  "#5DC863FF")+
  geom_line(aes(iter, EUR), color = "#21908CFF")+
  geom_line(aes(iter, SAS), color ="#440154FF")+
  geom_line(aes(iter, unk), color = "#3B528BFF") +
  geom_hline(aes(yintercept = pop1_N/100000), color = "#5DC863FF", linetype = "dashed")+
  geom_hline(aes(yintercept = pop2_N/100000), color = "#21908CFF", linetype = "dashed")+
  geom_hline(aes(yintercept = pop3_N/100000), color = "#440154FF", linetype = "dashed") +
  geom_hline(aes(yintercept = popu_N/100000), color = "#3B528BFF", linetype = "dashed") +
  labs(x = "Iterations",
       y = "Ancestry Proportion") +
  theme_bw()

