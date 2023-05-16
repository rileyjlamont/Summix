library(Summix)
library(BEDASSLE)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)
library(combinat)
library(gtools)
library(dplyr)
library(tidyr)
library(permute)

set.seed(2022)

chr22subset = read.delim("/Users/rileylamont/Desktop/Unknown Ancestry Simulations/Data Sets/chr22AFfilterednew.txt") 


#macrovariables
Ancestries = c('AFR', 'EAS', 'EUR') #CHANGE

#possible combos
AncDat = data.frame(permn(Ancestries))
AncDat <- AncDat[ -c(3,5,6) ]

prop_missing = .4  #fixed unknown ancestry, CHANGE
N = 100000 #total sample size
pop1_Nloop = c(1:11)*5000 #first fixed prop, CHANGE (based on prop-missing)

popu_N = N*prop_missing

subframe = data.frame(chr22subset)

#initialize empty data frames, 100 iterations
frame3 = data.frame(matrix(nrow = 100, ncol = 7))
frame2 = data.frame(matrix(nrow = 100, ncol = 6))

for(j in 1:ncol(AncDat)){
for(h in 1:length(pop1_Nloop)){
for(i in 1:100){
  
  
      ssample = subframe %>% 
        sample_n(100000)
      
      #simulate observed group 
      pop1_N = pop1_Nloop[h]
      pop1_name = AncDat[1,j] 
      
      pop1_geno = t(sapply(ssample[,paste("AF_", pop1_name, sep="")], function(x){x2<-as.numeric(x); rmultinom(1, pop1_N, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
      pop1_dat = data.frame(AC = 2*pop1_geno[,1] + pop1_geno[,2],
                            AN = 2*pop1_N) %>% 
        mutate(AF = AC / AN)
      
      #simulate observed group 2
      pop2_name = AncDat[2,j] 
      pop2_N = N - prop_missing*N - pop1_N
      
      pop2_geno = t(sapply(ssample[,paste("AF_", pop2_name, sep="")], function(x){x2<-as.numeric(x); rmultinom(1, pop2_N, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
      pop2_dat = data.frame(AC = 2*pop2_geno[,1] + pop2_geno[,2],
                            AN = 2*pop2_N) %>% 
        mutate(AF = AC / AN)
      
      
      #simulate observed UNKNOWN ancestry, 10%
      popu_name = AncDat[3,j]
      
      popu_geno = t(sapply(ssample[,paste("AF_", popu_name, sep="")], function(x){x2<-as.numeric(x); rmultinom(1, popu_N, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
      popu_dat = data.frame(AC = 2*popu_geno[,1] + popu_geno[,2],
                            AN = 2*popu_N) %>% 
        
        mutate(AF = AC / AN)
      

      
      #combine observed population
      combsubset = data.frame(AC = pop1_dat$AC + pop2_dat$AC + popu_dat$AC,
                              AN = pop1_dat$AN + pop2_dat$AN + popu_dat$AN) %>% 
        mutate(AF = AC / AN)
      
      #all together
      sumframet = data.frame(POS = ssample$POS,
                             REF = ssample$REF,
                             ALT = ssample$ALT,
                             pop1_ref = ssample[,paste("AF_", pop1_name, sep="")],
                             pop2_ref = ssample[,paste("AF_", pop2_name, sep="")],
                             popu_ref = ssample[,paste("AF_", popu_name, sep="")],
                             obs = combsubset$AF) 
      
      
      #remove Unknown
      sumframe2 = data.frame(POS = ssample$POS,
                             REF = ssample$REF,
                             ALT = ssample$ALT,
                             pop1_ref = ssample[,paste("AF_", pop1_name, sep="")],
                             pop2_ref = ssample[,paste("AF_", pop2_name, sep="")],
                             obs = combsubset$AF)
      
      
      #simulate reference - does not changed with other loops
      anc_list <- c("AFR", "EUR", "IAM", "SAS", "EAS")
      
      N_EUR = 500
      N_AFR = 500
      N_IAM = 500
      N_SAS = 500
      N_EAS = 500
      
      
      AncFrame <- c('AFR', 'EUR', 'IAM', 'SAS', 'EAS')
      ancvec = AncFrame
      ancnum = length(ancvec)
      
      ####Simulate new reference data based on N individuals per real populations in HGDP and 1KG###########
      
      #Real data sample sizes
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
      refdat_sim_flt <- refdat_sim %>% 
        filter(.99 > refdat_sim[,5:ncol(refdat_sim)] & refdat_sim[,5:ncol(refdat_sim)] > .01) 
      
      #clear
      rm(refdat_sim)
      gc()
      
      # Pull reference data
      refdatm = refdat_sim_flt %>% 
        select(POS, REF, ALT, all_of(anc_list))
      
      # Set Observed data
      obsvecm = sumframe2 %>% 
        select(POS, REF, ALT, obs)
      
      # Merge reference and observed
      mergeframe2 = merge(refdatm, obsvecm, by = c("POS","REF","ALT")) %>% 
        select(POS, REF, ALT, all_of(anc_list), obs)
      mergeframe2[,4:ncol(mergeframe2)] <- sapply(mergeframe2[,4:ncol(mergeframe2)],as.numeric)
      
      # Run Summix
      R_sum2 = summix(mergeframe2, reference = c(paste("", pop1_name, sep=""), paste("", pop2_name, sep="")), observed = "obs") 
      
      rm(mergeframe2)
      gc()
      
      
      # Set Observed data
      obsvecm = sumframet %>% 
        select(POS, REF, ALT, obs)
      
      # Merge reference and observed
      mergeframe3 = merge(refdatm, obsvecm, by = c("POS","REF","ALT")) %>% 
        select(POS, REF, ALT, all_of(anc_list), obs)
      mergeframe3[,4:ncol(mergeframe3)] <- sapply(mergeframe3[,4:ncol(mergeframe3)],as.numeric)
      
      R_sum3 = summix(mergeframe3, reference = c(paste("", pop1_name, sep=""), paste("", pop2_name, sep=""), paste("", popu_name, sep="")), observed = "obs")
      
      
      frame2[i, ] = R_sum2
      frame3[i, ] = R_sum3  
      
      colnames(frame2) = c("Objective", "Iterations", "Time", "Filtered", paste("AF_", pop1_name, sep=""), paste("AF_", pop2_name, sep=""))
      colnames(frame3) = c("Objective", "Iterations", "Time", "Filtered", paste("AF_", pop1_name, sep=""),paste("AF_", pop2_name, sep=""), paste("AF_", popu_name, sep=""))
      
      #create informative file name to keep track of outputs
tmp2<-paste("20221003Pop1", pop1_name, "_Npop1", pop1_N, "_Pop2", pop2_name, "_Npop2", pop2_N, "_PopU", popu_name, "_Npopu", popu_N, sep="")
tmp3<-paste("202201003_TEST_Pop1", pop1_name, "_Npop1", pop1_N, "_Pop2", pop2_name, "_Npop2", pop2_N, "_Pop3", popu_name, "_Npop3", popu_N, sep="")



}
  #save
  write.csv(frame2, file=paste("~/Desktop/Unknown Ancestry Simulations/AFREUREAS/", tmp2, ".csv", sep=""))
  write.csv(frame3, file=paste("~/Desktop/Unknown Ancestry Simulations/AFREUREAS/", tmp3, ".csv", sep=""))
  
  }
}




