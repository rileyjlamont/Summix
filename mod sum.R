library(nloptr)
library(Summix)
library(dplyr)

mod_summix_calc <- function (data, reference, observed, pi.start = c()) 
{
  if (!is(object = data, class2 = "data.frame")) {
    stop("ERROR: data must be a data.frame as described in the vignette")
  }
  if (typeof(observed) != "character") {
    stop("ERROR: 'observed' must be a character string for the column name of the observed ancestry in data")
  }
  if (!(observed %in% names(data))) {
    stop("ERROR: 'observed' must be the column name of the observed ancestry in data")
  }
  if (typeof(reference) != "character") {
    stop("ERROR: 'reference' must be a vector of column names from data to be used in the reference")
  }
  if (all(reference %in% names(data)) == FALSE) {
    stop("ERROR: 'reference' must be a vector of column names from data to be used in the reference")
  }
  filteredNA <- length(which(is.na(data[, observed] == TRUE)))
  observed.b <- as.data.frame(data[which(is.na(data[, observed]) == 
                                           FALSE), observed])
  refmatrix <- as.data.frame(data[which(is.na(data[, observed]) == 
                                          FALSE), reference])
  if (length(pi.start) != 0) {
    if (is.numeric(pi.start) == FALSE) {
      stop("ERROR: Please make sure pi.start is a positive numeric vector of length reference that sums to one")
    }
    if (length(pi.start) != length(reference)) {
      stop("ERROR: Please make sure pi.start is a positive numeric vector of length reference that sums to one")
    }
    if (all(pi.start > 0) == FALSE) {
      stop("ERROR: Please make sure pi.start is a positive numeric vector of length reference that sums to one")
    }
    if (sum(pi.start) != 1) {
      stop("ERROR: Please make sure pi.start is a positive numeric vector of length reference that sums to one")
    }
    starting = pi.start
  }
  else {
    starting = rep(1/ncol(refmatrix), ncol(refmatrix))
  }
  fn.ancmix = function(x) {
    expected = x %*% t(as.matrix(refmatrix))
    minfunc = sum((expected - observed.b)^2)
    return(minfunc)
  }
  gr.ancmix <- function(x) {
    gradfunc = x %*% t(as.matrix(refmatrix)) - observed.b
    gradvec <- apply(2 * refmatrix * t(gradfunc), 2, sum)
    return(gradvec)
  }
  
  #REMOVE HEQ - need equality <= 1 
  #heq.ancmix = function(x) {
    #equality = sum(x)
    #return(equality - 1)
  #}
  
  #hin returns hin>= 0 ... need equality = sum(x) -> 1 - equality 
  hin.ancmix <- function(x) {
    equality = 1 - sum(x) 
    h = c(x, equality)
    return(h)
  }
  start_time = Sys.time()
  S = suppressMessages(slsqp(starting, fn = fn.ancmix, gr = gr.ancmix, 
                             hin = hin.ancmix))
  end_time = Sys.time()
  ttime = end_time - start_time
  d <- data.frame(matrix(ncol = length(reference) + 4, nrow = 1))
  colnames(d) <- c("objective", "iterations", "time", "filtered", 
                   colnames(refmatrix))

  
  d[1] <- S$value
  d[2] <- S$iter
  d[3] <- ttime
  d[4] <- filteredNA
  d[5:(length(reference) + 4)] <- S$par[1:length(reference)]
  return(d)
}

calc_scaledObj_mod <- function(data, observed, reference) {
  start_time <- Sys.time()
  bins <- c("0-0.1", "0.1-0.3", "0.3-0.5")
  multiplier <- c(5, 1.5, 1)
  
  data$obs <- data[,observed]
  
  # create bins and run summix on each bin
  data <- data %>% rowwise() %>%
    mutate(bin = case_when(obs <= 0.1 ~ 1,
                           obs >0.1 & obs <= 0.3 ~ 2,
                           obs >0.3 & obs <= 0.5 ~ 3))
  for(b in 1:3) {
    # subset data to only SNPs in given bin
    subdata <- data[which(data$bin == b),]
    if(nrow(subdata) > 0) {
      res <- mod_summix_calc(subdata, observed = observed, reference = reference)
      res$obj_adj <-  res$objective/(nrow(subdata)/1000)
      res$nSNPs <- nrow(subdata)
      res$bin <- bins[b]
      if(b == 1) {
        sum_res <- res
      } else {
        sum_res <- rbind(sum_res, res)
      }
    } else {
      res <- mod_summix_calc(data, observed = observed, 
                            reference = reference)
      res$obj_adj <- NA
      res$nSNPs <- 0
      res$bin <- bins[b]
      if(b == 1) {
        sum_res <- res
      } else {
        sum_res <- rbind(sum_res, res)
      }
    }
  }
  objective_scaled <- 0
  nonNA <- 0
  for(b in 1:3) {
    if(!is.na(sum_res[b, ]$obj_adj)) {
      nonNA <- nonNA + 1
      objective_scaled <- objective_scaled + 
        (sum_res[b, ]$obj_adj*multiplier[b])
    }
  }
  objective_scaled <- objective_scaled/nonNA
  end_time <- Sys.time()
  #print(difftime(end_time, start_time, units = "auto"))
  
  return(objective_scaled)
}

mod_summix <- function(data, reference, observed, pi.start = NA) {
    if(!is.na(pi.start)) {
      sum_res <- mod_summix_calc(data = data, reference = reference, 
                                observed = observed, pi.start = pi.start)
    } else {
      sum_res <- mod_summix_calc(data = data, reference = reference, observed = observed)
    }
    new_obj <- calc_scaledObj_mod(data = data, observed = observed, reference = reference)
    sum_res[1] <- new_obj
    if(new_obj >= 0.5 & new_obj < 1.5) {
      print(paste0("CAUTION: Objective/1000SNPs = ", round(new_obj, 4), 
                   " which is within the moderate fit range"))
    } else if (new_obj >= 1.5) {
      print(paste0("WARNING: Objective/1000SNPs = ", round(new_obj, 4), 
                   " which is above the poor fit threshold"))
    }
    return(sum_res)
  }



#test this out with a loop!
#empty data frame
frame_mod = data.frame(matrix(nrow = 100, ncol = 6))
frame_test = data.frame(matrix(nrow = 100, ncol = 6))

#macrovariables
N = 100000 #total sample size
pop1_Nloop = c(1:7)*2500 #first fixed prop, adjust based on prop_missing
prop_missing = .9 #can also loop over this, but test run


for(h in 1:length(pop1_Nloop)){
for(i in 1:100){

ssample = subframe %>% 
  sample_n(100000)

#create observed European, 45%
eur_name = "AF_EUR"
eur_pop = pop1_Nloop[h]

eur_geno = t(sapply(ssample[[eur_name]], function(x){x2<-as.numeric(x); rmultinom(1, eur_pop, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
eur_dat = data.frame(AC = 2*eur_geno[,1] + eur_geno[,2],
                     AN = 2*eur_pop) %>% 
  mutate(AF = AC / AN)

#create observed EAS, 45%
afr_name = "AF_AFR"
afr_pop =  N - prop_missing*N - eur_pop

afr_geno = t(sapply(ssample[[afr_name]], function(x){x2<-as.numeric(x); rmultinom(1, afr_pop, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
afr_dat = data.frame(AC = 2*afr_geno[,1] + afr_geno[,2],
                     AN = 2*afr_pop) %>% 
  
  mutate(AF = AC / AN)


#create observed African, 10%
eas_name = "AF_EAS"
eas_pop = prop_missing*N

eas_geno = t(sapply(ssample[[eas_name]], function(x){x2<-as.numeric(x); rmultinom(1, eas_pop, prob=(c(x2**2, 2*x2*(1-x2), (1-x2)**2)))}))
eas_dat = data.frame(AC = 2*eas_geno[,1] + eas_geno[,2],
                     AN = 2*eas_pop) %>% 
  
  mutate(AF = AC / AN)

#combined observed population
combsubset = data.frame(AC = afr_dat$AC + eur_dat$AC + eas_dat$AC,
                        AN = afr_dat$AN + eur_dat$AN + eas_dat$AN) %>% 
  mutate(AF = AC / AN)

#remove AFR
sumframe2 = data.frame(POS = ssample$POS,
                       REF = ssample$REF,
                       ALT = ssample$ALT,
                       ref_eur = ssample$AF_EUR,
                       ref_eas = ssample$AF_EAS,
                       obs = combsubset$AF)


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

refdat = ssample %>% select(CHROM, POS, REF, ALT, paste(rep("AF_", each = length(anc_list)), anc_list, sep="")) 
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

R_sum_mod = mod_summix(mergeframe2, reference = c(paste("EUR"), paste("AFR")), observed = "obs") 
R_sum_test = summix(mergeframe2, reference = c(paste("EUR"), paste("AFR")), observed = "obs") 

##((commented out below, used to need to find loss/1k SNPs but new Summix loss does that for me))
##R_sum_mod$objective <- R_sum_mod$objective*(1000/(nrow(mergeframe2))) 
##R_sum_test$objective <- R_sum_test$objective*(1000/(nrow(mergeframe2)))

frame_mod[i, ] = R_sum_mod
frame_test[i, ] = R_sum_test

colnames(frame_mod) = c("Objective", "Iterations", "Time", "Filtered", "EUR", "AFR")
colnames(frame_test) = c("Objective", "Iterations", "Time", "Filtered", "EUR", "AFR")

tmp2<-paste("20221230_Mod_Sum_Pop1EUR", eur_pop, "_Pop2AFR", afr_pop, "_POPuEAS", eas_pop, sep="")
tmp3<-paste("20221230_Test_Sum_POP1EUR", eur_pop, "_Pop2AFR", afr_pop, "_POPuEAS", eas_pop, sep="")

}

write.csv(frame_mod, file=paste("~/Desktop/Unknown Ancestry Simulations/Mod_Summix_Sims/", tmp2, ".csv", sep=""))
write.csv(frame_test, file=paste("~/Desktop/Unknown Ancestry Simulations/Mod_Summix_Sims/", tmp3, ".csv", sep=""))
}

