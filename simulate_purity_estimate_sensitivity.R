#FUNCTION
mean_no_outliers <- function(x, remove_na = TRUE) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = remove_na)
  val <- 1.5 * IQR(x, na.rm = remove_na)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  z = mean(y, na.rm = remove_na)
  return(z)
}

#RANGE OF TESTED PARAMETER VALUES
subs_range <- seq(50, 400, 50)
depth_range <- seq(10, 200, 10)
tumour_inf_range <- seq(0, 0.1, 0.01)

#1000 simulations
nRand = 1000

#Create array
res_arr = array(dim = c(length(subs_range), length(depth_range), length(tumour_inf_range), nRand), 
                dimnames = list(n_subs = subs_range, depth = depth_range, purity = tumour_inf_range, sim = 1:nRand))

#RUN
set.seed(42)
for(n_subs in subs_range){
  
  print(n_subs)
  
  for(depth in depth_range){
    
    for(purity in tumour_inf_range){
      
      for(sim in 1:nRand){
        
        #generate simulated data
        NR = rpois(n_subs, depth)
        NV = sapply(NR, function(x) rbinom(1, x, purity / 2)) #create situation like PD50297 and PD51123
        
        NR[NR == 0] = 1 #should 0 be generated for the distribution
        
        VAF = NV / NR
        
        res_arr[paste0(n_subs), paste0(depth), paste0(purity), paste0(sim)] = mean_no_outliers(VAF) * 2
        
      }
      
    }
    
  }
  
}

saveRDS(res_arr, "tumour_infiltration_sensitivity_sim.rds")