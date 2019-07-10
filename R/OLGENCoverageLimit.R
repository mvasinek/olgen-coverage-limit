
coverage.relax <- function(vaf, e, pfp, ptp){
  coverage = 0
  vf = 0
  
  cov_start = 30
  cov_end = 10000
  
  for(cov_aux in cov_start:cov_end){
    k = seq(0,cov_aux,by=1)
    e_binom = pbinom(k, cov_aux, e) 
    
    se = 0
    
    for(e_id in 1:cov_aux){
      if(1 - e_binom[e_id] >= pfp){
        se = se + 1
      }
      else{
        break
      }
    }
    
    
    vf_binom = pbinom(k, cov_aux, vaf)
    if(1 - vf_binom[se] >= ptp){
      vf = se
      coverage = cov_aux
      break
    }
  }
  c(coverage, vf)
}

coverage.relax.minvar <- function(vaf,e,pfp,ptp,minvar){
  coverage = 0
  vf = 0
  
  cov_start = 30
  cov_end = 10000
  
  for(cov_aux in cov_start:cov_end){
    k = seq(0,cov_aux,by=1)
    e_binom = pbinom(k, cov_aux, e) 
    
    se = 0
    for(e_id in 1:cov_aux){
      if(1 - e_binom[e_id] >= pfp){
        se = se + 1
      }
      else{
        break
      }
    }
    
    vf_binom = pbinom(k, cov_aux, vaf)
    if(1 - vf_binom[se] >= ptp && minvar <= se){
      vf = se
      coverage = cov_aux
      break
    }
  }
  c(coverage,vf)
}

#Usage
#coverage.relax(0.1,0.01,0.001,0.999)
#coverage.relax.minvar(0.1,0.01,0.001,0.999,10)