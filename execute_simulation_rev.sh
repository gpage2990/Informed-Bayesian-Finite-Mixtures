#!/usr/bin/env Rscript                                                                                                                                                                                                                                                                                                                                                                                                                                         



jobnum <- 1
for (nobs  in c(100) ) {
  for (K in c(25) ) {
    for(U in c(4) ){
      for( datatype in c(1) ){
        for(ndx in c(7, 10, 20)){
          for(Akap in c(0.1, 0.5, 1.0)){
            cat(sprintf("nohup ./simulationStudy_asym_multi.R 100 %s %s %s %s %s %s> out%s.txt \n", nobs, K, U, datatype, ndx, Akap, jobnum))
            jobnum <- jobnum + 1
          }    
        }        
      }
    }
  }
}