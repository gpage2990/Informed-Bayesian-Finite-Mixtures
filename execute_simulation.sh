#!/usr/bin/env Rscript                                                                                                                                                                                                                                                                                                                                                                                                                                         


jobnum <- 1
for ( nobs  in c(1000, 100) ) {
  for ( K in c(25) ) {
    for(U in c(1,2,5,10) ){
      for( datatype in c(1,2) ){
        cat(sprintf("nohup ./simulationStudy_asym2.R 100 %s %s %s %s > out%s.txt \n", nobs, K, U, datatype, jobnum))
        jobnum <- jobnum + 1
        
      }
    }
  }
}