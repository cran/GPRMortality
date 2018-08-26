#' @export

GPRMortalitySummary <- function(model,percentile=c(0.025,0.5,0.975)){

if(is.list(model)){
          d <- model$simulation
          out <- list()
          head(d)
          location <- unique( d[,1] )
          i=0
          t = matrix(NA,0,length(percentile)+1)
          for( i in location){
            t <- rbind(t,
              cbind(i,
                      t(apply(d[d[,1]==i,-c(1,2)], 2,function(x) quantile(x,percentile) )))
            )

          }
          colnames(t)<-c("location",paste0(percentile*100,"%"))
          out$simulation <- t



          d <- model$variance
          head(d)
          t = list()
          name = names(d)
          i=1
          for( i in 1:length(d) ){
            t[[ name[i] ]] <-  t(apply(d[[i]],2,function(x) quantile(x,percentile)  ))
          }
          out$variance <- t
}

if(is.matrix(model)){
          d <- model

          head(d)
          i=1 ; ii = 0 ; iii=3

          out = matrix(NA,0,length(percentile)+3)

          location <- unique(d[,2])
          sex <- unique(d[,3])
          age_cat <- unique(d[,4])

          for( i in location){
          for( ii in sex){
          for( iii in age_cat){
            out <- rbind(out,
                         cbind(i,ii,iii,
                               t(apply(d[  d[,2]==i & d[,3]==ii & d[,4]==iii,-c(1,2,3,4)], 2,function(x) quantile(x,percentile) )))
            )

          }
          }
          }
          out = cbind( as.numeric(rownames(out)),out)
          colnames(out)<-c("year","location","sex","age_cat",paste0(percentile*100,"%"))
          rownames(out) <- NULL
}
    return(out)
}





