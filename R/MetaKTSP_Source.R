######################################################################################################################
######################################################################################################################
######################################################################################################################
#
# MetaKTSP: A Meta-Analytic Top Scoring Pair Method for Robust Cross-Study Validation of Omics Prediction Analysis
# Version : 1.1.1
# Authors : SungHwan Kim, Chien-Wei Lin and George C. Tseng
# latest update : 08/10/2016
#
######################################################################################################################
######################################################################################################################
######################################################################################################################
##' Fit a prediction rule of meta top scoring pair
##'
##'
##' @title A Meta-Analytic Top Scoring Pair Method for Robust Cross-Study Validation of Omics Prediction Analysis
##' @param DList Input variable matrix (a list of multiple datasets; row=features, column=samples).
##' @param Method "Fisher" = fisher's meta-analysis, "Stouffer" = Stouffer's meta-analysis, NULL = mean score method.
##' @param K.max The maximum number of top score pairs.
##' @param is.VO Logical value indicating whether K selection is Variance optimization (VO).
##' @param Para.num.gene # of top scoring pairs that are initially considered.
##' @export
##' @return meta.tspobj = a list of meta top scoring pairs results
##' meta.tspobj$model = Meta top scoring pairs (n by k-1) meta.tspobj$num.K = #
##' of K meta top scoring pairs
##' @author Chien-Wei Lin <masaki396@gmail.com>
##' @example
MetaTSP.Pvalue <- function(DList,  Method = "Fisher", K.max, is.VO=TRUE, Para.num.gene = 10000) {

  # Description
  #     Fit a prediction rule of meta top scoring pair
  # Usage
  #     S = MetaTSP(DList,  Method = "Fisher", K.max, is.VO=TRUE, Para.num.gene = 10000)
  # Input
  #     X = Input variable matrix (a list of multiple datasets; row=features, column=samples)
  #     Method = "Fisher" = fisher's meta-analysis, "Stouffer" = Stouffer's meta-analysis, NULL = mean score method
  #     K.max = # of top score pairs
  #     is.VO = Logical value indicating whether K selection is Variance optimization (VO).
  #     Para.num.gene = # of top scoring pairs that are initially considered.
  # Output
  #     meta.tspobj = a list of meta top scoring pairs results
  #        meta.tspobj$model = Meta top scoring pairs (n by k-1)
  #        meta.tspobj$num.K = # of K meta top scoring pairs

  # Computing the scores of all gene pairs
  out <- foreach( j=1:length(DList)) %dopar% {
    dat  <- t(DList[[j]])
    n <- dim(dat)[1]
    m <- dim(dat)[2]
    grp <- as.numeric(rownames(dat))
    labels <- as.character(unique(grp))
    xx <- as.matrix(cbind(rep(1,n), grp))
    hatxx <- solve(t(xx) %*% xx) %*% t(xx)
    tmp <- foreach(i = c(1:(m-1))) %dopar% {
      yy <- dat[ ,-c(1:i)] < dat[,i]
      t(t((hatxx %*% yy)[2,]))
    }
    tmp <-do.call(rbind, tmp)
    return(tmp)
  }
  out <- do.call(cbind, out)
  #column refer to different study

  #gene pair index
  m <- dim(t(DList[[1]]))[2]
  ind <- foreach(i=1:(m-1)) %dopar% {
    dat  <- t(DList[[1]])
    cbind(rep(i, length(i:m)-1 ), (i+1):m,  rep(colnames(dat)[i], length(i:m)-1), colnames(dat)[(i+1):m])
  }
  ind <- do.call(rbind, ind)

  # Standardize the score
  out <- scale(out)
  Score.pvalue <- 1 - pnorm(out, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
  rownames(Score.pvalue) = NULL

  if(Method == "Fisher")
  {
  ##  Opposite direction p-value is required.
    Score.pvalue <- rbind(Score.pvalue, 1 - Score.pvalue)
  }

  tmp.metaTSP <- switch( Method,
                         Stouffer = {
                           Score_quantile <- -qnorm(Score.pvalue, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
                           Score_stouffer <- rowSums(Score_quantile) / sqrt(ncol(Score_quantile))
                           upper.q <- quantile(Score_stouffer, probs = seq(0, 1, 0.05))[20]
                           lower.q <- quantile(Score_stouffer, probs = seq(0, 1, 0.05))[2]
                           stouffer.ind <- (Score_stouffer > upper.q) | (Score_stouffer < lower.q)
                           data.frame(ind[stouffer.ind, ], Score_stouffer = Score_stouffer[stouffer.ind], Score_quantile[stouffer.ind,], stringsAsFactors = F)
                         },
                         Fisher = {
                           tmp <-  -rowSums(log(Score.pvalue))
                           .ind <- rbind(ind, ind)
                           score.direction <- c(rep(1,nrow(ind)),rep(-1,nrow(ind)))
                           data.frame(.ind, score.direction, tmp, -log(Score.pvalue), stringsAsFactors = F)
                         }
  )

  tspobj <- list()
  if(Method =="Fisher"){
    tmp <- tmp.metaTSP[ ,6]
    o <-  order(tmp, decreasing=TRUE)
    tspobj$m.index <- tmp.metaTSP[o[1:Para.num.gene], ]
  } else {
    tspobj$m.index <- tmp.metaTSP[order(abs(as.numeric(tmp.metaTSP[ ,5])), decreasing=TRUE ), ]
  }

  #remove the genes that appeared more than one time
  #note that, even only one gene in a pair appear multiple time, the entire pair was removed

  tmp.index <- c()
  comp.tmp <- c()
  o.chs <- dim(tspobj$m.index)[1]

  i=1
  repeat{
    if(sum(comp.tmp %in% c(as.character(tspobj[[1]][i,1]),as.character(tspobj[[1]][i,2]))) == 0){
      tmp.index[i] <- sum(comp.tmp %in% c(as.character(tspobj[[1]][i,1]),as.character(tspobj[[1]][i,2]))) == 0
      comp.tmp <- c(comp.tmp, c(as.character(tspobj[[1]][i,1]),as.character(tspobj[[1]][i,2])))
    }else{
      tmp.index[i] <- sum(comp.tmp %in% c(as.character(tspobj[[1]][i,1]),as.character(tspobj[[1]][i,2]))) == 0
    }
    i <- i + 1
    if(sum(tmp.index) == K.max) {break}
  }

  tspobj$m.index <- tspobj$m.index[c(tmp.index, rep( FALSE, o.chs-length(tmp.index)))   ,]

  if(Method == "Fisher"){
    colnames(tspobj$m.index) = c("GeneIndex1", "GeneIndex2", "Gene1", "Gene2", "Direction", "Score_overall",
                                 paste("Score_study", 1:length(DList), sep = ""))

  } else {
    colnames(tspobj$m.index) = c("GeneIndex1", "GeneIndex2", "Gene1", "Gene2", "Score_overall",
                                 paste("Score_study", 1:length(DList), sep = ""))

  }
  gene.pair.table = tspobj$m.index[1:min(nrow(tspobj$m.index), K.max),]

  K <- K.max
  if(is.VO){
    k.try = 2:K.max
    Total_tau <- foreach(k =  k.try ,.combine=c) %dopar% {
      Var.0 <- var(apply(foreach(ii=1:k, .combine=rbind)%do%{
        as.numeric(foreach( j = 1:length(DList) ,.combine=c) %do% {
          label <- colnames(DList[[j]])
          DList[[j]][tspobj$m.index[ii, 3], label=="0" ]  >  DList[[j]][tspobj$m.index[ii, 4], label== "0" ]
        })},2,sum))

      Var.1 <- var(apply(foreach(ii=1:k, .combine=rbind)%do%{
        as.numeric(foreach( j = 1:length(DList) ,.combine=c) %do% {
          label <- colnames(DList[[j]])
          DList[[j]][tspobj$m.index[ii, 3], label=="1" ]  >  DList[[j]][tspobj$m.index[ii, 4], label== "1" ]
        })},2,sum))

      sum(abs(as.numeric(tspobj$m.index[1:k, 5]))) / sqrt(Var.0 + Var.1)
    }

    i <- which(Total_tau==max(Total_tau))[1]
    tspobj$m.index <- tspobj$m.index[1:k.try[i], ]
  } else {

    if(K==1){
      tspobj$m.index <- tspobj$m.index
    } else{
      ## When "is.VO=FALSE", just simple K collection.
      tspobj$m.index <- tspobj$m.index[1:K,]
    }
  }

  rownames(tspobj$m.index) <- rownames(gene.pair.table) <- NULL
  meta.tspobj <- list()
  meta.tspobj$model <- tspobj$m.index
  meta.tspobj$num.K <- dim(tspobj$m.index)[1]
  meta.tspobj$gene.pair.table <- gene.pair.table
  if(is.VO) meta.tspobj$VO <- Total_tau

  return(meta.tspobj)
}

##' MetaTSP.mean is a function of a meta-analytic top scoring pair (MetaTSP) algorithm via Mean score method.
##'
##'
##' @title A Meta-Analytic K Top Scoring Pair via Mean Score Method
##' @param DList Input variable matrix (a list of multiple datasets; row=features, column=samples).
##' @param K.max The maximum number of top score pairs.
##' @param is.VO Logical value indicating whether K selection is Variance optimization (VO)
##' @param quantile.index A quantile threshold for selecting top scoring pair. Default is 10.
##' @return An object of class list is returned with components as following:
##' model includes meta k top scoring pairs (n by k-1).
##' num.K is the number of K meta top scoring pairs.
##' @author Chien-Wei Lin
##' @export
##' @examples
MetaTSP.mean <- function(DList, K.max, is.VO=TRUE, quantile.index=10){
  # Description
  #     Fit a prediction rule of meta top scoring pair
  # Usage
  #     S = MetaTSP(DList,  Method = "Fisher", K.max, is.VO=TRUE, Para.num.gene = 10000)
  # Input
  #     X = Input variable matrix (a list of multiple datasets; row=features, column=samples)
  #     Method = "Fisher" = fisher's meta-analysis, "Stouffer" = Stouffer's meta-analysis, NULL = mean score method
  #     K.max = # of top score pairs
  #     is.VO = Logical value indicating whether K selection is Variance optimization (VO).
  #     Para.num.gene = # of top scoring pairs that are initially considered.
  # Output
  #     meta.tspobj = a list of meta top scoring pairs results
  #        meta.tspobj$model = Meta top scoring pairs (n by k-1)
  #        meta.tspobj$num.K = # of K meta top scoring pairs

  # Computing the scores of all gene pairs

  out <- foreach( j=1:length(DList)) %dopar% {
    dat  <- t(DList[[j]])
    n <- dim(dat)[1] #number of samples
    m <- dim(dat)[2] #number of features
    grp <- as.numeric(rownames(dat))
    labels <- as.character(unique(grp))
    xx <- as.matrix(cbind(rep(1,n), grp))
    hatxx <- solve(t(xx) %*% xx) %*% t(xx) #2 x n matrix
    tmp <- foreach(i = c(1:(m-1))) %dopar% {
      yy <- dat[ ,-c(1:i)] < dat[,i] # n x m matrix
      t(t((hatxx %*% yy)[2,]))
    }
    tmp <-do.call(rbind, tmp)
    return(tmp)
  }
  out <- do.call(cbind, out)

  score <- rowSums(out)

  tmp.thread <- abs(score) > quantile(abs(score), probs = seq(0, 1, 0.1))[quantile.index]
  score <- score[tmp.thread]
  out.tmp <- out[tmp.thread,]
  rownames(out.tmp) = NULL

  m <- dim(t(DList[[1]]))[2] #number of features
  ind <- foreach(i=1:(m-1)) %dopar% {
    dat  <- t(DList[[1]])
    cbind(rep(i, length(i:m)-1), (i+1):m, rep(colnames(dat)[i], length(i:m)-1), colnames(dat)[(i+1):m])
  }
  ind <- do.call(rbind, ind)
  ind <- ind[tmp.thread,]

  o.score <- order(abs(score), decreasing=T)
  o.chs = ifelse(dim(ind)[1]<1000, dim(ind)[1], 1000) #only use up to 1000 pairs
  tspobj <- list()
  tspobj$m.index <- data.frame(ind[o.score,][1:o.chs,], score[o.score][1:o.chs], out.tmp[o.score,][1:o.chs,], stringsAsFactors = F)
  #---------------------------------------------------------------------------
  # remove the genes that appeared more than one time
  # note that, even only one gene in a pair appear multiple time,
  # the entire pair was removed
  #---------------------------------------------------------------------------
  tmp.index <- c()
  comp.tmp <- c()
  i=1
  repeat{
    if(sum(comp.tmp %in% c(as.character(tspobj[[1]][i,1]),as.character(tspobj[[1]][i,2]))) == 0){
      tmp.index[i] <- sum(comp.tmp %in% c(as.character(tspobj[[1]][i,1]),as.character(tspobj[[1]][i,2]))) == 0
      comp.tmp <- c(comp.tmp, c(as.character(tspobj[[1]][i,1]),as.character(tspobj[[1]][i,2])))
    }else{
      tmp.index[i] <- sum(comp.tmp %in% c(as.character(tspobj[[1]][i,1]),as.character(tspobj[[1]][i,2]))) == 0
    }
    i <- i + 1
    if( i == dim(tspobj$m.index)[1]) {break}
  }
  tspobj$m.index <- tspobj$m.index[c(tmp.index, rep( FALSE, o.chs-length(tmp.index))) ,]
  colnames(tspobj$m.index) = c("GeneIndex1", "GeneIndex2", "Gene1", "Gene2", "Score_overall",
                               paste("Score_study", 1:length(DList), sep = ""))
  gene.pair.table = tspobj$m.index[1:min(nrow(tspobj$m.index), K.max),]
  ### Variation optimization approach ###
  if(is.VO){
    # k.try = seq(3, dim(tspobj$m.index)[1]  ,2)[seq(3, dim(tspobj$m.index)[1]  ,2) <= K.max]
    k.try = 2:K.max
    Total_tau <- foreach(k = k.try,.combine=c) %dopar% {
      Var.0 <- var(apply(foreach(ii=1:k, .combine=rbind)%do%{
        as.numeric(foreach( j = 1:length(DList) ,.combine=c) %do% {
          label <- colnames(DList[[j]])
          DList[[j]][tspobj$m.index[ii, 3], label=="0" ]  >  DList[[j]][tspobj$m.index[ii, 4], label== "0" ]
        })},2,sum))

      Var.1 <- var(apply(foreach(ii=1:k, .combine=rbind)%do%{
        as.numeric(foreach( j = 1:length(DList) ,.combine=c) %do% {
          label <- colnames(DList[[j]])
          DList[[j]][tspobj$m.index[ii, 3], label=="1" ]  >  DList[[j]][tspobj$m.index[ii, 4], label== "1" ]
        })},2,sum))

      sum(abs(as.numeric(tspobj$m.index[1:k, 5]))) / sqrt(Var.0 + Var.1)
    }

    i <- which(Total_tau==max(Total_tau))
    # tspobj$m.index <- tspobj$m.index[1:seq(3, dim(tspobj$m.index)[1]  ,2)[i], ]
    tspobj$m.index <- tspobj$m.index[1:k.try[i], ]
  } else {

    tspobj$m.index <- tspobj$m.index[1:K.max,]

  }

  rownames(tspobj$m.index) <- rownames(gene.pair.table) <- NULL
  meta.tspobj <- list()
  meta.tspobj$model <- tspobj$m.index
  meta.tspobj$num.K <- dim(tspobj$m.index)[1]
  meta.tspobj$gene.pair.table <- gene.pair.table
  if(is.VO) meta.tspobj$VO <- Total_tau
  return(meta.tspobj)
}

##' Function for predicting subjects class based on given MetaKTSP model, and summarize the results
##'
##'
##' @title Prediction based on a given MetaKTSP model
##' @param test.dat Test dataset.
##' @param tspobj MetaKTSP model.
##' @param K The maximum number of top score pairs.
##' @param is.youden Used youden index in result. Default TRUE.
##' @return An object of class list is returned with following three components:
##' @return "mul.rule" is the prediction result
##' @return "accuracy" is the prediction accuracy
##' @return "youden" is Youden index
##' @author Chien-Wei Lin
##' @export
meta.mul.rule.predict <- function(test.dat, tspobj, K , is.youden=TRUE){
  test.grp <- as.numeric(colnames(test.dat))

  if(K != 1){
    mul.rule <- foreach(i = 1: length(test.grp) ) %dopar% {
      cum.rule <- predict.mtsp(tspobj, test.dat, K=K)[ ,i]
      if(sum(cum.rule)/K == cum.rule[1]){
        cum.rule[1]
      }else{
        ifelse(table(cum.rule)[2] > (K/2) ,  1 ,  0)
      }
    }
    mul.rule <- do.call(cbind, mul.rule)
  }else{
    mul.rule <- predict.mtsp(tspobj, test.dat, K=K)
  }

  youden <- c()
  result<-list()
  if(is.youden){
    .tab <- table( mul.rule, test.grp)
    if(length(rownames(.tab)) == 1) {
      if(rownames(.tab)=="1"){ .tab <- rbind( c(0,0),.tab)
      } else {
        .tab <- rbind(.tab,c(0,0))
      }
    }
    if((length(rownames(.tab)) == 2 & length(colnames(.tab)) == 2)) {youden <- (.tab[1,1]  / sum(.tab[,1]))  + (.tab[2,2] / sum(.tab[,2])) - 1}
    result$youden <- youden
  }
  accuracy <- mean(as.numeric(mul.rule) == test.grp)
  result$mul.rule <- as.numeric(mul.rule)
  result$accuracy <- accuracy
  return(result)
}

##' Function for predicting subjects class based on given MetaKTSP model. For internal use in function meta.mul.rule.predict.
##'
##'
##' @title Prediction based on a given MetaKTSP model
##' @param object MetaKTSP model.
##' @param test.dat A matrix with features in row and samples in column (p x n).
##' @param K The maximum number of top score pairs.
##' @return A matrix of prediction result (K x n) with samples in column and different row correspond to different K.
##' @author Chien-Wei Lin <masaki396@gmail.com>
predict.mtsp <- function(object, test.dat, K){
  tspobj <- object
  test.dat <- t(test.dat)
  predict <- foreach(i = 1 : K, .combine=rbind)%dopar% {
    if( as.numeric(tspobj$model[i, 5]) > 0 ){
      as.numeric(test.dat[ ,colnames(test.dat) == tspobj$model[i, 3]]  >  test.dat[ ,colnames(test.dat) == tspobj$model[i, 4]])
    } else {
      as.numeric(test.dat[ ,colnames(test.dat) == tspobj$model[i, 3]]  <=  test.dat[ ,colnames(test.dat) == tspobj$model[i, 4]])
    }
  }
  return(predict)
}


