binom.prot.pval <- function(A,B,odds=1) {
   pval <- rep(NA,length(A))
   for (i in 1:length(A)) {
    r <- A[i]
    n <- A[i] + B[i]
    p.hyp <- 1 - odds/(odds + 1)
    if (r < 0) stop(paste("numerator less than 0 for case",i))
    if (r > n) stop(paste("numerator exceeds denominator for case",i))
    if ((round(r) != r) | (round(n) != n)) stop(paste("non-integer for case",i))
    if (r <= n/2) pval[i] <- pbinom(r,n,p.hyp,lower.tail=T)
    if (r > n/2) pval[i] <- pbinom(n-r,n,p.hyp,lower.tail=T)
    }
   pval
   }

fisher.prot <- function(A,B,or=1) {
  A.total <- sum(A)
  B.total <- sum(B)
  fish.mat <- matrix(NA,nrow=length(A),ncol=4)
  for (i in 1:length(A)) {
   # i=16
   mat <- matrix(c(A[i],A.total - A[i],B[i],B.total - B[i]),ncol=2)
   log.or.test <- log(mat[1,1]) + log(mat[2,2]) - log(mat[2,1]) - log(mat[1,2])
   or.test <- exp(log.or.test)
   #or.test <- mat[1,1]*mat[2,2]/(mat[2,1]*mat[1,2])
   if (or.test >= 1) {
      or.alt <- "greater"
      or.hyp <- or
      }
   if (or.test < 1) {
     or.alt <- "less"
     or.hyp <- 1/or
     }

   fisher.out.AB <- fisher.test(mat,or=or.hyp,alternative=or.alt)
   fisher.out.AB.twosided <- fisher.test(mat)
   psi.AB <- as.numeric(fisher.out.AB.twosided$estimate)
   psi.ci.AB <- as.numeric(fisher.out.AB.twosided$conf.int)
   psi.pval.AB <- as.numeric(fisher.out.AB$p.value)

   fish.mat[i,] <- c(psi.AB,psi.ci.AB,psi.pval.AB)
  #result[is.infinite(result)] <- 200   # max is 200
   }
  result <- data.frame(fish.mat)
  names(result) <- c("or","low.or","high.or","pval.or")
  result
  }




proteinCountCompare <- function(A, B) {
  wilson.AB <- binconf(A,A+B,method="wilson")
  wilson.AB <- data.frame(wilson.AB)
  names(wilson.AB) <- c("prop","Wilson.low","Wilson.high")

  n.odds <- length(odds.list)

  pval.AB.mat <- matrix(NA,nrow=length(A),ncol=n.odds)
  qval.AB.mat <- matrix(NA,nrow=length(A),ncol=n.odds)
  for (i in 1:n.odds) {
   pval.AB.i <- binom.prot.pval(A,B,odds=odds.list[i])
   pval.AB.mat[,i] <- pval.AB.i
   qval.AB.i <- qvalue(pval.AB.i,lambda=0)$qvalues
   qval.AB.mat[,i] <- qval.AB.i
   }
  pval.AB <- data.frame(pval.AB.mat)
  qval.AB <- data.frame(qval.AB.mat)

  for (i in 1:n.odds) {
    names.i <- paste("pval.",odds.list[i],sep="")
    names(pval.AB)[i] <- names.i
    names.q.i <- paste("qval.",odds.list[i],sep="")
    names(qval.AB)[i] <- names.q.i
    }


  pval.AB.f.mat <- matrix(NA,nrow=length(A),ncol=n.odds)
  qval.AB.f.mat <- matrix(NA,nrow=length(A),ncol=n.odds)
  for (i in 1:n.odds) {
    # i=16
    fisher.out <- fisher.prot(A,B,or=odds.list[i])
    pval.AB.f.i <- fisher.out[,4]
    qval.AB.f.i <- qvalue(pval.AB.f.i,lambda=0)$qvalues
    pval.AB.f.mat[,i] <- pval.AB.f.i
    qval.AB.f.mat[,i] <- qval.AB.f.i
    }
  pval.AB.f <- data.frame(pval.AB.f.mat)
  qval.AB.f <- data.frame(qval.AB.f.mat)

  for (i in 1:n.odds) {
    names.f.i <- paste("pval.fish.",odds.list[i],sep="")
    names(pval.AB.f)[i] <- names.f.i
    names.f.q.i <- paste("qval.fish.",odds.list[i],sep="")
    names(qval.AB.f)[i] <- names.f.q.i
    }
  fisher.odds <- fisher.out[,1:3]
  result <- list(wilson.AB=wilson.AB, pval.AB=pval.AB, qval.AB=qval.AB,
                 fisher.odds=fisher.odds, pval.AB.f=pval.AB.f, qval.AB.f=qval.AB.f)
  result
  }
