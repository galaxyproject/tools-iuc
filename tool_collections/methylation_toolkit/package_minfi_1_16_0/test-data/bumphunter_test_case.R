
library(bumphunter)


###################################################
### code chunk number 3: clustermakerdata
###################################################
pos <- list(pos1=seq(1,1000,35),
            pos2=seq(2001,3000,35),
            pos3=seq(1,1000,50))
chr <- rep(paste0("chr",c(1,1,2)), times = sapply(pos,length))
pos <- unlist(pos, use.names=FALSE)


###################################################
### code chunk number 4: clustermaker
###################################################
cl <- clusterMaker(chr, pos, maxGap = 300)
table(cl)


###################################################
### code chunk number 5: clusterplot
###################################################
ind <- which(chr=="chr1")


###################################################
### code chunk number 6: simulatedbumps
###################################################
Indexes <- split(seq_along(cl), cl)
beta1 <- rep(0, length(pos))
for(i in seq(along=Indexes)){
    ind <- Indexes[[i]]
    x <- pos[ind]
    z <- scale(x, median(x), max(x)/12)
    beta1[ind] <- i*(-1)^(i+1)*pmax(1-abs(z)^3,0)^3 ##multiply by i to vary size
}


###################################################
### code chunk number 7: getSegments
###################################################
segs <- getSegments(beta1, cl, cutoff=0.05)



###################################################
### code chunk number 9: regionFinder
###################################################
tab <- regionFinder(beta1, chr, pos, cl, cutoff=0.05)


###################################################
### code chunk number 10: simulationOfReps
###################################################
beta0 <- 3*sin(2*pi*pos/720)
X <- cbind(rep(1,20), rep(c(0,1), each=10))
error <- matrix(rnorm(20*length(beta1), 0, 1), ncol=20)
y <- t(X[,1])%x%beta0 + t(X[,2])%x%beta1 + error


###################################################
### code chunk number 11: bumphunter
###################################################
tab <- bumphunter(y, X, chr, pos, cl, cutoff=.5)


###################################################
### code chunk number 12: load-foreach
###################################################
library(doParallel)
registerDoParallel(cores = 2)


###################################################
### code chunk number 13: parallel-bumphunter
###################################################
tab <- bumphunter(y, X, chr, pos, cl, cutoff=.5, B=250, verbose = TRUE)
bumps = tab$tab

save(bumps,file = "bumps.RData",compress=TRUE)
write.csv(bumps,file = "bumps.csv")

###################################################
### code chunk number 14: closeConnetions
###################################################
bumphunter:::foreachCleanup()



