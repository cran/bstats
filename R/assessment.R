
.get.counts <- function(x)
    {
        x <- tolower(as.character(x))
        n <- nchar(x)
        tmp <- NULL
        for(i in 1:n){
            tmp <- c(tmp, substr(x,i,i))
        }
        table(tmp)
    }

.get.score <- function(x,ans="A"){
    ans <- tolower(as.character(ans))
    z <- match(ans, names(x))
    out <- 0
    if(!is.na(z)){
        out <- x[z]/sum(x)
    }
    return(out);
}


assessment <- function(x, keys, cutoff=0.60){
    if(missing(keys))
        stop("Solution keys are missing")
    keys <- tolower(as.character(keys))
    nkeys <- length(keys)
    if(ncol(x) != nkeys)
        stop("Answers and solution keys have different lengths")
    scores <- NULL
    counts <- NULL
    for(i in 1:nkeys){
        tmp <- NULL
        xa <- 0; xb <- 0; xc <- 0; xd <- 0;
        for(j in 1:nrow(x)){
            tmp1 <- .get.counts(x[j,i])
            tmp2 <- .get.score(tmp1, keys[i])            
            tmp <- c(tmp, tmp2)
            tmp2 <- .get.score(tmp1, "a"); xa <- xa + tmp2;
            tmp2 <- .get.score(tmp1, "b"); xb <- xb + tmp2;
            tmp2 <- .get.score(tmp1, "c"); xc <- xc + tmp2;
            tmp2 <- .get.score(tmp1, "d"); xd <- xd + tmp2;            
        }
        scores <- cbind(scores, tmp)
        counts <- cbind(counts, c(xa, xb, xc, xd))
    }
    scores <- as.data.frame(scores)
    ##    names(scores) <- names(x)
    tmp1 <- apply(scores, 2, sum) # counts for different Qs
    tmp2 <- apply(scores, 2, mean)*100 # percentage for Qs
    qwr <- rbind(counts, tmp1, tmp2)
    qwr <- round(qwr,2)
    qwr <- as.data.frame(qwr)
    row.names(qwr) <- c("A", "B", "C", "D", "Correct","Rate(%)")
    names(qwr) <- names(x)
    tmp1 <- mean(apply(scores,1,mean) > cutoff)
    list(summary = qwr, pass.rate=tmp1, correct.rate=mean(tmp2))
}
