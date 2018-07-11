calculate.false.positive.rate <- function(clust.ref, clust.test)
{
    nmb.false.positive <- length(clust.test) - length(which(clust.test %in% clust.ref))
    nmb.elements.test <- length(clust.test)

    return(nmb.false.negative / nmb.elements.test)

}


calculate.true.positive.rate <- function(clust.ref, clust.test)
{
    nmb.positive <- length(clust.ref)
    nmb.true.positive <- length(which(clust.test %in% clust.ref))

    return(nmb.true.positive / nmb.positive)
}


calculate.precision.recall.for.clust <- function(clust.ref, clust.test)
{
    nmb.real.elements <- length(clust.ref)
    nmb.well.attributed.elements <- length(which(clust.test %in% clust.ref))
    nmb.attributed.elements <- length(clust.test)
    nmb.non.attributed.elements.ref <- length(clust.ref)  -  nmb.well.attributed.elements
    nmb.non.attributed.elements.test <- length(clust.test)  -  nmb.well.attributed.elements

    precision.coef <- nmb.well.attributed.elements / nmb.attributed.elements
    recall.coef <- nmb.well.attributed.elements / nmb.real.elements

    return(list(precision.coef,recall.coef,nmb.non.attributed.elements.ref, nmb.non.attributed.elements.test))

}


calculate.F.beta.coef <- function(precision.coef, recall.coef, beta.coef = 1, threshold = 0)
{
    F.beta.coef <- sapply(1:length(precision.coef), function(i)
    {
        p <- as.numeric(precision.coef[i])
        r <- as.numeric(recall.coef[i])
        Fs <- trunc((p*r)/(p+r)*2*10^4)/10^4
        return(Fs)
    })
    
    F.beta.coef <- F.beta.coef[F.beta.coef>=threshold]
    if (length(F.beta.coef) == 0)
    {
        F.beta.coef <- 0
    }

    return(F.beta.coef)
}


calculate.G.coef <- function(precision.coef, recall.coef, threshold = 0)
{
    G.coef <- sqrt(precision.coef * recall.coef)
    G.coef <- sapply(1:length(precision.coef), function(i)
    {
        return(trunc(sqrt(precision.coef[i] * recall.coef[i])*10^4)/10^4)
    })
    G.coef <- G.coef[G.coef>=threshold]
    if (length(G.coef) == 0)
    {
        G.coef <- 0
    }

    return(G.coef)
}


calculate.FG.fiability.coef <- function(F.coef, G.coef, alpha.coef = 1)
{
    f <- mean(F.coef[!is.na(F.coef)])
    g <- mean(G.coef[!is.na(G.coef)])
    FG.coef <- (alpha.coef*f+g) / (1+alpha.coef)

    return(trunc(FG.coef*10^4)/10^4)
}


calculate.parameter.influence.factor <- function(param, FG.pos.coefs, FG.neg.coefs)
{
    If.pos <- sum(FG.pos.coefs)/length(FG.pos.coefs)
    If.neg <- sum(FG.neg.coefs)/length(FG.neg.coefs)

    return(If.pos - If.neg)
}
