library(flowCore)
library(Biobase)

add.keyword.to.fcs <- function(fcs, added.keyword, added.keyword.name)
{
    fcs.out <- fcs
    fcs.out@description <- c(fcs.out@description,added.keyword)
    names(fcs.out@description) <- c(names(fcs@description),added.keyword.name)
    
    return(fcs.out)
}

write.FCS.CIPHE <- function(fcs, fcs.path)
{
    metaData <- data.frame(name = dimnames(fcs@exprs)[[2]], desc = dimnames(fcs@exprs)[[2]])
    nmb.dim <- ncol(fcs@exprs)
    metaData$range <- lapply(1:(nmb.dim), function(i)
    {
        return(diff(range(fcs@exprs[,i])))
    })
    metaData$minRange <- sapply(1:(nmb.dim), function(i)
    {
        return(min(fcs@exprs[,i]))
    })
    metaData$maxRange <- sapply(1:(nmb.dim), function(i)
    {
        return(max(fcs@exprs[,i]))
    })
    rownames(metaData) <- sapply(1:(nmb.dim), function(i)
    {
        return(paste("$P",i,sep=""))
    })
    
    fcs.out <- fcs
    parameters(fcs.out) <- AnnotatedDataFrame(metaData)
    
    lapply(c(1:ncol(fcs.out@exprs)),function(x)
    {
        fcs.out@description[[paste0("$P",x,"R")]] <<- metaData$maxRange[[x]]
        fcs.out@description[[paste0("$P",x,"B")]] <<- "32"
        fcs.out@description[[paste0("$P",x,"S")]] <<- colnames(fcs.out@exprs)[x]
    })
    
    write.FCS(fcs.out, fcs.path, delimiter = '#')
}


keyword.exists.FCS <- function(fcs, key.part)
{
    key.exists <- F
    i <- 1
    while(!key.exists && i <= length(fcs@description))
    {
        if(grepl(key.part, names(fcs@description)[i], fixed = T))
        {
           key.exists <- T
        }
        i <- i+1
    }
    
    return(key.exists)
}

get.keywords.with.keypart.FCS <- function(fcs, key.part)
{
    keywords.list <- c()
    i <- 1
    while(i <= length(fcs@description))
    {
        if(grepl(key.part, names(fcs@description)[i], fixed = T))
        {
            keywords.list <- c(keywords.list,fcs@description[[i]])
            names(keywords.list)[length(keywords.list)] <- names(fcs@description)[i]
        }
        i <- i+1
    }
    
    return(keywords.list)
}