library(flowCore)
library(Biobase)

add.keyword.to.fcs <- function(fcs, added.keyword, added.keyword.name)
{
    fcs.out <- fcs
    descR <- description(fcs.out)
    descR[[added.keyword.name]] <- added.keyword
    fcs.out <- flowFrame(fcs.out@exprs,description=descR)
    
    return(fcs.out)
} 

write.FCS.CIPHE <- function(fcs, fcs.path)
{
    descR <- description(fcs)
    for(x in 1:ncol(fcs@exprs))
    {
		if(!is.null(descR[[paste0("$P",x,"R")]]))
		{
			descR[[paste0("$P",x,"R")]] <- 262144
		}
		else
		{
			descR <- c(descR, 262144)
			names(descR)[length(descR)] <- paste0("$P",x,"R")
		}
    }
    fcs.out <- flowFrame(fcs@exprs, description = descR)
    fcs.out@description <- descR
    
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