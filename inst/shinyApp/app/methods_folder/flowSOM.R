library(FlowSOM)
library(Biobase)

fct.parameters <- list("xgrid"=10,"ygrid"=10)

BRP_BM.flowSOM.execute <- function(fcs.file, directory, params = list(10,10), markers_col)
{
    fcs.out.SOM <- ReadInput(fcs.file, transform = FALSE, scale = FALSE)
    fcs.out.SOM <- BuildSOM(fcs.out.SOM, colsToUse = markers_col, 
							xdim=as.numeric(params[[1]]), 
							ydim=as.numeric(params[[2]]))
    fcs.labels <- fcs.out.SOM$map$mapping[,1]
    fcs.out.mat <- cbind(fcs.file@exprs,fcs.labels)
    colnames(fcs.out.mat) <- c(colnames(fcs.file),"cluster_FlowSOM")
	
	metaData <- data.frame(name = dimnames(fcs.out.mat)[[2]], desc = dimnames(fcs.out.mat)[[2]])
    nmb.dim <- ncol(fcs.out.mat)-1
    metaData$range <- lapply(1:(nmb.dim+1), function(i)
    {
        return(diff(range(fcs.out.mat[,i])))
    })
    metaData$minRange <- sapply(1:(nmb.dim+1), function(i)
    {
        return(min(fcs.out.mat[,i]))
    })
    metaData$maxRange <- sapply(1:(nmb.dim+1), function(i)
    {
        return(max(fcs.out.mat[,i]))
    })
    rownames(metaData) <- sapply(1:(nmb.dim+1), function(i)
    {
        return(paste("$P",i,sep=""))
    })
    fcs.out <- flowFrame(exprs=fcs.out.mat, parameters=AnnotatedDataFrame(metaData))

    lapply(c(1:dim(fcs.out.mat)[2]),function(k)
    {
        fcs.out@description[[paste0("$P",k,"R")]] <<- metaData$maxRange[[k]]
        fcs.out@description[[paste0("$P",k,"B")]] <<- "32"
        fcs.out@description[[paste0("$P",k,"S")]] <<- colnames(fcs.out.mat)[k]
    })
    parameters(fcs.out)$maxRange <- sapply(1:ncol(fcs.out.mat), function(i)
    {
        return(max(fcs.out.mat[,i]))
    })
    
    # #PLOT===================================================================================================================================================
    # p <- ncol(fcs.out.mat) - 1
    # n <- as.integer(sqrt(p)) + 1
    # png(output.plot.name)
    # par(mfrow=c(n,n))
    # for (i in 1:p)
    # {
    #     hist(fcs.out.mat[,i], main=paste("PARAM",i))
    # }
    # dev.off()
    # print(colnames(fcs.out))

    return(fcs.out)
}
