library(flowCore)
library(Biobase)

#EVERY FCS FILE MUST HAVE A LAST COLUMN WITH THE CLUSTER ID OF EACH EVENT
FPH.load.folder <- function(path.to.fcs.folder)
{
    folder.content <- list.files(path.to.fcs.folder, pattern=".fcs")
    fcs.files <- lapply(folder.content, function(f.path)
    {
        return(read.FCS(paste0(path.to.fcs.folder,f.path), emptyValue = FALSE))
    })

    return(fcs.files)
}

FPH.get.file.name <- function(fcs.file)
{
    f.name <- strsplit(fcs.file@description$FILENAME, "/")
    f.name <- unlist(f.name[[1]][[length(f.name[[1]])]])
    f.name <- substr(f.name,1,nchar(f.name)-4)

    return(f.name)
}

FPH.create.files.folders <- function(fcs.files, results.path)
{
    lapply(fcs.files, function(f)
    {
        dir.create(paste0(results.path, FPH.get.file.name(f), "/"))
    })
}

FPH.get.file.information <- function(fcs.file, clust.col)
{
    f <- fcs.file
    f.mat <- f@exprs
    clust.dim <- min(ncol(f.mat),clust.col)
    clust.names <- unique(f@exprs[,clust.dim])
    clust.names <- clust.names[clust.names!=0]
    nmb.clust <- length(clust.names)
    
    clusters.data <- lapply(1:nmb.clust, function(c)
    {
        events.id <- which(f.mat[,clust.dim] %in% clust.names[c])
        nmb.events <- length(events.id)
        
        out.list <- NULL
        if(nmb.events > 0)
        {
            out.list <- list(nmb.events, events.id)
        }
        return(out.list)
    })
    clusters.data[ sapply(clusters.data, is.null) ] <- NULL
    return(clusters.data)
}

FPH.get.clusters.associations <- function(clusters.ref, clusters.test)
{
    l <- length(clusters.test)
    mal.theta.local <- 0
    pr.local.coefs <- 0
    clust.rel.ids <- 0
    clust.rel.coefs <- NULL
    nmb.events <- 0
    if(length(clusters.ref) > 0)
    {
        nmb.events <- sum(sapply(c(1:length(clusters.ref)), function(i)
        {
            return( clusters.ref[[i]][[1]] )
        }))
    }
    
    pr.coefs.mat <- matrix(nrow = length(clusters.test), ncol = length(clusters.ref))
    rec.coefs.mat <- matrix(nrow = length(clusters.test), ncol = length(clusters.ref))
    out.el.mat.test <- matrix(nrow = length(clusters.test), ncol = length(clusters.ref))
    if(length(clusters.ref) > 0)
    {
        lapply(1:length(clusters.test), function(r)
        {
            lapply(1:length(clusters.ref), function(c)
            {
                t <- calculate.precision.recall.for.clust(clusters.ref[[c]][[2]], clusters.test[[r]][[2]])
                pr.coefs.mat[r,c] <<- t[[1]]
                rec.coefs.mat[r,c] <<- t[[2]]
                out.el.mat.test[r,c] <<- t[[4]]
            })
        })
    }


    id <- rep(0,l)
    clust.rel.coefs <- list(NULL,NULL)
    clust.rel.id <- 2
    if(l > 0)
    {
        clust.rel.coefs <- c(2,lapply(c(1:l), function(k)
        {
            maxVal <- 0
            t <- 1
            while(t <= length(clusters.ref))
            {
                if (maxVal <= sqrt(pr.coefs.mat[k,t]*rec.coefs.mat[k,t]))
                {
                    maxVal <- sqrt(pr.coefs.mat[k,t]*rec.coefs.mat[k,t])
                    id[k] <<- t
                }
                t <- t+1
            }
            mal.theta.local <<- mal.theta.local + out.el.mat.test[k,id[k]]
            return(list(pr.coefs.mat[k,id[k]],rec.coefs.mat[k,id[k]]))
        }))
        clust.rel.ids <- c(2,id)
    }

    theta.coef <- exp(-7/abs( mal.theta.local / nmb.events ))
    return(list(theta.coef, clust.rel.ids, clust.rel.coefs))
}

FPH.get.clusters.resemblance <- function(file.1, file.2, m1.cols, m2.cols, cl.col.1, cl.col.2, nmb.ev.col.1, nmb.ev.col.2)
{
    clust.res.ids <- 0
    clust.res.coefs <- NULL
    l <- min( nrow(file.1@exprs), nrow(file.2@exprs))
    size.max <- 0
    if(l < nrow(file.2@exprs))
    {
        size.max <- max(unlist(lapply(1:l, function(i){return(file.1@exprs[i,nmb.ev.col.1])})))
    }
    else
    {
        size.max <- max(unlist(lapply(1:l, function(i){return(file.2@exprs[i,nmb.ev.col.2])})))
    }
    
    dist.1.mat <- matrix(2,nrow = nrow(file.1@exprs), ncol = nrow(file.1@exprs))
    dist.2.mat <- matrix(2,nrow = nrow(file.2@exprs), ncol = nrow(file.2@exprs))
    lapply(2:nrow(file.1@exprs), function(c)
    {
        lapply(1:(c-1), function(r)
        {
            dist.1.mat[r,c] <<- sqrt(sum( (file.1@exprs[r,m1.cols] - file.1@exprs[c,m1.cols])^2 ))
        })
    })
    lapply(2:nrow(file.2@exprs), function(c)
    {
        lapply(1:(c-1), function(r)
        {
            dist.2.mat[r,c] <<- sqrt(sum( (file.2@exprs[r,m2.cols] - file.2@exprs[c,m2.cols])^2 ))
        })
    })
    
    smfi.mat <- matrix(nrow = nrow(file.2@exprs), ncol = nrow(file.1@exprs))
    sent.mat <- matrix(nrow = nrow(file.2@exprs), ncol = nrow(file.1@exprs))
    s.mat <- matrix(nrow = nrow(file.2@exprs), ncol = nrow(file.1@exprs))
    lapply(1:nrow(file.1@exprs), function(c)
    {
        c0 <- which(dist.1.mat[c,] == min(dist.1.mat[c,]))[[1]]
        lapply(1:nrow(file.2@exprs), function(r)
        {
            r0 <- which(dist.2.mat[r,] == min(dist.2.mat[r,]))[[1]]
            smfi.mat[r,c] <<- sqrt(sum( (file.1@exprs[r,m1.cols] - file.2@exprs[c,m2.cols])^2 ))
            sent.mat[r,c] <<- sqrt(sum( (file.1@exprs[r0,m1.cols] - file.2@exprs[c0,m2.cols])^2 ))
        })
    })
    d0 <- sqrt(length(m1.cols) * (4.5^2))
    smfi.mat <- exp(-(smfi.mat/0.2/d0))
    sent.mat <- exp(-(sent.mat/0.3/d0))
    s.mat <- (1.4*sent.mat + smfi.mat)/2.4
    
    
    id <- rep(0,l)
    s.ret <- NULL
    s.mfi.ret <- NULL
    s.ent.ret <- NULL
    if(l < nrow(file.2@exprs))
    {
        lapply(1:l, function(i)
        {
            k = which(s.mat[i,] == max(s.mat[i,]))
            if(length(k) > 1)
            {
                k <- which(!(k %in% id))[[1]]
            }
            id[i] <<- k
            s.ret <<- c(s.ret,s.mat[i,k])
            s.mfi.ret <<- c(s.mfi.ret,smfi.mat[i,k])
            s.ent.ret <<- c(s.ent.ret,sent.mat[i,k])
        })
        id <- c(1,id)
    }
    else
    {
        lapply(1:l, function(i)
        {
            k = which(s.mat[,i] == max(s.mat[,i]))
            if(length(k) > 1)
            {
                k <- which(!(k %in% id))[[1]]
            }
            id[i] <<- k
            s.ret <<- c(s.ret,s.mat[k,i])
            s.mfi.ret <<- c(s.mfi.ret,smfi.mat[k,i])
            s.ent.ret <<- c(s.ent.ret,sent.mat[k,i])
        })
        id <- c(2,id)
    }
    
    
    s.correction <- (length(which(!(1:l %in% id[-1])))/l)^6
    s.ret <- sapply(1:l, function(i){return(max(0,s.ret[i] - s.correction))})
    s.mfi.ret  <- sapply(1:l, function(i){return(max(0,s.mfi.ret[i] - s.correction))})
    s.ent.ret <- sapply(1:l, function(i){return(max(0,s.ent.ret[i] - s.correction))})
    
    return(list(id, s.ret, s.mfi.ret, s.ent.ret))
}

FPH.create.artificials.clusters <- function(fcs.test, purity.matrix, cl.col, threshold = 0)
{
    fcs.out <- fcs.test
    clusters.test <- FPH.get.file.information(fcs.test, cl.col)
    if(nrow(purity.matrix)>ncol(purity.matrix))
    {
        fcs.out.mat <- cbind(fcs.test@exprs,rep(0,nrow(fcs.test@exprs)))
        nc <- ncol(fcs.out.mat)
        lapply(1:length(clusters.test), function(i)
        {
            clust.id <- which(purity.matrix[i,] == max(purity.matrix[i,]))[[1]]
            lapply(clusters.test[[i]][[2]], function(ev)
            {
                fcs.out.mat[ev,nc] <<- 0
                if(purity.matrix[i,clust.id] >= threshold)
                {
                    fcs.out.mat[ev,nc] <<- clust.id
                }
            })
        })
        
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
    }
    
    return(fcs.out)
}

FPH.get.purity.matrix <- function(fcs.ref, fcs.test, clust.col.1, clust.col.2)
{
    clusters.ref <- FPH.get.file.information(fcs.ref, clust.col.1)
    clusters.test <- FPH.get.file.information(fcs.test, clust.col.2)
    purity.matrix <- matrix(0,nrow = length(clusters.test), ncol = length(clusters.ref))
    
    lapply(c(1:length(clusters.test)), function(l)
    {
        lapply(c(1:length(clusters.ref)), function(c)
        {
            purity.matrix[l,c] <<- length(which(clusters.test[[l]][[2]] %in% clusters.ref[[c]][[2]])) / clusters.test[[l]][[1]]
        })
    })
    
    return(purity.matrix)
}
    



