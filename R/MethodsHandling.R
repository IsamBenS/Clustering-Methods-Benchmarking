library(microbenchmark)

benchmark.method <- function(method.name, fcs.file, directory, parameters = list(), markers_cols)
{
    fct.name <- paste0("BRP_BM.",method.name,".execute")
    fct <- match.fun(fct.name, descend = FALSE)
    fcs.output <- NULL
    time.output <- mean(microbenchmark({    fcs.output <- fct(fcs.file, directory, parameters, markers_cols)    }, times = 1, unit = "ns")$time)

    method.output <- list(fcs.output, time.output)

    return(method.output)
}

save.exec.time <- function(method.name, directory, time.output)
{
    times.vector <- NULL
    if("times.csv" %in% list.files(directory))
    {
        times.vector <- c(as.matrix(read.csv(paste0(directory,"times.csv")))[1,])
        if(!(method.name %in% names(times.vector)))
        {
            times.vector[[length(times.vector) + 1]] <- time.output
            names(times.vector)[length(times.vector)] <- method.name
        }
        else
        {
            times.vector[[method.name]] <- time.output
        }

    }
    else
    {
        times.vector <- time.output
        names(times.vector) <- method.name
    }
    mat <- matrix(times.vector, nrow = 1)
    colnames(mat) <- names(times.vector)
    if(file.exists(paste0(directory,"times.csv")))
    {
        file.remove(paste0(directory,"times.csv"))
    }
    write.csv(mat, paste0(directory,"times.csv"), row.names = FALSE)
}

read.exec.time <- function(method.name, directory)
{
    time.output <- NULL
    if("times.csv" %in% list.files(directory))
    {
        temp <- as.vector(read.csv(paste0(directory,"times.csv")))
        if(method.name %in% names(temp))
        {
            time.output <- temp[method.name]
        }
    }
    
    return(time.output)
}
