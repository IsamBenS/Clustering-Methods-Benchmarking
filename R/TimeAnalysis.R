calculate.speed.factor <- function(exec.times, exec.time.test, mean=0.5)
{
    Sf <- NULL
    if(length(exec.times)>1)
    {
        x.m <- 1/max(unlist(exec.times))
        x.M <- 1/min(unlist(exec.times))
        g <- sum(sapply(exec.times,function(e){return(1/e)})) / length(exec.time.test)
        
        A <- t(matrix(c(x.m^2,x.m,1,x.M^2,x.M,1,g^2,g,1), nrow = 3))
        B <- c(0,1,mean)
        X <- solve(A,B, tol = 1e-30)
        
        Sf <- X[1]*(1/exec.time.test)^2 + X[2]*(1/exec.time.test) + X[3]
    }
    
    return(Sf)
}