calculate.efficiency.coef <- function(FG.coef, Sf.coef, FG.weight = 1.5, Sf.weight = 1)
{
    E.coef <- ( FG.coef*FG.weight + Sf.coef*Sf.weight ) / ( FG.weight+Sf.weight)
    
    return(E.coef)
}