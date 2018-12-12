## File Name: starts_uni_estimate_variance_proportions_pml.R
## File Version: 0.17

starts_uni_estimate_variance_proportions_pml <- function( coef, vcov, vars )
{
    NV <- length(vars)
    parm_fct <- function(x){
        tot <- sum(x)
        res <- x / tot
        return(res)
    }
    # define point for evaluating partial derivatives
    par <- coef[vars]
    val <- parm_fct(x=par)

    #--- compute gradient
    A <- matrix( 0, nrow=NV, ncol=NV)
    rownames(A) <- vars
    colnames(A) <- vars
    indices1 <- match( vars, colnames(A) )
    vcov <- vcov[indices1, indices1]
    indices <- 1:NV
    for (cc in indices){
        res <- CDM::numerical_Hessian_partial( par=par[vars], FUN=parm_fct, coordinate=cc )
        # A[,cc] <- res$grad
        A[cc,] <- res$grad
    }
    vcov_prop <- A %*% vcov %*% t(A)
    var_prop <- data.frame( parm=vars, est=val, se=sqrt( diag(vcov_prop) ) )
    var_prop$t <- var_prop$est / var_prop$se
    var_prop$p <- 2*stats::pnorm( - abs( var_prop$t ) )
    return(var_prop)
}
