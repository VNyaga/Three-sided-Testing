
#========================================================================================================
#   Testing
#========================================================================================================

z = function(x){
    n11 = x[1]
    n10 = x[2] 
    n01 = x[3] 
    n00 = x[4] 
    delta = x[5]
    
    n = sum(n11, n01, n10, n00)
    
    p1_hat = (n11 + n10)/n
    p0_hat = (n11 + n01)/n
    p01_hat = n01/n
    p10_hat = n10/n
    p00_hat = n00/n
    
    p10_tilde = (-p1_hat + delta^2 *(p0_hat + 2*p10_hat) + sqrt((p1_hat - delta^2*p0_hat)^2 + 4* delta^2 * p10_hat * p01_hat) )/(2*delta*(delta + 1))
    p01_tilde = delta * p10_tilde - (delta - 1)*(1 - p00_hat)
    
    z = (sqrt(n)*(p1_hat - delta*p0_hat))/sqrt(delta*(p10_tilde + p01_tilde));
    z
}

#========================================================================================================
#   Three sided Confidence Intervals
#========================================================================================================
CI_3_sided = function(df, alpha=0.05, dp=3) {
    delta_l=df[5] 
    delta_u=df[6]

    cml_ci = function(df, alpha=alpha) {
        n11 = df[1]
        n10 = df[2]
        n01 = df[3]
        n00 = df[4]
        
        n = sum(n11, n01, n10, n00)

        p1_hat = (n11 + n10)/n
        p0_hat = (n11 + n01)/n
        p01_hat = n01/n
        p10_hat = n10/n
        p00_hat = n00/n

        RR = p1_hat/p0_hat
        #Sample based CI
        #Undefined for off-diagnols = 0
        # 0 for pt = 0
        inits = exp(log(RR) + c(-1, 1)*qnorm(1 - alpha/2)*sqrt((n01 + n10)/((n11 + n10)*(n11 + n01))));

        # CML based CI

        fun = function(RR){
            p10_tilde =(-p1_hat + (RR^2)*(p0_hat + 2*p10_hat) + sqrt((p1_hat - (RR^2)*p0_hat)^2 + 4*(RR^2)*p10_hat*p01_hat))/(2*RR*(RR + 1))
            p01_tilde = RR*p10_tilde - (RR - 1)*(1 - p00_hat)
            z_alpha = qnorm(alpha/2)^2
            n*((p1_hat - RR*p0_hat)^2)/(RR*(p10_tilde + p01_tilde)) - z_alpha
        }

        ci = uniroot.all(fun, c(inits[1]*0.5, inits[2]*1.5))

        return(ci)
    }
    
    ci_alpha = cml_ci(df=c(df[1:4]), alpha*2)
    ci_alpha_2 = cml_ci(df=c(df[1:4]), alpha)
    
    #Lower
    if (ci_alpha[1] < delta_l) {
        
        newl = paste('(', digits(ci_alpha[1], digits = dp, format = "f"), sep = '')
        
    } else if ((ci_alpha_2[1] <= delta_l) & (delta_l <= ci_alpha[1])){
        
        newl = paste('[', digits(delta_l, digits = dp, format = "f"), sep = '')
        
    } else if((delta_l < ci_alpha_2[1]) & (ci_alpha_2[1] < delta_u)){
        
        newl = paste('(', digits(ci_alpha_2[1], digits = dp, format = "f"), sep = '')
        
    } else if (delta_u <= ci_alpha_2[1]){
        
        newl= paste('(', digits(delta_u, digits = dp, format = "f"), sep = '')
        
    }
    
    #Upper
    if (ci_alpha_2[2] <= delta_l){
        
        newu = paste(digits(delta_l, digits = dp, format = "f"), ')', sep='')
        
    } else if((delta_l < ci_alpha_2[2]) & (ci_alpha_2[2] < delta_u)){
        
        newu= paste(digits(ci_alpha_2[2], digits = dp, format = "f"), ')', sep='')
        
    } else if((ci_alpha[2] <= delta_u) & (delta_u <= ci_alpha_2[2])){
        
        newu= paste(digits(delta_u, digits = dp, format = "f"), ']', sep='')
        
    } else if(delta_u < ci_alpha[2]){
        
        newu = paste(digits(ci_alpha[2], digits = dp, format = "f"), ')', sep='')
    }
    
    new.ci = paste(newl, newu, sep = ', ')
    return(new.ci)
}

#============================================================================
#      P-values
#============================================================================

pval = function(x, dp=3) {
    delta_l = x[5]
    delta_u = x[6]
    
    noninf_z = z(x=x[1:5])
    super_z = z(x=x[c(1:4, 6)])
    
    noninf_p = digits(1 - pnorm(noninf_z), digits = dp, format = "f")
    super_p = digits(1-pnorm(super_z), digits = dp, format = "f")
    equiv_p2 = digits(pnorm(super_z), digits = dp, format = "f")
    equiv_p = digits(max(noninf_p, equiv_p2), digits = dp, format = "f")
    return(list(noninf_p=noninf_p, equiv_p= equiv_p, super_p=super_p))
}


#============================================================================
#      Classical confidence intervals
#============================================================================

cml_ci = function(df, alpha=0.05) {
    n11 = df[1] 
    n10 = df[2]
    n01 = df[3]
    n00 = df[4]

    n = sum(n11, n01, n10, n00)
    
    p1_hat = (n11 + n10)/n
    p0_hat = (n11 + n01)/n
    p01_hat = n01/n
    p10_hat = n10/n
    p00_hat = n00/n
    
    RR = p1_hat/p0_hat 
    #Sample based CI
    #Undefined for off-diagnols = 0
    # 0 for pt = 0
    inits = exp(log(RR) + c(-1, 1)*qnorm(1 - alpha/2)*sqrt((n01 + n10)/((n11 + n10)*(n11 + n01))));
    
    # CML based CI
    
    fun = function(RR){
        p10_tilde =(-p1_hat + (RR^2)*(p0_hat + 2*p10_hat) + sqrt((p1_hat - (RR^2)*p0_hat)^2 + 4*(RR^2)*p10_hat*p01_hat))/(2*RR*(RR + 1))
        p01_tilde = RR*p10_tilde - (RR - 1)*(1 - p00_hat)
        z_alpha = qnorm(alpha/2)^2
        n*((p1_hat - RR*p0_hat)^2)/(RR*(p10_tilde + p01_tilde)) - z_alpha
    }
    
    ci = uniroot.all(fun, c(inits[1]*0.5, inits[2]*1.5))
    
    return(ci)
}


#============================================================================
#      Nice printout
#============================================================================

Three_sided_test = function(a, b, c, d, delta_l, delta_u, alpha=0.05, dp=3) {
    x = matrix(c(a, b, c, d), nrow=2, byrow = TRUE)
    rownames(x) = c('New Test = 1', 'New Test = 0')
    colnames(x) = c('Comparator Test = 1', 'Comparator Test = 0')
    
    z_l = z(x=c(a, b, c, d, delta_l))
    z_u = z(x=c(a, b, c, d, delta_u))
    
    ci = CI_3_sided(df=c(a, b, c, d, delta_l, delta_u), alpha=alpha, dp=dp)
    ci_split = unlist(strsplit(ci, ','))
    ci_lower = substr(ci_split[1], 2, 20)
    ci_upper = substr(ci_split[2], 2, nchar(ci_split[2])-1)
    
    p = pval(x=c(a, b, c, d, delta_l, delta_u), dp=dp)
    
    cat("Three-sided testing procedure\n\n")
    cat("Data:")
    print(x)
    cat("\n")
    cat(paste("H0: RR < ", digits(delta_l, digits = 3, format = "f"), " or RR > ", digits(delta_u, digits = 3, format = "f"),
              " versus H1: ", digits(delta_l, digits = 3, format = "f"), " <= RR <= ", digits(delta_u, digits = 3, format = "f"), "\n", sep=''))
    cat("Test statistic (1) = "); message(digits(z_l, digits = 3, format = "f"));
    cat("Test statistic (2) = "); message(digits(z_u, digits = 3, format = "f")); 
    cat('P-value = '); message(p$equiv_p)
    
    cat("\n")
    cat(paste("H0: RR < ", digits(delta_l, digits = 3, format = "f"), " versus H1:  RR >= ", digits(delta_l, digits = 3, format = "f"), "\n"))
    cat("Test statistic = "); message(digits(z_l, digits = 3, format = "f")); cat('P-value = '); message(p$noninf_p)
    
    cat("\n")
    cat(paste("H0: RR <= ", digits(delta_u, digits = 3, format = "f"), " versus H1:  RR > ", digits(delta_u, digits = 3, format = "f"), "\n"))
    cat("Test statistic = "); message(digits(z_u, digits = 3, format = "f")); cat('P-value = '); message(p$super_p)
    
    cat("\n")
    cat(paste((1-alpha)*100, "% CI: ")); message(ci)
    
    cat("\n")
    cat("Conclusion: ")
    if(ci_upper < delta_l) message('New test is inferior') 
    if(ci_lower > delta_u) message('New test is superior')
    if(ci_lower >= delta_l & ci_upper <= delta_u) message('New test is similar')
    if(ci_upper < delta_u & ci_lower < delta_l) message('New test is non-superior') 
    if(ci_upper > delta_u & ci_lower > delta_l) message('New test is non-inferior')
    
}
