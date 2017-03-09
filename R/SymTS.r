dCTS <- function(x, alpha, c=1, ell=1, mu=0){
    if(c<= 0){
        stop("c must be greater than 0")
    }
    if(ell<= 0){
        stop("ell must be greater than 0")
    }
    if(alpha < 0 | alpha >= 2){
        stop("alpha must be greater than or equal to 0 and less than or equal to 2")
    }
    if(alpha == 0){
       result = .C(C_dCTS0, as.double(x), as.integer(length(x)), as.double(mu), as.double(alpha), as.double(c), as.double(ell), as.double(vector("double", length(x))))[[7]]
        if(c <= 0.5){
            result[abs(x-mu)<.00001] = Inf
        }

    }
    else
        .C(C_dCTS, as.double(x), as.integer(length(x)), as.double(mu), as.double(alpha), as.double(c), as.double(ell), as.double(vector("double", length(x))))[[7]]
}

pCTS <- function(x, alpha, c=1, ell=1, mu=0){
    if(c<= 0){
        stop("c must be greater than 0")
    }
    if(ell<= 0){
        stop("ell must be greater than 0")
    }
    if(alpha < 0 | alpha >= 2){
        stop("alpha must be greater than or equal to 0 and less than or equal to 2")
    }
    .C(C_pCTS, as.double(x), as.integer(length(x)), as.double(mu), as.double(alpha), as.double(c), as.double(ell), as.double(vector("double", length(x))))[[7]]
}

qCTS <- function(x, alpha, c=1, ell=1, mu=0){
    if(sum(x<=0)>0 | sum(x>=1)>0){
        stop("x must be between 0 and 1")
    }
    if(c<= 0){
        stop("c must be greater than 0")
    }
    if(ell<= 0){
        stop("ell must be greater than 0")
    }
    if(alpha < 0 | alpha >= 2){
        stop("alpha must be greater than or equal to 0 and less than or equal to 2")
    }
    .C(C_qCTS, as.double(x), as.integer(length(x)), as.double(mu), as.double(alpha), as.double(c), as.double(ell), as.double(vector("double", length(x))))[[7]]
}
dPowTS <- function(x, alpha, c=1, ell=1, mu=0){
    if(c<= 0){
        stop("c must be greater than 0")
    }
    if(ell<= 0){
        stop("ell must be greater than 0")
    }
    if(alpha < 0 | alpha >= 2){
        stop("alpha must be greater than or equal to 0 and less than or equal to 2")
    }
    if(alpha <= .01 & c <= 0.5*(1+ell)){
        stop("when alpha is close to 0, c/(1+ell) must be greater than 0.5")
    }
    .C(C_dPowTS, as.double(x), as.integer(length(x)), as.double(mu), as.double(alpha), as.double(c), as.double(ell), as.double(vector("double", length(x))))[[7]]
}
pPowTS <- function(x, alpha, c=1, ell=1, mu=0){
    if(c<= 0){
        stop("c must be greater than 0")
    }
    if(ell<= 0){
        stop("ell must be greater than 0")
    }
    if(alpha < 0 | alpha >= 2){
        stop("alpha must be greater than or equal to 0 and less than or equal to 2")
    }
    .C(C_pPowTS, as.double(x), as.integer(length(x)), as.double(mu), as.double(alpha), as.double(c), as.double(ell), as.double(vector("double", length(x))))[[7]]
}

qPowTS <- function(x, alpha, c=1, ell=1, mu=0){
    if(sum(x<=0)>0 | sum(x>=1)>0){
        stop("x must be between 0 and 1")
    }
    if(c<= 0){
        stop("c must be greater than 0")
    }
    if(ell<= 0){
        stop("ell must be greater than 0")
    }
    if(alpha < 0 | alpha >= 2){
        stop("alpha must be greater than or equal to 0 and less than or equal to 2")
    }
    .C(C_qPowTS, as.double(x), as.integer(length(x)), as.double(mu), as.double(alpha), as.double(c), as.double(ell), as.double(vector("double", length(x))))[[7]]
}

dSaS <- function(x, alpha, c=1, mu=0){
    if(c<= 0){
        stop("c must be greater than 0")
    }
    if(alpha <= 0 | alpha >= 2){
        stop("alpha must be greater than 0 and less than 2")
    }
    if(alpha == 1)
        dcauchy(x, location = mu, scale = c)
    else
        .C(C_dSaS, as.double(x), as.integer(length(x)), as.double(mu), as.double(alpha), as.double(c), as.double(vector("double", length(x))))[[6]]
}

pSaS <- function(x, alpha, c=1, mu=0){
    if(c<= 0){
        stop("c must be greater than 0")
    }
    if(alpha <= 0 | alpha >= 2){
        stop("alpha must be greater than 0 and less than 2")
    }
    if(alpha == 1)
        pcauchy(x, location = mu, scale = c)
    else
        .C(C_pSaS, as.double(x), as.integer(length(x)), as.double(mu), as.double(alpha), as.double(c), as.double(vector("double", length(x))))[[6]]
}

qSaS <- function(x, alpha, c=1, mu=0){
    if(sum(x<=0)>0 | sum(x>=1)>0){
        stop("x must be between 0 and 1")
    }
    if(c<= 0){
        stop("c must be greater than 0")
    }
    if(alpha <= 0 | alpha >= 2){
        stop("alpha must be greater than 0 and less than 2")
    }
    if(alpha == 1)
       qcauchy(x, location = mu, scale = c)
    else
        .C(C_qSaS, as.double(x), as.integer(length(x)), as.double(mu), as.double(alpha), as.double(c), as.double(vector("double", length(x))))[[6]]
}

rPowTS <- function(r, alpha, c=1, ell=1, mu=0){
    qPowTS(runif(r), alpha, c, ell, mu)
}

rCTS <- function(r, alpha, c=1, ell=1, mu=0){
    qCTS(runif(r), alpha, c, ell, mu)
}

rSaS <- function(r, alpha, c=1, mu=0){
    g <- runif(r)*pi-pi/2
    W <- rexp(r)
    c^(1/alpha) * sin(alpha*g) * (cos(g))^(-1/alpha) *(cos((1-alpha)*g)/W)^(-1+1/alpha) + mu
}
