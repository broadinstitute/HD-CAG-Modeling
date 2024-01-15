
simulateOneStep <- function(cag, expMatrix) {
    #cat(sprintf("Sim %d:\n", cag))
    pdf = expMatrix[cag,]
    #print(pdf)
    cdf = cumsum(pdf)
    #print(cdf)
    p = runif(1)
    idx = max(which(cdf <= p))
    #cat(sprintf("Sim %d: p=%s -> %d\n", cag, p, idx+1))
    return(idx+1)
}

simulateCells <- function(model, ncells, nyears, initialCAG, maxCAG) {
    cat(sprintf("Initial CAG:\n"))
    print(initialCAG)
    result = matrix(0, ncol=ncells, nrow=nyears+1)
    result[1,] = rep(initialCAG, ncells)
    qMatrix = model$generate_Qmatrix(model$fitted_parameters, maxCAG)
    expMatrix = expm::expm(qMatrix)
    for (t in seq(1,nyears)) {
        v = sapply(result[t,], simulateOneStep, expMatrix)
        #cat(sprintf("%s Simulate t=%d:\n", date(), t))
        #print(v)
        result[t+1,] = v
    }
    return(result)
}

