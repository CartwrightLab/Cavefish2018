#!/usr/bin/Rscript --vanilla

# fecundity selection

# Simulate a k-loci model
# G: Generations
# N: Number of diploid fish
# k: Number of loci/qtl
# m: immigration rate
# s: selection strength.
#     Individual with 2k blindness alleles has a fitness of 1+s
# h: dominance
#     1 blindness alleles at a locus adds h*s/2k
#     2 blindness alleles at a locus add s/2k
# u: mutation rate from sightedness to blindness
# Q: migration rate
# q0: initial allele frequency

# Initialize population based on the expected allele frequencies
sim0 = function(N,k,q) {
    matrix(rbinom(N*k,2,q),nrow=k)
}

# Calculate the fitness of each individual in the population
# based on their multi-locus genotype
fitness = function(x,k,s,h) {
    z = s/(2*k)
    w = z*(x - ifelse(x==1, (1-2*h),0))
    1+colSums(w)    
}

# Do one generation of simulaiton
simnext = function(x,N,k,m,s,h,u,Q) {
    # calculate the fitness of every individual
    w = fitness(x,k,s,h)
    
    # calculate allele frequencies after selection
    f = sweep(x,2,w,"*")
    q = rowSums(f)/(2*sum(w))

    # migration
    q = q*(1-m) + m*Q

    # mutation
    q = q+(1-q)*u

    # drift, meiosis, fertilization
    sim0(N,k,q)
}

# Run a simulation and keep track of population every 1000 generations
sim = function(G1,G2,N,k,m,s,h,u,Q,q0) {
    ret = list()
    x = sim0(N,k,q0)
    ret[["0"]] = x
    g = 1
    while(g <= G1) {
        x = simnext(x,N,k,0,s,h,u,Q)
        if(g %% 1000 == 0) {
            ret[[as.character(g)]] = x
        }
        g = g + 1
    }
    while(g <= G1+G2) {
        x = simnext(x,N,k,m,s,h,u,Q)
        if(g %% 1000 == 0) {
            ret[[as.character(g)]] = x
        }
        g = g + 1
    }
    ret
}


# simulate and calculate statistics
main = function(p) {
    # simulate
    x = do.call(sim, p)

    # calculate average fitness
    w = sapply(x, fitness, p$k,p$s,p$h)
    w = matrix(w,ncol=length(x),nrow=p$N)
    a = apply(w,2,mean)

    # calculate average phenotype
    f = sapply(x, fitness, p$k,1,p$h)
    f = matrix(f,ncol=length(x),nrow=p$N)
    f = apply(f,2,mean)-1

    # calculate average allele frequency
    q = sapply(x, sum)/(2*p$k*p$N)
    # calculate genetic load
    l = (1+p$s-a)/(1+p$s)

    d = data.frame(g=names(x),q=q,a=a,f=f,l=l)
    d
}

# If run as a script, do work
if(!interactive()) {
    # specify parameter order
    params_list = names(formals(sim))
    params = vector("list", length(params_list))
    names(params) = params_list

    ARGS = commandArgs(trailingOnly = TRUE)
    if(length(ARGS) == 0 || ARGS[1] == "") {
        # if no args passed, print header
        params[params_list] = c(0,0,1,1,0,0,0,0,0,0)
        d = main(params)
        d = cbind(params,d)
        d = names(d)
        cat(paste0(d,collapse=","))
        cat("\n")
        quit(status=0)
    }
    # setup parameters
    for( s in strsplit(ARGS,"=")) {
        params[[s[1]]] = as.numeric(s[2])
    }
    # remove unknown parameters
    params[!(names(params) %in% params_list)] = NULL
    d = main(params)
    d = cbind(params,d)
    write.table(d,col.names=FALSE,sep=",",quote=FALSE,row.names=FALSE)
}

