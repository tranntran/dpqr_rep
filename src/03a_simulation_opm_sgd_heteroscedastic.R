# This code is 1 out of 3 files that compare the accuracy of different DP methods
# in the simulation setting with heteroskedastic data. For each replication, we 
# estimate the coefficients on a newly simulated dataset and calculate the L2 
# difference to the true coefficients.
# This is an implementation of the algorithm proposed by Chen and Chua (2023),
# which uses the convolution smoothing of quantile loss within OPM and SGD. 
# The other two file has name simulation_[method]_heteroscedastic.R.
# The output is stored in folder output/simulation_heteroscedastic and it can be plotted 
# using file plot_simulation_heteroscedastic.R.
# Last edited: 3/20/24

select_kernel = function(kernel = c('gaussian', 'logistic', 'uniform')) {
    kernel = match.arg(kernel)
    if (kernel == 'gaussian') {
        k_bar = 1/sqrt(2*pi)
        integral_fn = function(u, h) {
            ans = sqrt(2/pi)*exp(-(u/h)^2/2) +(u/h)*(1-2*pnorm(-u/h))
            return(ans)
        }
        cdf_k_h = function(u, h) {
            return(pnorm(u/h))
        }
    } else if (kernel == 'logistic') {
        k_bar = 1/4
        integral_fn = function(u, h) {
            ans = u/h + 2*log(1+ exp(-u/h))
            return(ans)
        }
        cdf_k_h = function(u, h) {
            return(1/(1 + exp(-u/h)))
        }
    } else if (kernel == 'uniform') {
        k_bar = 1/2
        integral_fn = function(u, h) {
            ans = ifelse(abs(u/h) <= 1, u^2/(2*h^2) + 1/2, abs(u/h))
            return(ans)
        }
        cdf_k_h = function(u, h) {
            return(min((u/h + 1)/2, 1)*(u/h >= -1))
        }
    }
    return(list(k_bar, integral_fn, cdf_k_h))
}

opm_convsm = function(theta, x, y, tau, h, epsilon = 1, Cx, 
                      kernel = c('gaussian', 'logistic', 'uniform')) {
    kernel = match.arg(kernel)
    kernel_info = select_kernel(kernel)
    k_bar = kernel_info[[1]]
    integral_fn = kernel_info[[2]]
    
    theta = matrix(theta, ncol = 1)
    n = nrow(x)
    p = ncol(x)
    xi = n*max(tau, 1-tau)*Cx
    lambda = n*k_bar*Cx^2/h
    ans = 0
    for (i in 1:n) {
        u = y[i] - x[i, ]%*%theta
        ans = ans + h/2*integral_fn(u, h) + (tau -1/2)*u
    }
    delta = 2*lambda/epsilon
    gamma = n*k_bar/h # parameter for strong convexity
    
    # sample from Gamma distribution to satisfy epsilon-DP
    norm = rgamma(n = 1, shape = p, rate = epsilon/(2*xi))
    vector = rnorm(n = p, m = 0, s = 1)
    b = vector*norm/sqrt(sum(vector^2))
    
    ans = ans + (delta - gamma)*norm(theta, "F")^2/(2*n) + b%*%theta/n
    
    return(ans)
}

sgd_convsm = function(theta, x, y, tau, h, epsilon = 1, Cx, eta,
                      kernel = c('gaussian', 'logistic', 'uniform')) {
    kernel = match.arg(kernel)
    kernel_info = select_kernel(kernel)
    k_bar = kernel_info[[1]]
    cdf_k_h = kernel_info[[3]]
    
    L = max(tau, 1-tau)*Cx
    beta = k_bar*Cx^2/h
    
    n = nrow(x)
    p = ncol(x)
    theta_all = matrix(0, nrow = n^2, ncol = p)
    theta_all[1, ] = theta
    
    for(i in 1:(n^2-1)) {
        samp = sample(1:n, 1)
        x_t = x[samp, ]
        y_t = y[samp]
        grad = (cdf_k_h(c(x_t%*%as.matrix(theta_all[i,]) - y_t), h) - tau) * x_t
        
        # sample from Gamma distribution to satisfy epsilon-DP
        # norm = rgamma(n = 1, shape = p, rate = epsilon/(2*L))
        # vector = rnorm(n = p, m = 0, s = 1)
        # b = vector*norm/sqrt(sum(vector^2))
        
        delta = 1e-6
        sigma2 = 8*L^2*log(1/delta)/epsilon^2
        w = MASS::mvrnorm(n = 1, mu = rep(0, p), Sigma = sigma2*diag(p))
        
        theta_all[i+1, ] = theta_all[i, ] - 1/sqrt(i)*(grad + w)
    }
    
    ans = apply(theta_all, 2, mean)
    return(ans)
}


l2_diff = function(diff) {
    return(sqrt(sum(diff)^2))
}

# solution from stackoverflow
#https://stats.stackexchange.com/questions/14481/quantiles-from-the-combination-of-normal-distributions

F = function(x,w,u,s) sum(w*pnorm(x, mean=u, sd=s))
# provide an initial bracket for the quantile. default is c(-1000,1000). 
F_inv = function(p, w, u, s, br=c(-1000,1000)) {
    G = function(x) F(x, w, u, s) - p
    return(uniroot(G,br)$root) 
    }


dp_qr3 = function(x, y, theta, taus, epsilon, eta = 0.01, Cx = 1) {
    out_diff = NULL
    out_time = matrix(NA, ncol = 12, nrow = 6)
    x[x > Cx] = Cx
    p = ncol(x)
    for (t in 1:length(taus)) {
        tau = taus[t]
        print(paste('Tau', tau))
        out_par = rep(list(NULL), 12)
        true_theta = c(theta)
        true_theta[1] = true_theta[1] + F_inv(tau, rep(1/n, n), rep(0, n), x[, 2])
        
        out_time[t, 1] = system.time({tmp = optim(rep(0, p), opm_convsm, x = x, y = y,  
                                                  Cx = Cx, tau = tau, h = 0.5, 
                                                  epsilon = epsilon, kernel = 'gaussian')$par})[3]
        out_par[[1]] = rbind(out_par[[1]], tmp)
        
        out_time[t, 2] = system.time({tmp = optim(rep(0, p), opm_convsm, x = x, y = y,  
                                                  Cx = Cx, tau = tau, h = 3, 
                                                  epsilon = epsilon, kernel = 'gaussian')$par})[3]
        out_par[[2]] = rbind(out_par[[2]], tmp)
        
        out_time[t, 3] = system.time({tmp = optim(rep(0, p), opm_convsm, x = x, y = y,  
                                                  Cx = Cx, tau = tau, h = 0.5, 
                                                  epsilon = epsilon, kernel = 'logistic')$par})[3]
        out_par[[3]] = rbind(out_par[[3]], tmp)
        
        out_time[t, 4] = system.time({tmp = optim(rep(0, p), opm_convsm, x = x, y = y,  
                                                  Cx = Cx, tau = tau, h = 3, 
                                                  epsilon = epsilon, kernel = 'logistic')$par})[3]
        out_par[[4]] = rbind(out_par[[4]], tmp)
        
        out_time[t, 5] = system.time({tmp = optim(rep(0, p), opm_convsm, x = x, y = y,  
                                                  Cx = Cx, tau = tau, h = 0.5, 
                                                  epsilon = epsilon, kernel = 'uniform')$par})[3]
        out_par[[5]] = rbind(out_par[[5]], tmp)
        
        out_time[t, 6] = system.time({tmp = optim(rep(0, p), opm_convsm, x = x, y = y,  
                                                  Cx = Cx, tau = tau, h = 3, 
                                                  epsilon = epsilon, kernel = 'uniform')$par})[3]
        out_par[[6]] = rbind(out_par[[6]], tmp)
        
        out_time[t, 7] = system.time({tmp = sgd_convsm(theta = rep(0, p), x = x, y = y, 
                                                       tau = tau, h = 0.5, epsilon = epsilon, 
                                                       Cx = Cx, eta = eta, kernel = 'gaussian')})[3]
        out_par[[7]] = rbind(out_par[[7]], tmp)
        
        out_time[t, 8] = system.time({tmp = sgd_convsm(theta = rep(0, p), x = x, y = y, 
                                                       tau = tau, h = 3, epsilon = epsilon, 
                                                       Cx = Cx, eta = eta, kernel = 'gaussian')})[3]
        out_par[[8]] = rbind(out_par[[8]], tmp)
        
        out_time[t, 9] = system.time({tmp = sgd_convsm(theta = rep(0, p), x = x, y = y, 
                                                       tau = tau, h = 0.5, epsilon = epsilon, 
                                                       Cx = Cx, eta = eta, kernel = 'logistic')})[3]
        out_par[[9]] = rbind(out_par[[9]], tmp)
        
        out_time[t, 10] = system.time({tmp = sgd_convsm(theta = rep(0, p), x = x, y = y, 
                                                        tau = tau, h = 3, epsilon = epsilon, 
                                                        Cx = Cx, eta = eta, kernel = 'logistic')})[3]
        out_par[[10]] = rbind(out_par[[10]], tmp)
        
        out_time[t, 11] = system.time({tmp = sgd_convsm(theta = rep(0, p), x = x, y = y, 
                                                        tau = tau, h = 0.5, epsilon = epsilon, 
                                                        Cx = Cx, eta = eta, kernel = 'uniform')})[3]
        out_par[[11]] = rbind(out_par[[11]], tmp)
        
        out_time[t, 12] = system.time({tmp = sgd_convsm(theta = rep(0, p), x = x, y = y, 
                                                        tau = tau, h = 3, epsilon = epsilon, 
                                                        Cx = Cx, eta = eta, kernel = 'uniform')})[3]
        out_par[[12]] = rbind(out_par[[12]], tmp)
        
        
        out_diff = rbind(out_diff, sapply(out_par, function(l) l2_diff(l - true_theta)))
        
    }
    
    return(list(out_diff, out_time))
    
}


library(doParallel)
num_cores=detectCores()-1 # use all available core but 1
workers=makeCluster(num_cores,type="SOCK",outfile="log.txt")
registerDoParallel(workers)


iter = as.numeric(commandArgs(trailingOnly = TRUE))
start = 25*(iter-1) + 1
end = 25*iter 


# run 100 reps, the ith rep is run with seed i
oper = foreach(i=start:end, .combine=rbind, .multicombine=TRUE, #change
               .init=list()) %dopar% {
                   set.seed(i)
                   n = 2500  
                   x1 = runif(n, 0, 2)
                   x = cbind(1, x1)
                   theta = matrix(c(1, 2))
                   y = x%*% theta + x1*rnorm(n, 0, 1)
                   
                   out = dp_qr3(x = x, y = y, theta = theta, epsilon = 1, Cx = 2,
                                taus = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) #change
               }

stopCluster(workers)

save(oper, file = paste0("./output/simulation_heteroscedastic/opm_sgd_n2500_eps1_", iter, ".Rdata"))



