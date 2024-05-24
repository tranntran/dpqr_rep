# This code is 1 out of 3 files that compare the accuracy of different DP methods
# in the simulation setting with normal noise. For each replication, we estimate
# the coefficients on a newly simulated dataset and calculate the L2 difference to
# the true coefficients.
# This file focuses only on KNG + convolution smoothing method, check loss, and 
# check loss with regularization. The other two file has name simulation_[method]_normal.R.
# The output is stored in folder output/simulation_normalnoise and it can be plotted 
# using file plot_simulation_normal.R.
# Last edited: 2/8/24

kng_qr = function(theta, x, y, tau, epsilon = 1, Cx = 1) {
    n = nrow(x) 
    tmp = 0
    for (i in 1:n) {
        tmp = tmp - tau%*%x[i, ] + I(y[i]<= x[i, ]%*%theta)%*%x[i,]
    }
    delta = 2*(1-tau)*Cx
    ans = -epsilon / (2*delta) * norm(tmp, "M")
    return(ans)
}

kng_qr_reg = function(theta, x, y, tau, epsilon = 1, Cx = 1) {
    n = nrow(x) 
    tmp = 0
    for (i in 1:n) {
        tmp = tmp - tau%*%x[i, ] + I(y[i]<= x[i, ]%*%theta)%*%x[i,]
    }
    delta = 2*(1-tau)*Cx
    ans = -epsilon / (2*delta) * norm(tmp, "M") - 1/2*theta%*%theta
    return(ans)
}

kng_qr_gaussian = function(theta, x, y, tau, h, epsilon = 1, Cx = 1) {
    n = nrow(x)
    tmp = 0
    for (i in 1:n) {
        tmp = tmp + (pnorm((x[i, ] %*% theta - y[i])/h) - tau) %*% x[i, ]
    }
    delta = 2*max(1-tau, tau)*Cx
    ans = -epsilon / (2*delta) * norm(tmp, "M")
    return(ans)
}

kng_qr_gaussian_del = function(theta, x, y, tau, h, epsilon = 1, Cx = 1) {
    n = nrow(x)
    tmp = 0
    for (i in 1:n) {
        tmp = tmp + (pnorm((x[i, ] %*% theta - y[i])/h) - tau) %*% x[i, ]
    }
    delta = 2*max(1-tau, tau)*Cx
    k_bar = 1/sqrt(2*pi)
    gamma = n*k_bar/h
    lambda = n*k_bar*Cx^2/h
    del = 2*lambda/epsilon
    ans = -epsilon / (2*delta) * norm(tmp, "M") - max(del - gamma, 0)*norm(as.matrix(theta), "F")^2/(2*n)
    return(ans)
}

kng_qr_logistic = function(theta, x, y, tau, h, epsilon = 1, Cx = 1) {
    n = nrow(x)
    tmp = 0
    for (i in 1:n) {
        exp_val = exp(-(x[i,] %*% theta - y[i])/h)
        tmp = tmp + (1/(1 + exp_val) - tau) %*% x[i, ]
    }
    delta = 2*max(1-tau, tau)*Cx
    ans = -epsilon / (2*delta) * norm(tmp, "M")
    return(ans)   
}

kng_qr_logistic_del = function(theta, x, y, tau, h, epsilon = 1, Cx = 1) {
    n = nrow(x)
    tmp = 0
    for (i in 1:n) {
        exp_val = exp(-(x[i,] %*% theta - y[i])/h)
        tmp = tmp + (1/(1 + exp_val) - tau) %*% x[i, ]
    }
    delta = 2*max(1-tau, tau)*Cx
    k_bar = 1/4
    gamma = n*k_bar/h
    lambda = n*k_bar*Cx^2/h
    del = 2*lambda/epsilon
    ans = -epsilon / (2*delta) * norm(tmp, "M") - max(del - gamma, 0)*norm(as.matrix(theta), "F")^2/(2*n)
    return(ans)   
}

kng_qr_uniform = function(theta, x, y, tau, h, epsilon = 1, Cx = 1) {
    n = nrow(x)
    tmp = 0
    for (i in 1:n) {
        tmp_val =  min(((x[i,] %*% theta - y[i])/h + 1)/2, 1)*((x[i,] %*% theta - y[i])/h >= -1)
        tmp = tmp + (tmp_val - tau) %*% x[i, ]
    }
    delta = 2*max(1-tau, tau)*Cx
    ans = -epsilon / (2*delta) * norm(tmp, "M")
    return(ans)   
}

kng_qr_uniform_del = function(theta, x, y, tau, h, epsilon = 1, Cx = 1) {
    n = nrow(x)
    tmp = 0
    for (i in 1:n) {
        tmp_val =  min(((x[i,] %*% theta - y[i])/h + 1)/2, 1)*((x[i,] %*% theta - y[i])/h >= -1)
        tmp = tmp + (tmp_val - tau) %*% x[i, ]
    }
    delta = 2*max(1-tau, tau)*Cx
    k_bar = 1/2
    gamma = n*k_bar/h
    lambda = n*k_bar*Cx^2/h
    del = 2*lambda/epsilon
    ans = -epsilon / (2*delta) * norm(tmp, "M") - max(del - gamma, 0)*norm(as.matrix(theta), "F")^2/(2*n)
    return(ans)   
}

l2_diff = function(diff) {
    return(sqrt(sum(diff)^2))
}

dp_qr2 = function(x, y, theta, taus, epsilon, delta = 1e-8, nsteps = 100000, Cx = 1) {
    out_diff = NULL
    out_time = matrix(NA, ncol = 8, nrow = 6)
    x[x > Cx] = Cx
    for (t in 1:length(taus)) {
        tau = taus[t]
        print(paste('Tau', tau))
        out_par = rep(list(NULL), 8)
        true_theta = c(theta)
        true_theta[1] = true_theta[1] + qnorm(tau, 0, 2)
        
        start_points = matrix(0, ncol = 3, nrow = 4)
        
        out_time[t, 1] = system.time({tmp = MCMC(start_points, kng_qr_gaussian_del, 
                                                 nsteps = nsteps, nchains = 4, x = x,  Cx = Cx,
                                                 y = y, tau = tau, epsilon = epsilon, 
                                                 h = 0.5, conv_checker = convergence_gelman(), 
                                                 kernel = kernel_ram())})[3]
        out_par[[1]] = rbind(out_par[[1]], tail(as.matrix(tmp[[1]]), 1))
        
        out_time[t, 2] = system.time({tmp = MCMC(start_points, kng_qr_gaussian_del, 
                                                 nsteps = nsteps, nchains = 4, x = x,  Cx = Cx,
                                                 y = y, tau = tau, epsilon = epsilon, 
                                                 h = 3, conv_checker = convergence_gelman(), 
                                                 kernel = kernel_ram())})[3]
        out_par[[2]] = rbind(out_par[[2]], tail(as.matrix(tmp[[1]]), 1))
        
        out_time[t, 3] = system.time({tmp = MCMC(start_points, kng_qr_logistic_del, 
                                                 nsteps = nsteps, nchains = 4, x = x,  Cx = Cx,
                                                 y = y, tau = tau, epsilon = epsilon, 
                                                 h = 0.5, conv_checker = convergence_gelman(), 
                                                 kernel = kernel_ram())})[3]
        out_par[[3]] = rbind(out_par[[3]], tail(as.matrix(tmp[[1]]), 1))
        
        out_time[t, 4] = system.time({tmp = MCMC(start_points, kng_qr_logistic_del, 
                                                 nsteps = nsteps, nchains = 4, x = x,  Cx = Cx,
                                                 y = y, tau = tau, epsilon = epsilon, 
                                                 h = 3, conv_checker = convergence_gelman(), 
                                                 kernel = kernel_ram())})[3]
        out_par[[4]] = rbind(out_par[[4]], tail(as.matrix(tmp[[1]]), 1))
        
        out_time[t, 5] = system.time({tmp = MCMC(start_points, kng_qr_uniform_del, 
                                                 nsteps = nsteps, nchains = 4, x = x,  Cx = Cx,
                                                 y = y, tau = tau, epsilon = epsilon, 
                                                 h = 0.5, conv_checker = convergence_gelman(), 
                                                 kernel = kernel_ram())})[3]
        out_par[[5]] = rbind(out_par[[5]], tail(as.matrix(tmp[[1]]), 1))
        
        out_time[t, 6] = system.time({tmp = MCMC(start_points, kng_qr_uniform_del, 
                                                 nsteps = nsteps, nchains = 4, x = x,  Cx = Cx,
                                                 y = y, tau = tau, epsilon = epsilon, 
                                                 h = 3, conv_checker = convergence_gelman(), 
                                                 kernel = kernel_ram())})[3]
        out_par[[6]] = rbind(out_par[[6]], tail(as.matrix(tmp[[1]]), 1))
        
        out_time[t, 7] = system.time({tmp = MCMC(start_points, kng_qr, nsteps = nsteps, 
                                                 nchains = 4, x = x, y = y, tau = tau, Cx = Cx,
                                                 epsilon = epsilon, conv_checker = convergence_gelman(), 
                                                 kernel = kernel_ram())})[3]
        out_par[[7]] = rbind(out_par[[7]], tail(as.matrix(tmp[[1]]), 1))
        
        out_time[t, 8] = system.time({tmp = MCMC(start_points, kng_qr_reg, nsteps = nsteps, 
                                                 nchains = 4, x = x, y = y, tau = tau, Cx = Cx,
                                                 epsilon = epsilon, conv_checker = convergence_gelman(), 
                                                 kernel = kernel_ram())})[3]
        out_par[[8]] = rbind(out_par[[8]], tail(as.matrix(tmp[[1]]), 1))
        
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
                   library(fmcmc)
                   library(quantreg)
                   set.seed(i)
                   n = 2500  
                   x1 = rnorm(n, 0, 3)
                   x2 = rnorm(n, 0, 1)
                   x = cbind(1, x1, x2)
                   theta = matrix(c(3, 2, -1))
                   y = x%*% theta + rnorm(n, 0, 2)
                   
                   out = dp_qr2(x = x, y = y, theta = theta, epsilon = 1, Cx = 10,
                                taus = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) #change
               }

stopCluster(workers)

save(oper, file = paste0("./output/simulation_normalnoise/kng_n2500_eps1_", iter, ".Rdata"))

