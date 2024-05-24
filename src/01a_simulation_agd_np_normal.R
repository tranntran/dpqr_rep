# This code is 1 out of 3 files that compare the accuracy of different DP methods
# in the simulation setting with normal noise. For each replication, we estimate
# the coefficients on a newly simulated dataset and calculate the L2 difference to
# the true coefficients.
# This file focuses only on the DP adaptive gradient descent (DP-AGD) by Lee &
# Kifer (2018) together with check loss smoothed by convolution method as well
# as the non-private version from package quantreg.
# The other two file has name simulation_[method]_normal.R.
# The output is stored in folder output/simulation_normalnoise and it can be plotted 
# using file plot_simulation_normal.R.
# Last edited: 3/18/24

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

# function for noisy max algorithm
# return index of the largest values in a list after adding random laplace noise
noisy_max = function(candidates, delta_f, epsilon) {
    candidates = candidates + rmutil::rlaplace(n = length(candidates), m = 0, s = delta_f/epsilon)
    return(which.max(candidates))
}


# function for finding the average gradient descent (Algorithm 2 in Lee & Kifer, 2018)
# used after the adaptive GD fails to find the descent direction and triggers 
# an increase in privacy budget required for the update
grad_avg = function(rho_old, rho_h, g, g_tilde, c_grad) {
    n = length(g)
    mu = rep(0, n)
    sgm = c_grad^2 / (2*(rho_h - rho_old)) * diag(n)
    g_tilde2 = g + MASS::mvrnorm(n = 1, mu = mu, Sigma = sgm)
    s_tilde = (rho_old*g_tilde + (rho_h - rho_old)*g_tilde2) / rho_h
    return(s_tilde)
}


# function for converting (epsilon, delta) privacy into rho-zCDP
# solve and check validity using equation 5 in Lee & Kifer, 2018
find_rho = function(epsilon_tot, delta_tot, tol = 1e-10) {
    rho = (-sqrt(log(1/delta_tot)) + sqrt(log(1/delta_tot) + epsilon_tot))^2
    check = abs(epsilon_tot - rho - 2*sqrt(rho*log(1/delta_tot))) <= tol
    
    while (!check) {
        epsilon_tot = epsilon_tot*1.01
        rho = -sqrt(log(1/delta_tot)) + sqrt(log(1/delta_tot) + epsilon_tot)
        check = (epsilon_tot - rho - 2*sqrt(rho*log(1/delta_tot))) >= 0
    }
    return(rho)
}

# function for differentially private adaptive gradient descent 
dp_agd = function(x, y, tau, h, rho_nmax, rho_ng, epsilon_tot, delta_tot, gamma,
                  c_obj, c_grad, step_increase_rate, step_max_range, step_num, 
                  step_update_iter, kernel = c('gaussian', 'logistic', 'uniform')) {
    n = nrow(x)
    p = ncol(x)
    t = 1
    
    # initialize theta
    theta_t = matrix(0, p)
    
    # initialize list of possible step sizes, determined by step_max_range / step_num
    # this list will update every step_update_iter using 
    # new_step_max_range = (1 + step_increase_rate)*max(actual step size since last update)
    phi = seq(0, step_max_range, length.out = step_num)
    # initialize variable to save actual step size used
    actual_steps = NULL 
    
    # convert (epsilon, delta)-DP to rho-zCDP
    rho = find_rho(epsilon_tot = epsilon_tot, delta_tot = delta_tot)
    
    kernel = match.arg(kernel)
    kernel_info = select_kernel(kernel)
    integral_fn = kernel_info[[2]]
    cdf_k_h = kernel_info[[3]]
    
    loss_gradient = function(theta, x, y, tau, h) {
        ans = (cdf_k_h(x%*%theta - y, h) - tau) %*% x
        return(ans)
    }
    
    f = function(theta, x, y, tau, h) {
        ans = 0
        for (i in 1:nrow(x)) {
            u = y[i] - x[i, ]%*%theta
            ans = ans + h/2*integral_fn(u, h) + (tau - 1/2)*u
        }
        return(ans)
    }
    
    while (rho > 0) {
        print(t)
        i = 1
        g_t = 0
        for (j in 1:n) {
            tmp = loss_gradient(theta = as.matrix(theta_t[, t]), x = x[j, ], y = y[j], 
                                tau = tau, h = h)
            g_t = g_t + tmp / max(1, norm(tmp, "F") / c_grad)
        }
        g_t_tilde = g_t + MASS::mvrnorm(1, mu = rep(0, p), c_grad^2/(2*rho_ng)*diag(p))
        rho = rho - rho_ng
        g_t_tilde = g_t_tilde / norm(g_t_tilde, "F")
        
        while (i == 1) {
            omega = NULL
            for (a in phi) { # double check dim of g_t_tilde
                tmp = matrix(theta_t[, t] - a*g_t_tilde, ncol = 1)
                omega = c(omega, f(theta = tmp, x = x, y = y, tau = tau, h = h))
            }
            
            rho = rho - rho_nmax
            i = noisy_max(candidates = -omega, delta_f = c_obj, epsilon = sqrt(2*rho_nmax))
            print(i)
            if (i > 1) {
                if (rho > 0) {
                    theta_t = cbind(theta_t, t(theta_t[, t] - phi[i]*g_t_tilde))
                    actual_steps = c(actual_steps, phi[i])
                }
            } else {
                rho_old = rho_ng
                rho_ng = (1 + gamma)*rho_ng
                g_t_tilde = grad_avg(rho_old = rho_old, rho_h = rho_ng, g = g_t, 
                                     g_tilde = g_t_tilde, c_grad = c_grad)
                rho = rho - (rho_ng - rho_old)
                if (rho <= 0) break
                print(paste("rho", rho))
            }
        }
        
        if (t %% step_update_iter == 0) {
            step_max_range = (1 + step_increase_rate)*max(actual_steps)
            phi = seq(0, step_max_range, length.out = step_num)
            actual_steps = NULL
        }
        
        t = t+1
        
    }
    
    return(list("par" = tail(t(theta_t), 1), "all" = theta_t))
}

l2_diff = function(diff) {
    return(sqrt(sum(diff)^2))
}

dp_qr = function(x, y, theta, taus, epsilon, delta = 1e-8, Cx) { #nsteps = 100000, Cx = 1) {
    out_diff = NULL
    out_time = matrix(NA, ncol = 7, nrow = 6)
    x[x > Cx] = Cx
    for (t in 1:length(taus)) {
        tau = taus[t]
        print(paste('Tau', tau))
        out_par = rep(list(NULL), 7)
        true_theta = c(theta)
        true_theta[1] = true_theta[1] + qnorm(tau, 0, 2)
        
        start_points = matrix(0, ncol = 3, nrow = 4)
        
        rho_nmax = find_rho(epsilon/120, delta)
        rho_ng = rho_nmax
        
        out_time[t, 1] = system.time({tmp = dp_agd(x = x, y = y, tau = tau, h = 0.5, kernel = 'gaussian',
                                                   rho_nmax = rho_nmax, rho_ng = rho_ng, 
                                                   epsilon_tot = epsilon, delta_tot = delta, 
                                                   gamma = 0.3, c_obj = 15, c_grad = 15, 
                                                   step_increase_rate = 0.1, step_max_range = 2, 
                                                   step_num = 20, step_update_iter = 10)$par})[3]
        out_par[[1]] = rbind(out_par[[1]], tmp)
        
        out_time[t, 2] = system.time({tmp = dp_agd(x = x, y = y, tau = tau, h = 3, kernel = 'gaussian',
                                                   rho_nmax = rho_nmax, rho_ng = rho_ng, 
                                                   epsilon_tot = epsilon, delta_tot = delta, 
                                                   gamma = 0.3, c_obj = 15, c_grad = 15, 
                                                   step_increase_rate = 0.1, step_max_range = 2, 
                                                   step_num = 20, step_update_iter = 10)$par})[3]
        out_par[[2]] = rbind(out_par[[2]], tmp)
        
        out_time[t, 3] = system.time({tmp = dp_agd(x = x, y = y, tau = tau, h = 0.5, kernel = 'logistic',
                                                   rho_nmax = rho_nmax, rho_ng = rho_ng, 
                                                   epsilon_tot = epsilon, delta_tot = delta, 
                                                   gamma = 0.3, c_obj = 15, c_grad = 15, 
                                                   step_increase_rate = 0.1, step_max_range = 2, 
                                                   step_num = 20, step_update_iter = 10)$par})[3]
        out_par[[3]] = rbind(out_par[[3]], tmp)
        
        out_time[t, 4] = system.time({tmp = dp_agd(x = x, y = y, tau = tau, h = 3, kernel = 'logistic',
                                                   rho_nmax = rho_nmax, rho_ng = rho_ng, 
                                                   epsilon_tot = epsilon, delta_tot = delta, 
                                                   gamma = 0.3, c_obj = 15, c_grad = 15, 
                                                   step_increase_rate = 0.1, step_max_range = 2, 
                                                   step_num = 20, step_update_iter = 10)$par})[3]
        out_par[[4]] = rbind(out_par[[4]], tmp)
        
        out_time[t, 5] = system.time({tmp = dp_agd(x = x, y = y, tau = tau, h = 0.5, kernel = 'uniform',
                                                   rho_nmax = rho_nmax, rho_ng = rho_ng, 
                                                   epsilon_tot = epsilon, delta_tot = delta, 
                                                   gamma = 0.3, c_obj = 15, c_grad = 15, 
                                                   step_increase_rate = 0.1, step_max_range = 2, 
                                                   step_num = 20, step_update_iter = 10)$par})[3]
        out_par[[5]] = rbind(out_par[[5]], tmp)
        
        out_time[t, 6] = system.time({tmp = dp_agd(x = x, y = y, tau = tau, h = 3, kernel = 'uniform',
                                                   rho_nmax = rho_nmax, rho_ng = rho_ng, 
                                                   epsilon_tot = epsilon, delta_tot = delta, 
                                                   gamma = 0.3, c_obj = 15, c_grad = 15, 
                                                   step_increase_rate = 0.1, step_max_range = 2, 
                                                   step_num = 20, step_update_iter = 10)$par})[3]
        out_par[[6]] = rbind(out_par[[6]], tmp)
        
        out_time[t, 7] = system.time({out_par[[7]] = rbind(out_par[[7]], coef(rq(y ~ x[,-1], tau = tau)))})[3]
        
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
                   # create training and testing data with 5000 obs each
                   # x1 = unif(0, 1)
                   # x2 = unif(0, 1)
                   # x3 = 3 + 2x1 + 0.5x2 + exp(0.1)
                   # remember to update c_obj and c_grad based on Cx or rescale x
                   n = 2500  
                   x1 = rnorm(n, 0, 3)
                   x2 = rnorm(n, 0, 1)
                   x = cbind(1, x1, x2)
                   theta = matrix(c(3, 2, -1))
                   y = x%*% theta + rnorm(n, 0, 2)
                   
                   out = dp_qr(x = x, y = y, theta = theta, epsilon = 1, Cx = 10,
                               taus = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) #change
               }

stopCluster(workers)

save(oper, file = paste0("./output/simulation_normalnoise/agd_np_n2500_eps1_", iter, ".Rdata"))

