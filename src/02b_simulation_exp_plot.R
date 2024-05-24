# This code is used to plot the output from the simulation with exponential noise 
# (i.e., output from file simulation_[method]_exp.R).
# Last updated: 01/28/2024

library(tidyr)
library(ggplot2)
library(wesanderson)
library(ggpubr)

###### COMPARE PERFORMANCE OF THE BEST OF EACH METHOD ###### 
# This code can be used to plot the performance of the best kernel and bandwidth
# combination by each method/mech (KNG, OPM, DPSGD, and DPAGD)

plot_best_methods = function(n, eps) {
    out = rep(list(NULL), 27)
    out_time = rep(list(NULL), 27)
    
    if (n == 2500) {
        for (k in 1:20) {
            file_name = paste0("output/simulation_expnoise/kng_del_n", n, "_eps", eps, "_", k, ".Rdata")
            load(file_name)
            for (i in 1:25) {
                tmp = t(oper[[i]])
                for (j in 1:8) {
                    out[[j]] = rbind(out[[j]], tmp[j,])
                }
            }
        }
    } else if (n == 250) {
        for (k in 1:100) {
            file_name = paste0("output/simulation_expnoise/kng_del_n", n, "_eps", eps, "_", k, ".Rdata")
            load(file_name)
            for (i in 1:5) {
                tmp = t(oper[[i]])
                for (j in 1:8) {
                    out[[j]] = rbind(out[[j]], tmp[j,])
                }
            }
        }
        
    }
    
    for (k in 1:20) {
        file_name = paste0("output/simulation_expnoise/opm_sgd_n", n, "_eps", eps, "_", k, ".Rdata")
        load(file_name)
        for (i in 1:25) {
            tmp = t(oper[[i]])
            for (j in 1:12) {
                out[[j+8]] = rbind(out[[j+8]], tmp[j,])
            }
        }
    }
    
    for (k in 1:20) {
        file_name = paste0("output/simulation_expnoise/agd_np_n", n, "_eps", eps, "_", k, ".Rdata")
        load(file_name)
        for (i in 1:25) {
            tmp = t(oper[[i]])
            for (j in 1:7) {
                out[[j+20]] = rbind(out[[j+20]], tmp[j,])
            }
        }
    }
    
    out = lapply(out, log10)
    
    # find the mean/ median L2 distance to the truth
    tab = t(sapply(out, function(x) apply(x, 2, mean)))
    
    # prepare data for plotting
    tab = as.data.frame(tab)
    taus = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)#*100
    colnames(tab) = taus
    method_tmp = as.vector(t(outer(c('Gaussian', 'Logistic', 'Uniform'), c('(h=0.5)', '(h=3)'), paste)))
    tab$Method = c(paste('KNG', method_tmp), 'KNG Check Loss', 'KNG Check w/ Reg', 
                   paste('OPM', method_tmp), paste('DPSGD', method_tmp), 
                   paste('DPAGD', method_tmp), 'Non-Private')
    
    
    tab_best_best = tab[c(3, 11, 17, 23, 27),]
    tab_best_best$Method  = factor(tab_best_best$Method, 
                                   levels = c("Non-Private", 'KNG Logistic (h=0.5)',
                                              'OPM Logistic (h=0.5)', 'DPSGD Logistic (h=0.5)',
                                              'DPAGD Logistic (h=0.5)'))
    
    tab_best_best = gather(tab_best_best, Quantile, Value, 1:6)
    tab_best_best$Quantile = as.numeric(tab_best_best$Quantile) 
    
    colors_bestbest = c("#000000", "#3B9AB2", "#BDC881", "#E3B710", "#F11B00")
    names(colors_bestbest) = levels(tab_best_best$Method)
    my_linewidth = c(1.25, rep(1, 5))
    names(my_linewidth) = levels(tab_best_best$Method)
    
    yl = ifelse(n == 2500, 3, 3)
    
    e = ifelse(eps == '1', 1, 0.5)
    
    pl = ggplot(tab_best_best, aes(x = Quantile, y = Value, group = Method)) +
        geom_line(aes(color = Method, size = Method)) + 
        coord_cartesian(ylim=c(-1.5, yl)) + 
        labs(subtitle =  substitute(paste("n=", n, ", ", epsilon, "=", e, sep = '')), 
             y = expression(log[10]*"(||"*hat(beta)-beta*"||"[2]*")"))  +
        scale_x_continuous(breaks=taus, labels=taus) + 
        scale_color_manual(values = colors_bestbest) + 
        scale_size_manual(values = my_linewidth)  +
        theme_minimal() +
        theme(axis.text.x=element_text(size=13, angle=90), axis.text.y = element_text(size=13),
              axis.title = element_text(size=15), plot.subtitle = element_text(size=16), 
              legend.text = element_text(size = 14), legend.title = element_text(size = 15))
    return(list(pl, tab_best_best))
    
}



time_best_methods = function(n, eps) {
    out = rep(list(NULL), 27)
    if (n == 2500) {
        for (k in 1:20) {
            file_name = paste0("output/simulation_expnoise/kng_del_n", n, "_eps", eps, "_", k, ".Rdata")
            load(file_name)
            for (i in 26:50) {
                tmp = t(oper[[i]])
                for (j in 1:8) {
                    out[[j]] = rbind(out[[j]], tmp[j,])
                }
            }
        }
    } else if (n == 250) {
        for (k in 1:100) {
            file_name = paste0("output/simulation_expnoise/kng_del_n", n, "_eps", eps, "_", k, ".Rdata")
            load(file_name)
            for (i in 6:10) {
                tmp = t(oper[[i]])
                for (j in 1:8) {
                    out[[j]] = rbind(out[[j]], tmp[j,])
                }
            }
        }
        
    }
    
    for (k in 1:20) {
        file_name = paste0("output/simulation_expnoise/opm_sgd_n", n, "_eps", eps, "_", k, ".Rdata")
        load(file_name)
        for (i in 26:50) {
            tmp = t(oper[[i]])
            for (j in 1:12) {
                out[[j+8]] = rbind(out[[j+8]], tmp[j,])
            }
        }
    }
    
    for (k in 1:20) {
        file_name = paste0("output/simulation_expnoise/agd_np_n", n, "_eps", eps, "_", k, ".Rdata")
        load(file_name)
        for (i in 26:50) {
            tmp = t(oper[[i]])
            for (j in 1:7) {
                out[[j+20]] = rbind(out[[j+20]], tmp[j,])
            }
        }
    }
    
    
    # find the mean/ median L2 distance to the truth
    tab = t(sapply(out, function(x) apply(x, 2, mean)))
    tab_se = t(sapply(out, function(x) apply(x, 2, sd)/sqrt(500)))
    
    # prepare data for plotting
    tab = as.data.frame(tab)
    tab_se = as.data.frame(round(tab_se, 3))
    taus = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)*100
    colnames(tab) = taus
    colnames(tab_se) = taus
    method_tmp = as.vector(t(outer(c('Gaussian', 'Logistic', 'Uniform'), c('(h=0.5)', '(h=3)'), paste)))
    tab$Method = c(paste('KNG', method_tmp), 'KNG Check Loss', 'KNG Check w/ Reg', 
                   paste('OPM', method_tmp), paste('DPSGD', method_tmp), 
                   paste('DPAGD', method_tmp), 'Non-Private')
    tab_se$Method = c(paste('KNG', method_tmp), 'KNG Check Loss', 'KNG Check w/ Reg', 
                   paste('OPM', method_tmp), paste('DPSGD', method_tmp), 
                   paste('DPAGD', method_tmp), 'Non-Private')
    
    
    tab_best_best = tab[c(3, 11, 17, 23, 27),]
    tab_best_best$Method  = factor(tab_best_best$Method, 
                                   levels = c("Non-Private", 'KNG Logistic (h=0.5)',
                                              'OPM Logistic (h=0.5)', 'DPSGD Logistic (h=0.5)',
                                              'DPAGD Logistic (h=0.5)'))
    tab_se_best_best = tab_se[c(3, 11, 17, 23, 27),]
    tab_se_best_best$Method  = factor(tab_best_best$Method, 
                                   levels = c("Non-Private", 'KNG Logistic (h=0.5)',
                                              'OPM Logistic (h=0.5)', 'DPSGD Logistic (h=0.5)',
                                              'DPAGD Logistic (h=0.5)'))
    return(list(tab_best_best, tab_se_best_best))
    
    
}

#https://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots
p1 = plot_best_methods(2500, '1')
t1 = time_best_methods(2500, '1')

p2 = plot_best_methods(2500, '05')
t2 = time_best_methods(2500, '05')

p3 = plot_best_methods(250, '1')
t3 = time_best_methods(250, '1')

p4 = plot_best_methods(250, '05')
t4 = time_best_methods(250, '05')

ggarrange(p1[[1]], p2[[1]], p3[[1]], p4[[1]], ncol=2, nrow=2, common.legend = TRUE, legend="bottom")


save_date = Sys.Date()
save_date = paste0(substring(save_date, 3, 4), substring(save_date, 6, 7),
                   substring(save_date, 9, 10))
save_file = paste0('./plot/', save_date, '_compare_all_exp.png')
png(paste0(save_file), width=875, height=550)
ggarrange(p1[[1]], p2[[1]], p3[[1]], p4[[1]], ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()

###### COMPARE THE FULL PERFORMANCE BY METHOD ###### 

# This code can be used to plot the performance of all kernels and bandwiths
# by each method/mech (KNG, OPM, DPSGD, DPQGD)

rm(list = ls())

# load packages
library(tidyr)
library(ggplot2)
library(wesanderson)
library(ggpubr)

## function to plot all results

plot_all_methods = function(method, n, eps) {
        # combine all output into a list with 9 elements (one for each method)
        out = rep(list(NULL), 7)
        
        if (method == 'kng' & n == 2500) {
            for (k in 1:20) {
                file_name = paste0("output/simulation_expnoise/kng_del_n2500_eps", eps, "_", k, ".Rdata")
                load(file_name)
                for (i in 1:25) {
                    tmp = t(oper[[i]])
                    for (j in 1:6) {
                        out[[j]] = rbind(out[[j]], tmp[j,])
                    }
                }
            }
        } else if (method == 'kng' & n == 250) {
            for (k in 1:100) {
                file_name = paste0("output/simulation_expnoise/kng_del_n250_eps", eps, "_", k, ".Rdata")
                load(file_name)
                for (i in 1:5) {
                    tmp = t(oper[[i]])
                    for (j in 1:6) {
                        out[[j]] = rbind(out[[j]], tmp[j,])
                    }
                }
            }
        } else if (method == 'opm') {
            for (k in 1:20) {
                file_name = paste0("output/simulation_expnoise/opm_sgd_n", n, "_eps", eps, "_", k, ".Rdata")
                load(file_name)
                for (i in 1:25) {
                    tmp = t(oper[[i]])
                    for (j in 1:6) {
                        out[[j]] = rbind(out[[j]], tmp[j,])
                    }
                }
            }
        } else if (method == 'sgd') {
            for (k in 1:20) {
                file_name = paste0("output/simulation_expnoise/opm_sgd_n", n, "_eps", eps, "_", k, ".Rdata")
                load(file_name)
                for (i in 1:25) {
                    tmp = t(oper[[i]])
                    for (j in 7:12) {
                        out[[j-6]] = rbind(out[[j-6]], tmp[j,])
                    }
                }
            }
        } else if (method == 'agd') {
            for (k in 1:20) {
                file_name = paste0("output/simulation_expnoise/agd_np_n", n, "_eps", eps, "_", k, ".Rdata")
                load(file_name)
                for (i in 1:25) {
                    tmp = t(oper[[i]])
                    for (j in 1:6) {
                        out[[j]] = rbind(out[[j]], tmp[j,])
                    }
                }
            }
            
        }
        
        
        for (k in 1:20) {
            file_name = paste0("output/simulation_expnoise/agd_np_n", n, "_eps", eps, "_", k, ".Rdata")
            load(file_name)
            for (i in 1:25) {
                tmp = t(oper[[i]])
                out[[7]] = rbind(out[[7]], tmp[7,])
            }
        }
        
        out = lapply(out,  log10)
        # find the mean/ median L2 distance to the truth
        tab = t(sapply(out, function(x) apply(x, 2, mean)))
        tab_median = t(sapply(out, function(x) apply(x, 2, median)))
        tab_se = t(sapply(out, function(x) apply(x, 2, sd)/sqrt(500)))
        
        # prepare table for SE output
        tab_se = as.data.frame(tab_se)
        taus = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)#*100
        colnames(tab_se) = taus
        method_tmp = as.vector(t(outer(c('Gaussian', 'Logistic', 'Uniform'), c('(h=0.5)', '(h=3)'), paste)))
        
        tab_se$Method = c(method_tmp, "Non-Private") #, 'Check Loss w Reg') 
        tab_se$Method  = factor(tab_se$Method,
                             levels = c('Non-Private', method_tmp)) #, 'Check Loss w Reg') )
        
        # prepare data for plotting
        tab = as.data.frame(tab)
        taus = c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)#*100
        colnames(tab) = taus
        method_tmp = as.vector(t(outer(c('Gaussian', 'Logistic', 'Uniform'), c('(h=0.5)', '(h=3)'), paste)))
        
        tab$Method = c(method_tmp, "Non-Private") #, 'Check Loss w Reg') 
        tab$Method  = factor(tab$Method,
                         levels = c('Non-Private', method_tmp)) #, 'Check Loss w Reg') )
        
        # change data format from wide to long for ggplot
        # using all data
        tab_all = gather(tab, Quantile, Value, 1:6)
        tab_all$Quantile = as.numeric(tab_all$Quantile)
        
        my_colors = c("#000000", "#888888", "#3B9AB2", "#78B7C5", "#BDC881", "#E3B710", "#EC7A05")
        names(my_colors) = levels(tab$Method)
        
        my_linewidth = rep(1, 8)
        
        if (method == 'kng') {ylim = c(-1.5, 2)} 
        else if (method == 'opm') {ylim = c(-1.5, 1.75)}
        else if (method == 'sgd') {ylim = c(-1.5, 2)}
        else if (method == 'agd') {ylim = c(-1.5, 3)}
        
        m = toupper(method)
        e = ifelse(eps == '1', 1, 0.5)
        
        pl = ggplot(tab_all, aes(x = Quantile, y = Value, group = Method)) +
        geom_line(aes(color = Method), size = 1) + 
        coord_cartesian(ylim=ylim) + 
        labs(subtitle =  substitute(paste(m, ", n=", n, ", ", epsilon, "=", e, sep = '')),
             y = expression(log[10]*"(||"*hat(beta)-beta*"||"[2]*")"))  +
        scale_x_continuous(breaks=taus, labels=taus) +
        scale_color_manual(values = my_colors)   +
        theme_minimal() +
        theme(axis.text.x=element_text(size=13, angle=90), axis.text.y = element_text(size=13),
              axis.title = element_text(size=15), plot.subtitle = element_text(size=16), 
              legend.text = element_text(size = 14), legend.title = element_text(size = 15))
        return(list(pl, tab_se))
}


methods = c('kng', 'opm', 'sgd', 'agd')

for (method in methods) {
    print(method)
    p1 = plot_all_methods(method, 2500, '1')
    p2 = plot_all_methods(method, 2500, '05')
    p3 = plot_all_methods(method, 250, '1')
    p4 = plot_all_methods(method, 250, '05')
    
    save_date = Sys.Date()
    save_date = paste0(substring(save_date, 3, 4), substring(save_date, 6, 7),
                       substring(save_date, 9, 10))
    save_file = paste0('./plot/', save_date, '_compare_all_', method, '_exp.png')
    png(paste0(save_file), width=800, height=600)
    print(ggarrange(p1[[1]], p2[[1]], p3[[1]], p4[[1]], ncol=2, nrow=2, common.legend = TRUE, legend="bottom"))
    dev.off()
    
    print(list(p1[[2]], p2[[2]], p3[[2]], p4[[2]]))

}
