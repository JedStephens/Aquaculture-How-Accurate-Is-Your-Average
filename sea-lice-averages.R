#Margin of Error Sea Lice Calculations ----------------------------------------
#Authors: Jed Stephens, Alexes Mes
rm(list = ls())
gc()
library(PracTools)
library(foreach)
library(dplyr)
library(readxl)
library(sf)
library(tidyr)


#Figures
library(plotly)
library(latex2exp)
library(ggplot2)

#-----------------------------------------------
#Section 1: Sampling functions.
#-----------------------------------------------

Vybar_s <- function(n, N, Su2){
  (1-n/N)*(Su2/n) #eq. 3.1, Valliant
}
CV2_ybar_s <- function(n, N, Su2, y_bar_u2){
  (1/n - 1/N)*(Su2/y_bar_u2) #eq. 3.3, Valliant.
}

#------------------------------------------------
#Section 2: High incidence pens
#Poisson
#------------------------------------------------

garwood_poisson_confidence_interval <- function(n, y_bar, alpha){
  k <- n * y_bar #Mean times number of samples.
  
  #Implementing Garwood, xxx
  #Critical values.
  chi_lower = qchisq(alpha/2, 2*k, lower.tail = TRUE)
  chi_upper = qchisq(alpha/2, 2*(k+1), lower.tail = FALSE)
  
  #CI
  CI_lower <- chi_lower/(2*n)
  CI_upper <- chi_upper/(2*n)
  
  #e
  #l = y_bar - e
  #h = y_bar + e
  
  #Results.
  data.frame(n = n,
             y_bar = y_bar,
             CI_lower_bound = CI_lower,
             CI_upper_bound = CI_upper,
             el = y_bar - CI_lower,
             eh = CI_upper - y_bar)
}

#Simulations.
n_vector <- 2:500
#N_vector <- c(50000, 100000, 150000)
N_vector <- c(50000)
#lamda_vector <- seq(from = 2, to = 8, by = 0.5)
lamda_vector <- seq(from = 2, to = 5, by = 0.125)


poisson_experiments <- expand.grid(sample_n=n_vector,
                                  population_N = N_vector,
                                  mean_lice = lamda_vector)

poisson_experiments$Vybar_s <- Vybar_s(poisson_experiments$sample_n,
                                       poisson_experiments$population_N,
                                       poisson_experiments$mean_lice)

poisson_experiments$CV2_ybar_s <- CV2_ybar_s(
  poisson_experiments$sample_n,
  poisson_experiments$population_N,
  poisson_experiments$mean_lice,
  poisson_experiments$mean_lice)

#Determine Garwood.
garwood_poisson <- expand.grid(sample_n = n_vector,
                               mean_lice = lamda_vector)

garwood_poisson <- foreach(r = 1:nrow(garwood_poisson),
                           .combine = dplyr::bind_rows) %do%{
  garwood_poisson_confidence_interval(garwood_poisson[r, "sample_n"],
                                      garwood_poisson[r, "mean_lice"],
                                      0.05)
}

garwood_poisson %>% filter(n == 5)
garwood_poisson %>% filter(n == 10)

filter(poisson_experiments, sample_n == 5)
filter(poisson_experiments, sample_n == 10)

#End of Poisson Analysis ------------------------------------------------------
#Visualise Poisson pens.

gw_matrix <- matrix(data = garwood_poisson$eh, nrow = length(n_vector), ncol = length(lamda_vector))

fig <- plot_ly(garwood_poisson, x = ~lamda_vector, y = ~n_vector, z = ~gw_matrix, type = "surface")
# fig <- fig %>%
#   layout(xaxis = list(title = "Mean sea lice per fish"),
#          yaxis = list(title = "Fish sampled"))
#https://plotly.com/r/reference/#layout-scene-zaxis
fig

#Replaced with GGplot below.
#fig2 <- plot_ly(poisson_experiments, x = ~sample_n, y = ~mean_lice, z = ~CV2_ybar_s*100)
#fig2

#figures/poisson_CV2bar_ys
plot2 <- ggplot(data = poisson_experiments, aes(sample_n, CV2_ybar_s)) +
  geom_line(aes(group = mean_lice)) +
  geom_hline(yintercept = 0.1, color = "blue") + 
  theme_bw() + 
  theme(text=element_text(size=16,  family="serif")) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(xlim = c(0, 30)) + 
  ylab(TeX("${CV^{2}}(\\bar{y}_s)$")) +
  xlab("Sample (n)")

ggsave("plot2.pdf", plot = plot2,
       width = 16, height = 8, units = "cm", dpi = 650)


#There are hardly any changes, therefore we do not include this one.
#fig3 <- plot_ly(poisson_experiments, x = ~sample_n, y = ~population_N , z = ~CV2_ybar_s*100)
#fig3
#Changes are barely percentiable. 

#-------------------------------------------------------------------------------
#Section 3: Low incidence pens
#For Negative Binomial-
#-------------------------------------------------------------------------------
nb_variance <- function(mu, r){
  mu + mu^2/r
}

n_vector <- c(4:30, seq(from = 35, to = 500, by =1))
#N_vector <- c(100000) #Taken from the above section. 
mu_vector <- c(seq(from = 0.025, to = 0.5, by = 0.025),
               seq(from = 0.5+0.125, to = 2, by = 0.125))
dispersion_vector <- c(seq(from = 0.125, to = 5, by = 0.125))

nb_unconstrained_options = nrow(expand.grid(mu = mu_vector, dispersion = dispersion_vector))

nb_experiments <- expand.grid(sample_n=n_vector,
                                   population_N = N_vector,
                                   mu = mu_vector,
                                   dispersion = dispersion_vector)

nb_experiments <- nb_experiments %>% mutate(
  Su2 = nb_variance(mu, dispersion)
) %>%
mutate(CV2_ybar_s = CV2_ybar_s(sample_n, population_N, Su2, mu))

#Implement our two rules:
nb_experiments <- nb_experiments %>%
  filter(Su2 >= mu) %>% filter(Su2 <= 5)

sensible_nb_params <- nb_experiments %>%
  select(mu, dispersion) %>%
  unique()

sensible_nb_dists <- nb_experiments %>%
  select(mu, dispersion, Su2) %>%
  unique()

# Confidence interval ----------------------------------------------------------
nb_variance <- function(mu, r){
  mu + mu^2/r
}

p_from_r_and_mu <- function(r, mu){
  r/(mu+r)
}


mu_from_r_and_p <- function(r, p){
  #Note, p is the probability of successes!
  # mu is the mean of failure.
  r*(1-p)/p
}

#Colloborated this against another paper (Kismoodthy, 2005) for both bounds

#failures, successes, n_samples, confidence.
CI_CB <- function(x, r, n_samples = 1, alpha = 0.05){
  n = r * n_samples # number of successes before you stop counting failures
  x = x * n_samples #number of failures when you stop counting
  p = n/(x+n) #probability of success per trial
  alpha <- alpha/2 #confidence interval
  
  #CI bound on P
  lower_p <- n/(n + (x + 1) * qf(1 - alpha, 2 * (x + 1), 2 * n))
  upper_p <- n * qf(1 - alpha, 2 * n, 2 * x)/(n * qf(1 - alpha, 2 * n, 2 * x) + x)
  
  #Total number of trails remains the same. Whereas the new bound on p
  # means the number of successes change.
  total_trails <- x + n
  n_successes_lower <- total_trails * lower_p
  n_successes_upper <- total_trails * upper_p
  
  
  #CI bound on mu
  #a higher p gives a higher probability of success (n successes increases).
  #therefore a high p is associated with the lowesest number of failures.
  mu = mu_from_r_and_p(n, p) / n_samples
  upper_mu = mu_from_r_and_p(n_successes_lower, lower_p) / n_samples
  lower_mu = mu_from_r_and_p(n_successes_upper, upper_p) / n_samples
  
  
  #Return
  list(p = p,
       lower_p = lower_p,
       upper_p = upper_p,
       dif_p = upper_p - lower_p,
       mu = mu,
       lower_mu = lower_mu,
       upper_mu = upper_mu,
       dif_mu = upper_mu - lower_mu)
}


determine_CI_CB <- function(mu, dispersion, n_samples = 1){
  p_successes = p_from_r_and_mu(dispersion, mu)
  n_successes = dispersion #In the traditional frame: total_trials  * p_successes
  total_trials = n_successes / p_successes # Need not be integer in neg-binomial
  n_failures = total_trials * (1- p_successes)
  #Of course it follows that n_successes is total_trails * (p_successes)
  
  #We quickly validate the mean using the fact that the theoretical mean
  # from traditional parameterisation is:
  theoretical_mean <-  (n_successes * (1-p_successes) / p_successes)
  
  if(isFALSE(all.equal(theoretical_mean, mu))){
    print(theoretical_mean)
    print(mu)
    stop("Critical misalignment of frameworks.")
  }
  
  #We provide the main output.
  return(CI_CB(n_failures, n_successes, n_samples = n_samples))
}



# Generate Error Tables -------------------------------------------------------


nb_expers <- expand.grid(mu = mu_vector,
            dispersion = dispersion_vector,
            sample_n = n_vector,
            population_N = N_vector)


nb_expers <- nb_expers %>% mutate(
  sample_n = as.integer(sample_n),
  Su2 = nb_variance(mu, dispersion)
) %>%
  mutate(CV2_ybar_s = CV2_ybar_s(sample_n, population_N, Su2, mu))


#Implement two ecological rules, outlined in the paper.:
nb_expers <- nb_expers %>%
  filter(Su2 >= mu) %>% filter(Su2 <= 5)


all_CI_CB <- determine_CI_CB(nb_expers$mu, nb_expers$dispersion, nb_expers$sample_n)
nb_expers$CI_upper_bound <- all_CI_CB$upper_mu
nb_expers$CI_lower_bound <- all_CI_CB$lower_mu
nb_expers$eh_CB <- all_CI_CB$upper_mu - all_CI_CB$mu
nb_expers$el_CB <- all_CI_CB$mu - all_CI_CB$lower_mu

#openxlsx::write.xlsx(nb_expers, 'nb_experiments.xlsx')



#Numerical implementation of the same idea. -------------------------------------

#NOT RUN
# #Finally we can do the same thing numerically.
# nb_numerical_CI <- function(mu, dispersion,
#                             n_samples,
#                             quantiles = c(0.025, 0.05, 0.25, .75, .95, 0.975)){
# N = 100000 #True population.
# reps_sample = 2000 #Increasing reps reduces numerical error.
# pop <- rnbinom(N, size = dispersion, mu = mu)
# 
# #Repeated sampling.
# repeated_sampling_mean <- numeric(length = reps_sample) 
# for (r in 1:reps_sample) {
#   s <- sample(pop, n_samples)
#   repeated_sampling_mean[r] <- mean(s)
# }
# #Difference the the true mean and sample mean and report the order.
# return(quantile(mu - repeated_sampling_mean, quantiles))
# }



#NOT RUN: very time consuming. 
# percentage_bar <- as.integer(seq(from = 1, to = nrow(nb_expers), length.out = 1000))
# nb_numerical_CI_results <- foreach(i = 1:nrow(nb_expers)) %do%{
#   if(i %in% percentage_bar){
#     cat('i is', i, ':', round(i/nrow(nb_expers)*100,2), '%\n')
#   }
#   nb_numerical_CI(nb_expers[[i, "mu"]], nb_expers[[i, "dispersion"]], nb_expers[[i, "sample_n"]])
# }
# 
# save_data <- list(data = nb_expers,
#      results = nb_numerical_CI_results)
# 
# #saveRDS(save_data, 'nb_numerical_experiments.Rds')

#------------------------------------------------------------------------------
#Section 4: Combined table for comparison
#------------------------------------------------------------------------------
#Take the garwood_poisson table.
#Merge into this the CV error at a given fixed N.

p_tidy <- filter(poisson_experiments, population_N == 50000) %>%
  select(y_bar = mean_lice,
         n= sample_n, CV2_ybar_s) %>% #Took out Vybar_s
  right_join(garwood_poisson) %>%
  mutate(Su2 = y_bar) #Since Piosson.

#For the nb experiments take just one population and then tidy the output to align
#with p_tidy.
nb_tidy <- filter(nb_expers, population_N == 50000) %>%
  select(y_bar = mu,
        dispersion,
        Su2,
        n = sample_n,
        CV2_ybar_s,
        CI_lower_bound, 
        CI_upper_bound,
        el = el_CB,
        eh = eh_CB)

output_table1 <- dplyr::bind_rows(p_tidy, nb_tidy)

#openxlsx::write.xlsx(output_table1, 'output_table1.xlsx')
saveRDS(output_table1, "output_table1.Rds")

#Group the output table so that we can report on less.
output_table1 <- mutate(output_table1,
  true_mu_group = case_when(
  y_bar < 0.1 ~ "0-0.099",
  (y_bar >= 0.1) & (y_bar < 0.5) ~ "0.1-0.5",
  (y_bar >= 0.51) & (y_bar < 0.99) ~ "0.51-0.99",
  (y_bar >= 0.99) & (y_bar < 1.99) ~ "1.0-1.99",
  (y_bar >= 1.99) & (y_bar < 2.99) ~ "2.0-2.99",
  (y_bar >= 2.99) & (y_bar < 3.99) ~ "3.0-3.99",
  (y_bar >= 3.99) & (y_bar < 4.99) ~ "4.0-4.99",
  TRUE ~ "other"
 )
 )


sum_output_table2 <- output_table1 %>%
  filter(n %in% c(5, 10, 15, 30, 50, 100, 250)) %>%
  group_by(true_mu_group, n) %>%
  summarise(Su2_worst = quantile(Su2, 0.9, names = FALSE, type = 3),
            Su2_median = quantile(Su2, 0.5, names = FALSE, type = 3),
            Su2_best = quantile(Su2, 0.1, names = FALSE, type = 3)
            )

sum_output_table3 <- sum_output_table2 %>%
  pivot_longer(Su2_worst:Su2_best, values_to = "Su2")



ll3 <- sum_output_table3 %>% ungroup() %>%
  #mutate(Su2_text = as.character(round(Su2, 5))) %>%
  left_join(output_table1,  #%>% mutate(Su2_text = as.character(round(Su2, 5))),
            by = c("true_mu_group", "n", "Su2"),
            copy = TRUE)

num_tidy <- function(x, d = 3){
  format(round(x, d), n_small = d)
}

ll4 <- ll3 %>% mutate(out_text = paste0("[", num_tidy(el,3), ",", num_tidy(eh, 3), "]",
                                 "(", round(CV2_ybar_s*100, 0), "%", ")")
               ) %>%
  select(true_mu_group, n, name, out_text) %>%
  filter(true_mu_group != "other")

ll5 <- ll4 %>% pivot_wider(id_cols = c("true_mu_group", "n"), values_from = "out_text")
knitr::kable(ll5, "latex")
#Figures ---------------------------------------------------------------------


# #Great Plot!
# nb_dispersion <- nb_experiments %>% select(mu, dispersion, Su2) %>% unique()
# 
# nb_pivot <- nb_dispersion %>%
#   tidyr::pivot_wider(names_from = mu, values_from = Su2)
# nb_mu_vector <- as.numeric(colnames(nb_pivot)[-1])
# 
# nb_matrix <- as.matrix(select(nb_pivot, -dispersion))
# 
# fig6 <- plot_ly(nb_experiments, x = ~nb_mu_vector,
#                 y = ~nb_pivot$dispersion,
#                 z = ~nb_matrix,
#                 type = "surface")
# fig6

#figures/CV2bar_ys
plot3 <- ggplot(data = nb_experiments, aes(sample_n, CV2_ybar_s, colour = Su2)) + 
  geom_point() + 
  coord_cartesian(ylim = c(0, 1.50)) +
  geom_hline(yintercept = 0.1) + 
  coord_cartesian(xlim = c(0, 30)) + 
  scale_y_continuous(labels = scales::percent) +
  theme_bw() + 
  theme(text=element_text(size=16,  family="serif")) +
  ylab(TeX("${CV^{2}}(\\bar{y}_s)$")) +
  xlab("Sample (n)") +
  labs(color = TeX("$S_{u}^{2}$"))

ggsave("plot3.png", plot = plot3,
       width = 16, height = 10, units = "cm", dpi = 650)

#figures/nb-mu-dispersion-su2.pdf
plot4 <- ggplot(data = nb_experiments, aes(dispersion, Su2, colour = mu)) +
  geom_line(aes(group = mu)) + 
  theme_bw() + 
  theme(text=element_text(size=16,  family="serif")) +
  ylab(TeX("$S_{u}^{2}$")) +
  xlab("r") +
  labs(color = TeX("$\\mu$"))

ggsave("plot4.png", plot = plot4,
       width = 9, height = 30, units = "cm", dpi = 650)