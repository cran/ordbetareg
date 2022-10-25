## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center",
  fig.width=5,
  fig.height=3
)
options(rmarkdown.html_vignette.check_title = FALSE)

## ----setup, include = FALSE---------------------------------------------------

library(ordbetareg)
library(dplyr)
library(ggplot2)
library(haven)
library(brms)
library(tidyr)
library(stringr)
library(Hmisc)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----runmod-------------------------------------------------------------------

# whether to run models from scratch
run_model <- F


## ----load_data----------------------------------------------------------------

data("pew")

pew %>% 
  ggplot(aes(x=as.numeric(therm))) +
  geom_histogram(bins=100) +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  scale_x_continuous(breaks=c(0,25,50,75,100),
                     labels=c("0","Colder","50","Warmer","100")) +
  ylab("") +
  xlab("") +
  labs(caption=paste0("Figure shows the distribution of ",sum(!is.na(pew$therm))," non-missing survey responses."))


## ----munge_data---------------------------------------------------------------

model_data <- select(pew,therm,age="F_AGECAT_FINAL",
                        sex="F_SEX_FINAL",
                        income="F_INCOME_FINAL",
                        ideology="F_IDEO_FINAL",
                        race="F_RACETHN_RECRUITMENT",
                        education="F_EDUCCAT2_FINAL",
                     region="F_CREGION_FINAL",
                        approval="POL1DT_W28",
                       born_again="F_BORN_FINAL",
                       relig="F_RELIG_FINAL",
                        news="NEWS_PLATFORMA_W28") %>% 
    mutate_at(c("race","ideology","income","approval","sex","education","born_again","relig"), function(c) {
      factor(c, exclude=levels(c)[length(levels(c))])
    }) %>% 
    # need to make these ordered factors for BRMS
    mutate(education=ordered(education),
           income=ordered(income))


## ----run_ordbetareg-----------------------------------------------------------

if(run_model) {
  
  ord_fit_mean <- ordbetareg(formula=therm ~ mo(education)*mo(income) +
                               (1|region), 
                       data=model_data,
                cores=2,chains=2,iter=1000,
                refresh=0)
                # NOTE: to do parallel processing within chains
                # add the options below
                #threads=threading(5),
                #backend="cmdstanr"
                #where threads is the number of cores per chain
                # you must have cmdstanr set up to do so
                # see https://mc-stan.org/cmdstanr/
  
} else {
  
  data("ord_fit_mean")
  
}

## ----plot_cut-----------------------------------------------------------------

all_draws <- prepare_predictions(ord_fit_mean)

cutzero <- plogis(all_draws$dpars$cutzero)
cutone <- plogis(all_draws$dpars$cutzero + exp(all_draws$dpars$cutone))

pew %>% 
  ggplot(aes(x=therm)) +
  geom_histogram(bins=100) +
  theme_minimal() + 
  theme(panel.grid=element_blank()) +
  scale_x_continuous(breaks=c(0,25,50,75,100),
                     labels=c("0","Colder","50","Warmer","100")) +
  geom_vline(xintercept = mean(cutzero)*100,linetype=2) +
  geom_vline(xintercept = mean(cutone)*100,linetype=2) +
  ylab("") +
  xlab("") +
  labs(caption=paste0("Figure shows the distribution of ",sum(!is.na(pew$therm))," non-missing survey responses."))



## ----post_predict,message=F---------------------------------------------------

plots <- pp_check_ordbeta(ord_fit_mean,ndraws=100)

plots$discrete
plots$continuous


## ----coef_plot----------------------------------------------------------------

library(modelsummary)

modelsummary(ord_fit_mean,statistic = "conf.int",
                          metrics = "RMSE",
                          coef_map=c("b_Intercept"="Intercept",
                                     "bsp_moeducation"="Education",
                                     "bsp_moincome"="Income",
                                     "bsp_moeducation:moincome"="EducationXIncome"))

## ----marg_effect--------------------------------------------------------------

library(marginaleffects)

marg_effs <- marginaleffects(ord_fit_mean)

summary(marg_effs) %>% knitr::kable()


## ----run_brms_phi-------------------------------------------------------------

if(run_model) {
  
  ord_fit_phi <- ordbetareg(bf(therm ~ 0, 
                phi ~ age + sex),
                phi_reg = T,
                data=model_data,
                cores=2,chains=2,iter=1000,
                refresh=0)
                # NOTE: to do parallel processing within chains
                # add the options below
                #threads=threading(5),
                #backend="cmdstanr"
                #where threads is the number of cores per chain
                # you must have cmdstanr set up to do so
                # see https://mc-stan.org/cmdstanr/
  
} else {
  
  data("ord_fit_phi")
  
}



## ----phicoef------------------------------------------------------------------

summary(ord_fit_phi)


## ----plot_phi_sim-------------------------------------------------------------

# we can use some dplyr functions to make this really easy

female_data <- distinct(model_data, age) %>% 
  mutate(sex="Female")

male_data <- distinct(model_data, age) %>% 
  mutate(sex="Male")

to_predict <- bind_rows(female_data,
                        male_data) %>% 
  filter(!is.na(age))

pred_post <- posterior_predict(ord_fit_phi,
                               newdata=to_predict)

# better with iterations as rows

pred_post <- t(pred_post)
colnames(pred_post) <- 1:1000

# need to convert to a data frame

data_pred <- as_tibble(pred_post) %>% 
  mutate(sex=to_predict$sex,
         age=to_predict$age) %>% 
  gather(key="iter",value='estimate',-sex,-age)

data_pred %>% 
  ggplot(aes(x=estimate)) +
  geom_density(alpha=0.5,aes(fill=sex)) +
  theme(panel.background = element_blank(),
        panel.grid=element_blank())


## ----check_data,eval=FALSE----------------------------------------------------
#  
#  # NOT RUN IN THE VIGNETTE
#  
#    single_data <- sim_ordbeta(N=100,iter=1,
#                               return_data=T)
#  
#  # examine the first dataset
#  
#  knitr::kable(head(single_data$data[[1]]))
#  

## ----sim_data_full------------------------------------------------------------

if(run_model) {
  
  sim_data <- sim_ordbeta(N=c(250,500,750),
                          k=1,
                          beta_coef = .5,
                          iter=100,cores=10,
                          beta_type="binary",
                          treat_assign=0.3)
  
} else {
  
  data("sim_data")
  
}



## ----sim_plot-----------------------------------------------------------------

sim_data %>% 
    select(`Proportion S Errors`="s_err",N,Power="power",
         `M Errors`="m_err",Variance="var_marg") %>% 
  gather(key = "type",value="estimate",-N) %>%
  ggplot(aes(y=estimate,x=N)) +
  #geom_point(aes(colour=model),alpha=0.1) +
  stat_summary(fun.data="mean_cl_boot") + 
  ylab("") +
  xlab("N") +
  scale_x_continuous(breaks=c(250,500,750)) +
  scale_color_viridis_d() +
  facet_wrap(~type,scales="free_y",ncol = 2) +
  labs(caption=stringr::str_wrap("Summary statistics calculated as mean with bootstrapped confidence interval from simulation draws. M Errors  and S errors are magnitude of bias (+1 equals no bias) and incorrect sign of the estimated marginal effect respectively. Variance refers to estimated posterior variance (uncertainty) of the marginal effect(s).",width=50)) +
  guides(color=guide_legend(title=""),
         linetype=guide_legend(title="")) +
  theme_minimal() +
  theme(plot.caption = element_text(size=7),
        axis.text.x=element_text(size=8))


