
# Fit the random effects model with REMA
# (1) compare estimating biomass using the GOA-wide vs. stratified approach
# (i.e., split by EGOA, CGOA, and WGOA)
# (2) explore also fitting to the IPHC survey relative population numbers


# set up ----

library(tidyverse)
library(janitor)
library(rema) # pak::pkg_install('afsc-assessments/rema')
library(cowplot)

ggplot2::theme_set(cowplot::theme_cowplot(font_size = 12) +
                     cowplot::background_grid() +
                     cowplot::panel_border())

# data clean up ----

# current trawl survey biomass data
biomass_dat <- read_csv(here::here('data', 'raw_csv', 'all_skate_trawl_biomass.csv')) %>% 
  janitor::clean_names(case = 'snake') %>% # put column names into snake case format (lower case with underscores instead of spaces)
  mutate(species_group = dplyr::case_when(common_name == 'big skate' ~ 'big_skate', # case_when is a more concise alternative to ifelse
                                          common_name == 'longnose skate' ~ 'longnose_skate',
                                          .default = 'other_skates'),
         strata = dplyr::case_when(nmfs_reporting_area %in% c('Shumagin') ~ 'WGOA',
                                   nmfs_reporting_area %in% c('Chirikof', 'Kodiak') ~ 'CGOA',
                                   nmfs_reporting_area %in% c('Yakutat', 'Southeastern') ~ 'EGOA')) 

# biomass by EGOA, CGOA, WGOA
biomass_by_strata <- biomass_dat %>% 
  group_by(species_group, year, strata) %>% 
  # format needed for rema
  summarise(biomass = sum(area_biomass, na.rm = TRUE),
            cv = sqrt(sum(area_biomass_var, na.rm = TRUE)) / biomass) %>% 
  ungroup() %>% 
  write_csv(here::here('data', 'rema_inputs', 'goa_skates_strata_biomass.csv'))

# total GOA biomass
total_biomass <- biomass_dat %>% 
  mutate(strata = 'GOA') %>% 
  group_by(species_group, strata, year) %>% 
  # format needed for rema
  summarise(biomass = sum(area_biomass, na.rm = TRUE),
            cv = sqrt(sum(area_biomass_var, na.rm = TRUE)) / biomass) %>% 
  ungroup() %>% 
  write_csv(here::here('data', 'rema_inputs', 'goa_skates_total_biomass.csv'))

# the iphc longline survey relative population numbers, where rpns are
# area-weighted numbers of fish per hook. this *can* be used as a cpue index in
# the model that inform biomass TREND (where trawl survey biomass informs
# SCALE). the current model does NOT use rpns, but i'm showing you a potential
# alternative model to explore for future years if you're interested.
cpue_dat <- read_csv(here::here('data', 'raw_csv', 'all_skate_iphc_rpns.csv')) %>% 
  janitor::clean_names(case = 'snake') %>% 
  mutate(species_group = dplyr::case_when(species == 'Big skate' ~ 'big_skate',
                                          species == 'Longnose skate' ~ 'longnose_skate',
                                          .default = 'other_skates'),
         strata = dplyr::case_when(fmp_sub_area %in% c('EY/SE', 'WY') ~ 'EGOA',
                                   .default = fmp_sub_area)) 

cpue_by_strata <- cpue_dat %>% 
  group_by(species_group, year, strata) %>% 
  summarise(cpue = sum(boot_strata_rpn, na.rm = TRUE),
            cv = sqrt(sum(boot_sd^2, na.rm = TRUE)) / cpue) %>% 
  ungroup() %>% 
  write_csv(here::here('data', 'rema_inputs', 'goa_skates_strata_cpue.csv'))

total_cpue <- cpue_dat %>% 
  mutate(strata = 'GOA') %>% 
  group_by(species_group, strata, year) %>% 
  summarise(cpue = sum(boot_strata_rpn, na.rm = TRUE),
            cv = sqrt(sum(boot_sd^2, na.rm = TRUE)) / cpue) %>% 
  ungroup() %>% 
  write_csv(here::here('data', 'rema_inputs', 'goa_skates_total_cpue.csv'))

# current model (goa wide) ----

# currently the univariate RE.tpl model is fit multiple times for each species group... one as GOA-wide
# for each species groups and then again for each region/strata (EGOA, CGOA,
# WGOA). here, we just run it GOA wide

# bs = 'big_skate'
bs_input <- prepare_rema_input(biomass_dat = total_biomass %>% 
                                 filter(species_group == 'big_skate') %>% 
                                 select(-species_group),
                               end_year = 2023, # you can define the last year for model predictions, defaults to last year of data
                               model_name = 'big_skate_goa_wide')
bs_mod <- fit_rema(bs_input)
bs_out <- tidy_rema(bs_mod)
bs_plots <- plot_rema(bs_out, biomass_ylab = 'Biomass (t)')
(bs_p1 <- bs_plots$biomass_by_strata + ggtitle('Big skates'))

# ln = 'longnose_skate'
ln_input <- prepare_rema_input(biomass_dat = total_biomass %>% 
                                 filter(species_group == 'longnose_skate') %>% 
                                 select(-species_group),
                               end_year = 2023, # you can define the last year for model predictions, defaults to last year of data
                               model_name = 'longnose_skate_goa_wide')
ln_mod <- fit_rema(ln_input)
ln_out <- tidy_rema(ln_mod)
ln_plots <- plot_rema(ln_out, biomass_ylab = 'Biomass (t)')
(ln_p1 <- ln_plots$biomass_by_strata + ggtitle('Longnose skates'))

# os = 'other_skates'
os_input <- prepare_rema_input(biomass_dat = total_biomass %>% 
                                 filter(species_group == 'other_skates') %>% 
                                 select(-species_group),
                               end_year = 2023, # you can define the last year for model predictions, defaults to last year of data
                               model_name = 'other_skates_goa_wide')
os_mod <- fit_rema(os_input)
os_out <- tidy_rema(os_mod)
os_plots <- plot_rema(os_out, biomass_ylab = 'Biomass (t)')
(os_p1 <- os_plots$biomass_by_strata + ggtitle('Other skates'))

cowplot::plot_grid(bs_p1, ln_p1, os_p1, ncol = 1)
ggsave(here::here('results', 'base_goa_wide', 'fits_goa_wide.png'), units = 'in',
       height = 6, width = 7, dpi = 300, bg = 'white')

# save parameter estimates (one process error estimated for each species group)
(pe_goa_wide <- bind_rows(bs_out$parameter_estimates,
          ln_out$parameter_estimates,
          os_out$parameter_estimates) %>% 
  write_csv(here::here('results', 'base_goa_wide', 'process_error_goa_wide.png')))

# current model (by strata) ----

# currently the univariate RE.tpl model is fit multiple times... one as GOA-wide
# for each species groups and then again for each region/strata (EGOA, CGOA,
# WGOA). here, we just run it once for each species group

# bs = 'big_skate'
strata_bs_input <- prepare_rema_input(biomass_dat = biomass_by_strata %>% 
                                 filter(species_group == 'big_skate') %>% 
                                 select(-species_group),
                               end_year = 2023, # you can define the last year for model predictions, defaults to last year of data
                               model_name = 'big_skate_by_strata')
strata_bs_mod <- fit_rema(strata_bs_input)
strata_bs_out <- tidy_rema(strata_bs_mod)
strata_bs_plots <- plot_rema(strata_bs_out, biomass_ylab = 'Biomass (t)')
(strata_bs_p1 <- strata_bs_plots$biomass_by_strata + ggtitle('Big skates'))
(strata_bs_p2 <- strata_bs_plots$total_predicted_biomass + expand_limits(y = 0) + ggtitle('Big skates'))

# ln = 'longnose_skate'
strata_ln_input <- prepare_rema_input(biomass_dat = biomass_by_strata %>% 
                                 filter(species_group == 'longnose_skate') %>% 
                                 select(-species_group),
                               end_year = 2023, # you can define the last year for model predictions, defaults to last year of data
                               model_name = 'longnose_skate_by_strata')
strata_ln_mod <- fit_rema(strata_ln_input)
strata_ln_out <- tidy_rema(strata_ln_mod)
strata_ln_plots <- plot_rema(strata_ln_out, biomass_ylab = 'Biomass (t)')
(strata_ln_p1 <- strata_ln_plots$biomass_by_strata + ggtitle('Longnose skates'))
(strata_ln_p2 <- strata_ln_plots$total_predicted_biomass + expand_limits(y = 0) + ggtitle('Longnose skates'))

# os = 'other_skates'
strata_os_input <- prepare_rema_input(biomass_dat = biomass_by_strata %>% 
                                 filter(species_group == 'other_skates') %>% 
                                 select(-species_group),
                               end_year = 2023, # you can define the last year for model predictions, defaults to last year of data
                               model_name = 'other_skates_by_strata')
strata_os_mod <- fit_rema(strata_os_input)
strata_os_out <- tidy_rema(strata_os_mod)
strata_os_plots <- plot_rema(strata_os_out, biomass_ylab = 'Biomass (t)')
(strata_os_p1 <- strata_os_plots$biomass_by_strata + ggtitle('Other skates'))
(strata_os_p2 <- strata_os_plots$total_predicted_biomass + expand_limits(y = 0) + ggtitle('Other skates'))

cowplot::plot_grid(strata_bs_p1, strata_ln_p1, strata_os_p1, ncol = 1)
ggsave(here::here('results', 'base_by_strata', 'fits_by_strata.png'), units = 'in',
       height = 9, width = 9, dpi = 300, bg = 'white')

# save parameter estimates (one process error estimated for each strata and each species group)
(pe_by_strata <- bind_rows(strata_bs_out$parameter_estimates,
                           strata_ln_out$parameter_estimates,
                           strata_os_out$parameter_estimates) %>% 
    mutate(strata = rep(unique(biomass_by_strata$strata), 3)) %>% 
    select(model_name, strata, parameter, estimate, std_err, lci, uci) %>% 
    write_csv(here::here('results', 'base_by_strata', 'process_error_by_strata.csv')))

# compare base methods ----

# compare total biomass estimates when using the "goa wide" vs "by strata" approaches
bind_rows(bind_rows(bs_out$biomass_by_strata %>% mutate(model_name = 'biomass_by_strata'),
                    strata_bs_out$total_predicted_biomass %>% mutate(model_name = 'biomass_goa_wide')) %>% 
            mutate(species_group = 'Big skates'), 
          bind_rows(ln_out$biomass_by_strata %>% mutate(model_name = 'biomass_by_strata'),
                    strata_ln_out$total_predicted_biomass %>% mutate(model_name = 'biomass_goa_wide')) %>% 
            mutate(species_group = 'Longnose skates'),
          bind_rows(os_out$biomass_by_strata %>% mutate(model_name = 'biomass_by_strata'),
                    strata_os_out$total_predicted_biomass %>% mutate(model_name = 'biomass_goa_wide')) %>% 
            mutate(species_group = 'Other skates')) %>% 
ggplot(aes(x = year, y = pred,
           col = model_name)) +
  geom_ribbon(aes(ymin = pred_lci, ymax = pred_uci,
                  fill = model_name), col = NA,
              alpha = 0.25) +
  geom_line() +
  facet_wrap(~species_group, ncol = 1) +
  geom_point(aes(x = year, y = obs), col = 'black') +
  # did not include error bars because the underlying CVs are different under
  # the two methods and this way of plotting only shows the GOA wide
  # geom_errorbar(aes(x = year, ymin = obs_lci, ymax = obs_uci), col = 'black') +
  scale_y_continuous(labels = scales::comma, expand = c(0.01, 0), limits = c(0, NA)) +
  labs(x = NULL, y = 'Biomass (t)', 
       title = 'Compare total biomass predictions when estimated GOA-wide vs. by strata',
       fill = NULL, colour = NULL, shape = NULL) +
  ggplot2::scale_fill_viridis_d(direction = 1) +
  ggplot2::scale_colour_viridis_d(direction = 1)

ggsave(here::here('results', 'compare_by_strata_vs_goa_wide.png'), units = 'in',
       height = 9, width = 9, dpi = 300, bg = 'white')

# note that the fits are very similar but are in fact different (the underlying
# data going in is summarized differently and the stratified model has three
# process error parameters whereas the goa-wide model only has one per species
# group). i recommend just picking one instead of presenting both because it's
# otherwise confusing to the reader....

# explore using iphc rpn index ----

# here i'm only showing the stratified approach with the second survey but you could do something very
# similar using the goa-wide method

# bs = 'big_skate'
twosrv_bs_input <- prepare_rema_input(multi_survey = 1, # use two indices of abundance instead of one
                                      biomass_dat = biomass_by_strata %>% 
                                        filter(species_group == 'big_skate') %>% 
                                        select(-species_group),
                                      cpue_dat = cpue_by_strata %>% 
                                        filter(species_group == 'big_skate') %>% 
                                        select(-species_group),
                                      sum_cpue_index = TRUE, # is the CPUE index summable? yes, RPNs are summable because they're area-weighted
                                      end_year = 2023, # you can define the last year for model predictions, defaults to last year of data
                                      model_name = 'big_skate_two_surveys')
twosrv_bs_mod <- fit_rema(twosrv_bs_input)
twosrv_bs_out <- tidy_rema(twosrv_bs_mod)
twosrv_bs_plots <- plot_rema(twosrv_bs_out, biomass_ylab = 'Biomass (t)', cpue_ylab = 'Relative population numbers')
twosrv_bs_p1 <- twosrv_bs_plots$biomass_by_strata + ggtitle('Big skates bottom trawl survey biomass')
twosrv_bs_p2 <- twosrv_bs_plots$cpue_by_strata + ggtitle('Big skates IPHC relative populations numbers')
cowplot::plot_grid(twosrv_bs_p1, twosrv_bs_p2, ncol = 1)
ggsave(here::here('results', 'two_surveys', 'big_skate_fits_two_surveys.png'), units = 'in',
       height = 9, width = 9, dpi = 300, bg = 'white')

# ln = 'longnose_skate'
twosrv_ln_input <- prepare_rema_input(multi_survey = 1, # use two indices of abundance instead of one
                                      biomass_dat = biomass_by_strata %>% 
                                        filter(species_group == 'longnose_skate') %>% 
                                        select(-species_group),
                                      cpue_dat = cpue_by_strata %>% 
                                        filter(species_group == 'longnose_skate') %>% 
                                        select(-species_group),
                                      sum_cpue_index = TRUE, # is the CPUE index summable? yes, RPNs are summable because they're area-weighted
                                      end_year = 2023, # you can define the last year for model predictions, defaults to last year of data
                                      model_name = 'longnose_skate_two_surveys')
twosrv_ln_mod <- fit_rema(twosrv_ln_input)
twosrv_ln_out <- tidy_rema(twosrv_ln_mod)
twosrv_ln_plots <- plot_rema(twosrv_ln_out, biomass_ylab = 'Biomass (t)', cpue_ylab = 'Relative population numbers')
twosrv_ln_p1 <- twosrv_ln_plots$biomass_by_strata + ggtitle('Longnose skates bottom trawl survey biomass')
twosrv_ln_p2 <- twosrv_ln_plots$cpue_by_strata + ggtitle('Longnose skates IPHC relative populations numbers')
cowplot::plot_grid(twosrv_ln_p1, twosrv_ln_p2, ncol = 1)
ggsave(here::here('results', 'two_surveys', 'longnose_skate_fits_two_surveys.png'), units = 'in',
       height = 9, width = 9, dpi = 300, bg = 'white')

# os = 'other_skates'
twosrv_os_input <- prepare_rema_input(multi_survey = 1, # use two indices of abundance instead of one
                                      biomass_dat = biomass_by_strata %>% 
                                        filter(species_group == 'other_skates') %>% 
                                        select(-species_group),
                                      cpue_dat = cpue_by_strata %>% 
                                        filter(species_group == 'other_skates') %>% 
                                        select(-species_group),
                                      sum_cpue_index = TRUE, # is the CPUE index summable? yes, RPNs are summable because they're area-weighted
                                      end_year = 2023, # you can define the last year for model predictions, defaults to last year of data
                                      model_name = 'other_skates_two_surveys')
twosrv_os_mod <- fit_rema(twosrv_os_input)
twosrv_os_out <- tidy_rema(twosrv_os_mod)
twosrv_os_plots <- plot_rema(twosrv_os_out, biomass_ylab = 'Biomass (t)', cpue_ylab = 'Relative population numbers')
twosrv_os_p1 <- twosrv_os_plots$biomass_by_strata + ggtitle('Other skates bottom trawl survey biomass')
twosrv_os_p2 <- twosrv_os_plots$cpue_by_strata + ggtitle('Other skates IPHC relative populations numbers')
cowplot::plot_grid(twosrv_os_p1, twosrv_os_p2, ncol = 1)
ggsave(here::here('results', 'two_surveys', 'other_skates_fits_two_surveys.png'), units = 'in',
       height = 9, width = 9, dpi = 300, bg = 'white')

# compare all methods ----

p1 <- compare_rema_models(rema_models = list(bs_mod, strata_bs_mod, twosrv_bs_mod), biomass_ylab = 'Biomass (t)')$plots$total_predicted_biomass +
  ggtitle('Big skates')
p2 <- compare_rema_models(rema_models = list(ln_mod, strata_ln_mod, twosrv_ln_mod), biomass_ylab = 'Biomass (t)')$plots$total_predicted_biomass +
  ggtitle('Longnose skates')
p3 <- compare_rema_models(rema_models = list(os_mod, strata_os_mod, twosrv_os_mod), biomass_ylab = 'Biomass (t)')$plots$total_predicted_biomass +
  ggtitle('Other skates')
cowplot::plot_grid(p1, p2, p3, ncol = 1)
ggsave(here::here('results', 'compare_by_strata_vs_goa_wide_vs_two_surveys.png'), units = 'in',
       height = 9, width = 9, dpi = 300, bg = 'white')
