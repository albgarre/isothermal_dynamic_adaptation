
library(tidyverse)
library(readxl)
library(FME)
library(ggthemes)
library(wesanderson)

source("helpers.R")

## Load the data

data_wild_iso <- read_excel("./data/data_final/data_nopreadapt.xlsx", sheet = "for_R_wild") %>%
    mutate(preadapted = FALSE)

data_wild_iso_adapt <- read_excel("./data/data_final/data_preadapt.xlsx", sheet = "for_R_wild") %>%
    mutate(preadapted = TRUE)

data_mutant_iso <- read_excel("./data/data_final/data_nopreadapt.xlsx", sheet = "for_R_mutant") %>%
    mutate(preadapted = FALSE)

data_mutant_iso_adapt <- read_excel("./data/data_final/data_preadapt.xlsx", sheet = "for_R_mutant") %>%
    mutate(preadapted = TRUE)

## Analysis at 51

models_51 <- list(
    `wild-normal` = modFit(get_residuals,
                           c(D = 10, logN0 = 7),
                           data = filter(data_wild_iso, temperature == 51),
                           model = one_line
                    ),
    `wild-adapted` = modFit(get_residuals,
                            c(D1 = 10, D2 = 60, logN0 = 7, t_crit = 30),
                            lower = c(0, 0, 0, 0),
                            data = filter(data_wild_iso_adapt, temperature == 51)
                            ),
    `mutant-normal` = modFit(get_residuals,
                             c(D = 10, logN0 = 7),
                             data = filter(data_mutant_iso, temperature == 51),
                             model = one_line
                             ),
    `mutant-adapted` = modFit(get_residuals,
                           c(D1 = 10, D2 = 60, logN0 = 7, t_crit = 30),
                           lower = c(0, 0, 0, 0),
                           data = filter(data_mutant_iso_adapt, temperature == 51)
                           )
    
)

## Analysis at 52.5

models_52p5 <- list(
    `wild-normal` = modFit(get_residuals,
                           c(D1 = 1, D2 = 10, logN0 = 7, t_crit = 5),
                           data = filter(data_wild_iso, temperature == 52.5)
    ),
    `wild-adapted` = modFit(get_residuals,
                            c(D = 10, logN0 = 7),
                            data = filter(data_wild_iso_adapt, temperature == 52.5),
                            model = one_line
    ),
    `mutant-normal` = modFit(get_residuals,
                             c(D1 = 1, D2 = 10, logN0 = 7, t_crit = 5),
                             data = filter(data_mutant_iso, temperature == 52.5)
    ) ,
    `mutant-adapted` = modFit(get_residuals,
                              c(D = 1, logN0 = 7),
                              data = filter(data_mutant_iso_adapt, temperature == 52.5),
                              model = one_line
    )
    
)

## Analysis at 55

models_55 <- list(
    `wild-normal` = modFit(get_residuals,
                           c(D1 = .1, D2 = 1, logN0 = 7, t_crit = .1),
                           data = filter(data_wild_iso, temperature == 55)
    ) ,
    `wild-adapted` = modFit(get_residuals,
                            c(D = 6, logN0 = 7),
                            lower = c(0, 0),
                            data = filter(data_wild_iso_adapt, temperature == 55),
                            model = one_line
    ),
    `mutant-normal` = modFit(get_residuals,
                             c(D1 = .3, D2 = 1, logN0 = 7, t_crit = .1),
                             lower = c(0, 0, 5, 0),
                             upper = c(10, 10, 9, 2),
                             data = filter(data_mutant_iso, temperature == 55)
    ) ,
    `mutant-adapted` = modFit(get_residuals,
                              c(D = .1, logN0 = 7),
                              data = filter(data_mutant_iso_adapt, temperature == 55),
                              model = one_line
    )
    
)

## Fitted models - Figure 1

p1 <- bind_rows(data_wild_iso,
          data_wild_iso_adapt) %>%
    filter(temperature == 51) %>%
    make_plot(., list(models_51$`wild-normal`,
                      models_51$`wild-adapted`)) +
    theme(panel.grid.major = element_line(size = 0))

p2 <- bind_rows(data_wild_iso,
          data_wild_iso_adapt) %>%
    filter(temperature == 52.5) %>%
    make_plot(., list(models_52p5$`wild-normal`,
                      models_52p5$`wild-adapted`)) +
    theme(panel.grid.major = element_line(size = 0))

p3 <- bind_rows(data_wild_iso,
          data_wild_iso_adapt) %>%
    filter(temperature == 55) %>%
    make_plot(., list(models_55$`wild-normal`,
                      models_55$`wild-adapted`)) +
    theme(panel.grid.major = element_line(size = 0))

p <- cowplot::plot_grid(p1, p2, p3, labels = "AUTO",
                   nrow = 1)

ggsave(p,
       filename = "Figure_1.png",
       width = 15,
       height  = 5)

## Fitted models (mutant) - Sup. Figure 1


p1 <- bind_rows(data_wild_iso,
                data_mutant_iso) %>%
    filter(temperature == 51) %>%
    make_plot2(., list(models_51$`wild-normal`,
                       models_51$`mutant-normal`)) +
    theme(panel.grid.major = element_line(size = 0))


p2 <- bind_rows(data_wild_iso,
                data_mutant_iso) %>%
    filter(temperature == 52.5) %>%
    make_plot2(., list(models_52p5$`wild-normal`,
                       models_52p5$`mutant-normal`)) +
    theme(panel.grid.major = element_line(size = 0))

p3 <- bind_rows(data_wild_iso,
                data_mutant_iso) %>%
    filter(temperature == 55) %>%
    make_plot2(., list(models_55$`wild-normal`,
                       models_55$`mutant-normal`)) +
    theme(panel.grid.major = element_line(size = 0))

p <- cowplot::plot_grid(p1, p2, p3, labels = "AUTO",
                   nrow = 1)

ggsave(p,
       filename = "supp_Figure_1.png",
       width = 15,
       height  = 5)

## Fitted models (mutant adapted) - Figure 2

p1 <- bind_rows(data_mutant_iso,
                data_mutant_iso_adapt) %>%
    filter(temperature == 51) %>%
    make_plot3(., list(models_51$`mutant-normal`,
                       models_51$`mutant-adapted`)) +
    theme(panel.grid.major = element_line(size = 0))


p2 <- bind_rows(data_mutant_iso,
                data_mutant_iso_adapt) %>%
    filter(temperature == 52.5) %>%
    make_plot3(., list(models_52p5$`mutant-normal`,
                       models_52p5$`mutant-adapted`)) +
    theme(panel.grid.major = element_line(size = 0))


p3 <- bind_rows(data_mutant_iso,
                data_mutant_iso_adapt) %>%
    filter(temperature == 55) %>%
    make_plot3(., list(models_55$`mutant-normal`,
                       models_55$`mutant-adapted`)) +
    theme(panel.grid.major = element_line(size = 0))


p <- cowplot::plot_grid(p1, p2, p3, labels = "AUTO",
                   nrow = 1)

ggsave(p,
       filename = "Figure_2.png",
       width = 15,
       height  = 5)

## Figure 3

my_cols <- wes_palette("Moonrise2", 2)

p1 <- models_51 %>%
    map(~ summary(.)$par) %>%
    map(as.data.frame) %>%
    map(., ~ rownames_to_column(., "par")) %>%
    imap_dfr(., ~ mutate(.x, condition = .y)) %>%
    separate(condition, into = c("strain", "state"), remove = FALSE) %>%
    filter(grepl("D", par)) %>%
    filter(par != "D1") %>%
    mutate(strain = ifelse(strain == "wild", "Wild type", "sigB mutant")) %>%
    mutate(strain = factor(strain, levels = c("Wild type", "sigB mutant"))) %>%
    ggplot(aes(x = strain, y = Estimate, colour = state)) +
    geom_point(position = position_dodge(width = .25)) +
    geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                      ymax = Estimate + `Std. Error`),
                  width = .5,
                  position = position_dodge(width = .25)) +
    theme_clean(base_size = 16) +
    theme(legend.position = "none") +
    scale_color_manual(values = rev(my_cols)) +
    xlab("") + 
    ylab("Estimated D-value (min)") +
    theme(panel.grid.major = element_line(size = 0))

p2 <- models_52p5 %>%
    map(~ summary(.)$par) %>%
    map(as.data.frame) %>%
    map(., ~ rownames_to_column(., "par")) %>%
    imap_dfr(., ~ mutate(.x, condition = .y)) %>%
    separate(condition, into = c("strain", "state"), remove = FALSE) %>%
    filter(grepl("D", par)) %>%
    mutate(strain = ifelse(strain == "wild", "Wild type", "sigB mutant")) %>%
    mutate(strain = factor(strain, levels = c("Wild type", "sigB mutant"))) %>%
    ggplot(aes(x = strain, y = Estimate, colour = state, linetype = par)) +
    geom_point(position = position_dodge(width = .25)) +
    geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                      ymax = Estimate + `Std. Error`),
                  width = .5,
                  position = position_dodge(width = .25)) +
    theme_clean(base_size = 16) +
    theme(legend.position = "none") +
    scale_color_manual(values = rev(my_cols)) +
    xlab("") + 
    ylab("Estimated D-value (min)") +
    scale_linetype_manual(values = c(1, 2, 1)) +
    theme(panel.grid.major = element_line(size = 0))

p3 <- models_55 %>%
    map(~ summary(.)$par) %>%
    map(as.data.frame) %>%
    map(., ~ rownames_to_column(., "par")) %>%
    imap_dfr(., ~ mutate(.x, condition = .y)) %>%
    separate(condition, into = c("strain", "state"), remove = FALSE) %>%
    filter(grepl("D", par)) %>%
    mutate(strain = ifelse(strain == "wild", "Wild type", "sigB mutant")) %>%
    mutate(strain = factor(strain, levels = c("Wild type", "sigB mutant"))) %>%
    ggplot(aes(x = strain, y = Estimate, colour = state, linetype = par)) +
    geom_point(position = position_dodge(width = .25)) +
    geom_errorbar(aes(ymin = Estimate - `Std. Error`,
                      ymax = Estimate + `Std. Error`),
                  width = .5,
                  position = position_dodge(width = .25)) +
    theme_clean(base_size = 16) +
    theme(legend.position = "none") +
    scale_color_manual(values = rev(my_cols)) +
    xlab("") + 
    ylab("Estimated D-value (min)") +
    scale_linetype_manual(values = c(1, 2, 1)) +
    theme(panel.grid.major = element_line(size = 0))


# cowplot::plot_grid(p1, p2, p3, labels = "AUTO", nrow = 1)
p <- cowplot::plot_grid(p2, p3, labels = "AUTO", nrow = 1)

ggsave(p,
       filename = "Figure_3.png",
       width = 10,
       height  = 5)


## Secondary models for the wild type

get_pars <- function(x) {
    x %>%
        map(~ summary(.)$par) %>%
        map(as.data.frame) %>%
        map(., ~ rownames_to_column(., "par"))  %>%
        imap_dfr(., ~ mutate(.x, condition = .y)) %>%
        separate(condition, into = c("strain", "state"), remove = FALSE) 
}


## Model fits

lower_bound_wild <- list(
    `51` = models_51,
    `52.5` = models_52p5,
    `55` = models_55
    # `57.5` = models_57p5
) %>%
    map(get_pars) %>%
    imap_dfr(., ~ mutate(.x, temperature = as.numeric(.y))) %>%
    filter(grepl("D", par)) %>%
    mutate(logD = log10(Estimate)) %>%
    filter(strain == "wild", state == "normal") %>%
    # ggplot() +
    # geom_point(aes(x = temperature, y = logD))
    group_by(temperature) %>%
    summarize(min_logD = min(logD)) %>%
    nls(min_logD ~ logDref - (temperature - 53)/z,
        data = .,
        start = list(logDref = 1, z = 5)) 

summary(lower_bound_wild)

upper_bound_wild <- list(
    `51` = models_51,
    `52.5` = models_52p5,
    `55` = models_55
    # `57.5` = models_57p5
) %>%
    map(get_pars) %>%
    imap_dfr(., ~ mutate(.x, temperature = as.numeric(.y))) %>%
    filter(grepl("D", par)) %>%
    mutate(logD = log10(Estimate)) %>%
    filter(strain == "wild") %>%
    # ggplot() + geom_point(aes(temperature, logD))
    group_by(temperature) %>%
    summarize(
              max_logD = max(logD)) %>%
    nls(max_logD ~ logDref - (temperature - 53)/z,
        data = .,
        start = list(logDref = 1, z = 5)) 

summary(upper_bound_wild)

## Figure 4

my_cols <- wes_palette("Moonrise2", 2)
my_cols1 <- wes_palette("Moonrise1", 4)[c(4,3)]

p <- list(
    `51` = models_51,
    `52.5` = models_52p5,
    `55` = models_55
) %>%
    map(get_pars) %>%
    imap_dfr(., ~ mutate(.x, temperature = as.numeric(.y))) %>%
    filter(grepl("D", par)) %>%
    mutate(logD = log10(Estimate)) %>%
    filter(strain == "wild") %>%
    ggplot(aes(x = temperature, y = logD, colour = state)) +
    geom_point() +
    geom_label(aes(label = par, fill = state)) +
    geom_line(aes(x, y),
              inherit.aes = FALSE,
              colour = my_cols[2], 
              linetype = 2,
              data = tibble(x = c(51, 55),
                            y = coef(upper_bound_wild)[1] - (x - 53)/coef(upper_bound_wild)[2])
              ) +
    geom_line(aes(x, y),
              inherit.aes = FALSE,
              colour = my_cols[1], linetype = 2,
              data = tibble(x = c(51, 55),
                            y = coef(lower_bound_wild)[1] - (x - 53)/coef(lower_bound_wild)[2])
    ) +
    theme_clean(base_size = 14) +
    xlab("Treatment temperature (ºC)") +
    ylab("Logarithm of the D-value (log min)") +
    theme(panel.grid.major = element_line(size = 0)) +
    scale_colour_manual(values = rev(my_cols1)) +
    scale_fill_manual(values = my_cols1) +
    theme(legend.position = "none")

ggsave(p,
       filename = "Figure_4.png",
       width = 6,
       height  = 4)

## Supp. Table 1


list(
    `51` = models_51,
    `52.5` = models_52p5,
    `55` = models_55
) %>%
    map(get_pars) %>%
    imap_dfr(., ~ mutate(.x, temperature = as.numeric(.y))) %>%
    write_excel_csv2(., file = "pars_Subtilis.csv")

    

