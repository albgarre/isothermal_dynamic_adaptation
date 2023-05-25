
source("analysis_final.R")

library(bioinactivation)

my_cols <- wes_palette("Moonrise2", 2)

## Load dynamic data

data_10C <- read_excel("./data/data_may/Selected data for analysis(no preadaptation)_V3.xlsx", sheet = "Wild_10C_min") %>%
    mutate(strain = "wild") %>%
    select(-Count) %>%
    mutate(Dilution = as.numeric(Dilution))

data_5C <- read_excel("./data/data_final/data_nopreadapt.xlsx", sheet = "for_R_wild_5C") %>%
    rename(time = Time, logN = `log CFU`)  %>%
    mutate(strain = "wild")

data_1C <- read_excel("./data/data_may/Selected data for analysis(no preadaptation)_V3.xlsx", sheet = "Wild_1C_min") %>%
    mutate(strain = "wild") %>%
    select(-Count) %>%
    mutate(Dilution = as.numeric(Dilution))

## Get the parameters

par_list <- list(
    lower = c(D_R = 10^coef(lower_bound_wild)[["logDref"]],
              z = coef(lower_bound_wild)[["z"]],
              temp_ref = 53,
              N0 = 1e6
    ),
    upper = c(D_R = 10^coef(upper_bound_wild)[["logDref"]],
              z = coef(upper_bound_wild)[["z"]],
              temp_ref = 53,
              N0 = 1e6
    )
)

## Comparison for 1C/min

times <- seq(0, max(data_1C$time, na.rm = TRUE), length=100) 

temp_profile <- data_1C %>%
    select(time, temperature = Temp) %>%
    as.data.frame()

d_plot <- data_1C %>%
    mutate(rep = cumsum(data_1C$time == 0)) %>%
    group_by(rep) %>%
    mutate(logN0 = mean(ifelse(time == 0, logN, NA), na.rm = TRUE)) %>%
    mutate(logS = logN - logN0)

p3 <- par_list %>%
    map(.,
        ~ predict_inactivation("Bigelow", 
                               times,
                               ., 
                               temp_profile)
    ) %>%
    map(., ~ .$simulation) %>%
    imap_dfr(., ~ mutate(.x, bound = .y)) %>%
    ggplot() +
    geom_line(aes(x = time, y = logN - 6, colour = bound),
              linetype = 2, size = 1) +
    geom_point(aes(x = time, y = logS),
               data = d_plot, size = 2, shape = 1,
               inherit.aes = FALSE) +
    coord_cartesian(ylim = c(-4, 0)) +
    scale_colour_manual(values = my_cols) +
    theme_clean(base_size = 14) +
    # cowplot::theme_cowplot() +
    xlab("Treatment time (min)") +
    ylab("Number of log-reductions (·)") +
    theme(legend.position = "none") +
    theme(panel.grid.major = element_line(size = 0))

##

min_time <- 0
max_time <- max(d_plot$time, na.rm = TRUE)

max_count <- 0
# min_count <- min(d_plot$logS)
min_count <- -4

tt <- seq(min_time, max_time, length = 100)
min_temp <- min(d_plot$Temp)
max_temp <- max(d_plot$Temp)

slope <- (max_count - min_count)/(max_temp - min_temp)
intercept <- (min_temp * (max_count-min_count)/(max_temp - min_temp) - min_count)

l3 <- d_plot %>%
    mutate(fake_temp = Temp * slope - intercept) %>%
    geom_line(aes(x = time, y = fake_temp), data = .,
              size = 1)

p3 <- p3 + l3 +
    scale_y_continuous(
        sec.axis = sec_axis(~(. + intercept)/slope,
                            name = "Temperature (ºC)"))


## Comparison for 5C/min

times <- seq(0, max(data_5C$time, na.rm = TRUE), length=100) 

temp_profile <- data_5C %>%
    select(time, temperature = Temp) %>%
    as.data.frame()

d_plot <- data_5C %>%
    mutate(., rep = cumsum(.$time == 0)) %>%
    group_by(rep) %>%
    mutate(logN0 = mean(ifelse(time == 0, logN, NA), na.rm = TRUE)) %>% 
    mutate(logS = logN - logN0)

p2 <- par_list %>%
    map(.,
        ~ predict_inactivation("Bigelow", 
                               times,
                               ., 
                               temp_profile)
    ) %>%
    map(., ~ .$simulation) %>%
    imap_dfr(., ~ mutate(.x, bound = .y)) %>%
    ggplot() +
    geom_line(aes(x = time, y = logN - 6, colour = bound),
              linetype = 2, size = 1) +
    geom_point(aes(x = time, y = logS),
               data = d_plot, size = 2, shape = 1,
               inherit.aes = FALSE) +
    coord_cartesian(ylim = c(-4, 0)) +
    scale_colour_manual(values = my_cols) +
    theme_clean(base_size = 14) +
    # cowplot::theme_cowplot() +
    xlab("Treatment time (min)") +
    ylab("Number of log-reductions (·)") +
    theme(legend.position = "none") +
    theme(panel.grid.major = element_line(size = 0))

##

min_time <- 0
max_time <- max(data_5C$time, na.rm = TRUE)

max_count <- 0
# min_count <- min(d_plot$logS)
min_count <- -4

tt <- seq(min_time, max_time, length = 100)
min_temp <- min(d_plot$Temp)
max_temp <- max(d_plot$Temp)

slope <- (max_count - min_count)/(max_temp - min_temp)
intercept <- (min_temp * (max_count-min_count)/(max_temp - min_temp) - min_count)

l2 <- d_plot %>%
    mutate(fake_temp = Temp * slope - intercept) %>%
    geom_line(aes(x = time, y = fake_temp), data = .,
              size = 1)

p2 <- p2 + l2 +
    scale_y_continuous(
        sec.axis = sec_axis(~(. + intercept)/slope,
                            name = "Temperature (ºC)"))



## Comparison for 10C/min

times <- seq(0, max(data_10C$time, na.rm = TRUE), length=100) 

temp_profile <- data_10C %>%
    select(time, temperature = Temp) %>%
    as.data.frame()

d_plot <- data_10C %>%
    mutate(., rep = cumsum(.$time == 0)) %>%
    group_by(rep) %>%
    mutate(logN0 = mean(ifelse(time == 0, logN, NA), na.rm = TRUE)) %>% 
    mutate(logS = logN - logN0) 

p1 <- par_list %>%
    map(.,
        ~ predict_inactivation("Bigelow", 
                               times,
                               ., 
                               temp_profile)
        ) %>%
    map(., ~ .$simulation) %>%
    imap_dfr(., ~ mutate(.x, bound = .y)) %>%
    ggplot() +
    geom_line(aes(x = time, y = logN - 6, colour = bound),
              linetype = 2, size = 1) +
    geom_point(aes(x = time, y = logS),
               data = d_plot, size = 2, shape = 1,
               inherit.aes = FALSE) +
    coord_cartesian(ylim = c(-4, 0)) +
    scale_colour_manual(values = my_cols) +
    theme_clean(base_size = 14) +
    xlab("Treatment time (min)") +
    ylab("Number of log-reductions (·)") +
    theme(legend.position = "none") +
    theme(panel.grid.major = element_line(size = 0))

##

min_time <- 0
max_time <- max(data_10C$time, na.rm = TRUE)

max_count <- 0
# min_count <- min(d_plot$logS)
min_count <- -4

tt <- seq(min_time, max_time, length = 100)
min_temp <- min(d_plot$Temp)
max_temp <- max(d_plot$Temp)

slope <- (max_count - min_count)/(max_temp - min_temp)
intercept <- (min_temp * (max_count-min_count)/(max_temp - min_temp) - min_count)

l1 <- d_plot %>%
    mutate(fake_temp = Temp * slope - intercept) %>%
    geom_line(aes(x = time, y = fake_temp), data = .,
              size = 1)

p1 <- p1 + l1 +
    scale_y_continuous(
                       sec.axis = sec_axis(~(. + intercept)/slope,
                                           name = "Temperature (ºC)"))

## Figure 5

p <- cowplot::plot_grid(p1, p2, p3, labels = "AUTO", nrow = 1)

ggsave(p,
       filename = "Figure_5.png",
       width = 15,
       height  = 5)


