
## Functions for model fitting

two_lines <- function(time, pars) {
    
    pars <- as.list(pars)
    logN_crit <- pars$logN0 - pars$t_crit/pars$D1
    
    ifelse(time < pars$t_crit,
           pars$logN0 - time/pars$D1,
           logN_crit - (time-pars$t_crit)/pars$D2
    )
    
}

one_line <- function(time, pars) {
    
    pars <- as.list(pars)
    pars$logN0 - time/pars$D
    
}

get_residuals <- function(this_p, data, model = two_lines) {
    
    this_p <- as.list(this_p)
    
    pred <- tibble(time = unique(data$time),
                   logN = model(time, this_p)
    ) %>%
        as.data.frame()
    
    data <- select(data, time, logN) %>%
        as.data.frame()
    
    modCost(model = pred,
            obs = data)
    
}

make_plot <- function(d, models) {
    
    my_cols <- wes_palette("Moonrise2", 2)
    
    p <- d %>%
        ggplot() +
        geom_point(aes(x = time, y = logN, colour = preadapted, 
                       shape = preadapted), size = 3)
    
    pars1 <- coef(models[[1]])
    
    if ("D" %in% names(pars1)) {
        pred1 <- tibble(times = seq(0, max(d$time), length = 100),
                        logN = one_line(times, pars1),
                        preadapted = FALSE)
    } else {
        pred1 <- tibble(times = seq(0, max(d$time), length = 100),
                        logN = two_lines(times, pars1),
                        preadapted = FALSE)
    }
    
    pars2 <- coef(models[[2]])
    
    if ("D" %in% names(pars2)) {
        pred2 <- tibble(times = seq(0, max(d$time), length = 100),
                        logN = one_line(times, pars2),
                        preadapted = TRUE)
    } else {
        pred2 <- tibble(times = seq(0, max(d$time), length = 100),
                        logN = two_lines(times, pars2),
                        preadapted = TRUE)
    }
    
    
    
    l <- bind_rows(pred1, pred2) %>%
        geom_line(aes(times, logN, colour = preadapted,
                      linetype = preadapted), 
                  data = .,
                  size = 1)
    
    p + l + 
        xlab("Treatment time (min)") + 
        ylab("Microbial concentration (log CFU/mL)") +
        scale_shape_manual(values = c(1, 16)) +
        scale_linetype_manual(values = c(2, 1)) +
        scale_color_manual(values = my_cols) +
        # theme_minimal(base_size = 14) +
        # geom_rangeframe() +
        theme_clean(base_size = 16) +
        theme(legend.position = "none") +
        coord_cartesian(ylim = c(1, 7))
    
}



make_plot2 <- function(d, models) {
    
    my_cols <- wes_palette("Moonrise2", 4)[c(4,1)]
    
    p <- d %>%
        ggplot() +
        geom_point(aes(x = time, y = logN, colour = strain,
                       shape = strain), size = 3)
    
    pars1 <- coef(models[[1]])
    
    if ("D" %in% names(pars1)) {
        pred1 <- tibble(times = seq(0, max(d$time), length = 100),
                        logN = one_line(times, pars1),
                        strain = "wild")
    } else {
        pred1 <- tibble(times = seq(0, max(d$time), length = 100),
                        logN = two_lines(times, pars1),
                        strain = "wild")
    }
    
    pars2 <- coef(models[[2]])
    
    if ("D" %in% names(pars2)) {
        pred2 <- tibble(times = seq(0, max(d$time), length = 100),
                        logN = one_line(times, pars2),
                        strain = "mutant")
    } else {
        pred2 <- tibble(times = seq(0, max(d$time), length = 100),
                        logN = two_lines(times, pars2),
                        strain = "mutant")
    }
    
    
    
    l <- bind_rows(pred1, pred2) %>%
        geom_line(aes(times, logN, colour = strain, linetype = strain), 
                  size = 1, data = .)
    
    p + l + 
        xlab("Treatment time (min)") + 
        ylab("Microbial concentration (log CFU/mL)") +
        scale_shape_manual(values = c(1, 16)) +
        scale_linetype_manual(values = c(2, 1)) +
        scale_color_manual(values = my_cols) +
        # theme_minimal(base_size = 14) +
        # geom_rangeframe() +
        theme_clean(base_size = 16) +
        theme(legend.position = "none") +
        coord_cartesian(ylim = c(1, 7))
    
    
}

make_plot3 <- function(d, models) {
    
    my_cols <- wes_palette("Moonrise2", 2)
    
    p <- d %>%
        ggplot() +
        geom_point(aes(x = time, y = logN, colour = preadapted,
                       shape = preadapted), size = 3)
    
    pars1 <- coef(models[[1]])
    
    if ("D" %in% names(pars1)) {
        pred1 <- tibble(times = seq(0, max(d$time), length = 100),
                        logN = one_line(times, pars1),
                        preadapted = FALSE)
    } else {
        pred1 <- tibble(times = seq(0, max(d$time), length = 100),
                        logN = two_lines(times, pars1),
                        preadapted = FALSE)
    }
    
    pars2 <- coef(models[[2]])
    
    if ("D" %in% names(pars2)) {
        pred2 <- tibble(times = seq(0, max(d$time), length = 100),
                        logN = one_line(times, pars2),
                        preadapted = TRUE)
    } else {
        pred2 <- tibble(times = seq(0, max(d$time), length = 100),
                        logN = two_lines(times, pars2),
                        preadapted = TRUE)
    }
    
    
    
    l <- bind_rows(pred1, pred2) %>%
        geom_line(aes(times, logN, colour = preadapted,
                      linetype = preadapted), 
                  data = .,
                  size = 1)
    
    p + l + 
        xlab("Treatment time (min)") + 
        ylab("Microbial concentration (log CFU/mL)") +
        scale_shape_manual(values = c(1, 16)) +
        scale_linetype_manual(values = c(2, 1)) +
        scale_color_manual(values = my_cols) +
        # theme_minimal(base_size = 14) +
        # geom_rangeframe() +
        theme_clean(base_size = 16) +
        theme(legend.position = "none") +
        coord_cartesian(ylim = c(0, 6.5))
    
}
