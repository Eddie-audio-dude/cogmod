#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(tidyverse)
library(patchwork)
library(khroma)
library(deSolve)
library(shiny)

dt <- 1
max_t <- 500

scri_system <- function(t, y, parms) {
    # First half of "y" are the salience units, second half are the identification units
    N <- length(y) / 2
    vis <- y[1:N]
    id <- y[1:N + N]
    
    # Localization signal (x) depends only on the current time
    loc <- parms$strength_loc * dgamma(t, shape = parms$loc_shape, rate = parms$loc_rate)
    
    # Derivative of the salience units (v)
    d_vis <- (parms$saturation_vis - vis) * (parms$baseline + loc + id) -
        vis * (parms$leak_vis + parms$ff_loc * c(parms$ff_loc_sum %*% loc) + parms$ff_id * c(parms$ff_id_sum %*% id) + parms$lat_vis * c(parms$lat_vis_sum %*% vis))
    
    # Derivative of the identification units (z)
    d_id <- (parms$strength_id - id) * pgamma(t, shape = parms$loc_shape * (1 + parms$id_delay), rate = parms$loc_rate) * vis^parms$recurrent_gating -
        id * (parms$leak_id + parms$lat_id * c(parms$lat_id_sum %*% id))
    
    # Return a listified vector of both sets of derivatives
    return(list(c(d_vis, d_id)))
}

scri <- function(times = seq(0, 500),
                 strength_loc = rep(0.539, 8),
                 strength_id = 0.023 * c(1, rep(0.222, 7)),
                 loc_peak = 130, loc_spread = 35,
                 leak_vis = 0.328, leak_id = 0.071,
                 ff_loc = 0.168, ff_id = 20.689,
                 lat_vis = 1.217, lat_id = 0.445,
                 ff_loc_spread = Inf, ff_id_spread = Inf,
                 lat_vis_spread = Inf, lat_id_spread = Inf,
                 id_delay = 0.108, baseline = 0.004,
                 recurrent_gating = TRUE, saturation_vis = 1,
                 odeMethod = "bdf", return_format = "df", ...) {
    if (length(strength_loc) != length(strength_id)) {
        stop(paste0("Length of localization strength vector (", length(strength_loc), ") not equal to length of identification strength vector (", length(strength_id), ")."))
    }
    
    N <- length(strength_loc)
    
    # Convert peak and spread of localization signal to shape and rate parameters
    loc_rate <- (loc_peak + sqrt(loc_peak^2 + 4 * loc_spread^2)) / (2 * loc_spread^2)
    loc_shape <- 1 + loc_peak * loc_rate
    
    # Compute the distance between each receptive field location, assuming they are equidistant on a circle of radius 1 (so distances are relative to the radius of the search array)
    angle <- seq(0, 2 * pi, length.out=N+1)[1:N]
    distance <- 2 * sin(outer(angle, angle, FUN=function(a, b) abs(a - b)) / 2)
    
    # The following matrices represent the strength of feedforward (ff) and lateral (lat) inhibition between each unit as a function of their distance
    ff_loc_sum <- exp(-distance^2 / (2 * ff_loc_spread^2)) * (1 - diag(N))
    ff_id_sum <- exp(-distance^2 / (2 * ff_id_spread^2)) * (1 - diag(N))
    
    lat_vis_sum <- exp(-distance^2 / (2 * lat_vis_spread^2)) * (1 - diag(N))
    lat_id_sum <- exp(-distance^2 / (2 * lat_id_spread^2)) * (1 - diag(N))
    
    # Compute initial values for each salience unit
    if (baseline > 0) {
        if (lat_vis > 0) {
            A <- -lat_vis * rowSums(lat_vis_sum)
            B <- -(baseline + leak_vis)
            C <- baseline
            
            init_v <- (-B - sqrt(B^2 - 4 * A * C)) / (2 * A)
        } else {
            init_v <- rep(baseline / (baseline + leak_vis), N)
        }
    } else {
        init_v <- rep(0, N)
    }
    
    # Prepare arguments for "ode" function
    y0 <- c(init_v, rep(0, N))
    parms <- list(
        strength_loc = strength_loc,
        strength_id = strength_id,
        loc_shape = loc_shape,
        loc_rate = loc_rate,
        leak_vis = leak_vis,
        leak_id = leak_id,
        ff_loc = ff_loc,
        ff_id = ff_id,
        lat_vis = lat_vis,
        lat_id = lat_id,
        ff_loc_sum = ff_loc_sum,
        ff_id_sum = ff_id_sum,
        lat_vis_sum = lat_vis_sum,
        lat_id_sum = lat_id_sum,
        id_delay = id_delay,
        baseline = baseline,
        recurrent_gating = recurrent_gating,
        saturation_vis = saturation_vis
    )
    
    # Call the `ode` function to compute activations at each time point
    solved <- ode(y = y0, times = times, parms = parms, func = scri_system, method = odeMethod, ...)
    
    loc <- t(strength_loc %o% dgamma(solved[, "time"], shape = loc_shape, rate = loc_rate))
    vis <- solved[,1 + 1:N]
    id <- solved[,1 + N + 1:N]
    
    dimnames(loc) <- dimnames(vis) <- dimnames(id) <- list(
        "t" = solved[, "time"],
        "array_position" = 1:N
    )
    
    if (return_format == "df") {
        to_return <- rbind(
            cbind(unit_type = "localization", array2DF(loc, responseName = "act")),
            cbind(unit_type = "identification", array2DF(id, responseName = "act")),
            cbind(unit_type = "salience", array2DF(vis, responseName = "act"))
        ) %>%
            mutate(
                t = as.numeric(t),
                array_position = factor(array_position),
                unit_type = factor(unit_type, levels = c("localization", "identification", "salience"))
            )
    } else {
        to_return <- list(
            loc = loc,
            vis = vis,
            id = id
        )
    }
    
    return(to_return)
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Salience by Competitive and Recurrent Interactions"),

    # Sidebar with a slider input for number of bins 
    sidebarPanel(
        sliderInput("N", "Set size", min = 2, max = 8, step = 1, value = 4),
        checkboxInput("recurrent_gating", "Salience gates identification", value = TRUE),
        hr(),
        sliderInput("baseline", "Tonic baseline excitation (b)", min = 0, max = 0.1, value = 0.01, step = 0.01),
        sliderInput("strength_loc", "Localization strength (chi)", min = 0, max = 1, value = 0.5, step = 0.01),
        sliderInput("strength_id", "Identification strength (eta)", min = 0, max = 0.5, value = c(0.01, 0.03), step = 0.01),
        sliderInput("loc_peak", "Localization peak time (omega_p)", min = 0, max = 500, value = 130, step = 10),
        sliderInput("loc_spread", "Localization spread (omega_s)", min = 0, max = 100, value = 30, step = 5),
        sliderInput("id_delay", "Identification delay (kappa)", min = 0, max = 3, value = 0, step = 0.1),
        hr(),
        sliderInput("leak_vis", "Salience unit leakage (lambda_v)", min = 0, max = 1, value = 0.3, step = 0.01),
        sliderInput("leak_id", "Identification unit leakage (lambda_z)", min = 0, max = 1, value = 0.05, step = 0.01),
        hr(),
        sliderInput("ff_loc", "Localization feedforward inhibition (alpha_x)", min = 0, max = 1, value = 0.1, step = 0.05),
        sliderInput("ff_id", "Identification feedforward inhibition (alpha_z)", min = 0, max = 1, value = 0.1, step = 0.05),
        sliderInput("lat_vis", "Salience lateral inhibition (beta_v)", min = 0, max = 1, value = 0.1, step = 0.05),
        sliderInput("lat_id", "Identification lateral inhibition (beta_z)", min = 0, max = 1, value = 0.1, step = 0.05)
    ),
    mainPanel(
        plotOutput("scri_plot", width = "600px", height = "700px")
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$scri_plot <- renderPlot({
        if (input$N == 2) {
            input_mask <- c(1, 0, 0, 0, 1, 0, 0, 0)
        } else if (input$N == 3) {
            input_mask <- c(1, 0, 1, 0, 1, 0, 0, 0)
        } else if (input$N == 4) {
            input_mask <- c(1, 0, 1, 0, 1, 0, 1, 0)
        } else if (input$N == 5) {
            input_mask <- c(1, 1, 1, 0, 1, 0, 1, 0)
        } else if (input$N == 6) {
            input_mask <- c(1, 1, 1, 0, 1, 0, 1, 1)
        } else if (input$N == 7) {
            input_mask <- c(1, 1, 1, 1, 1, 0, 1, 1)
        } else if (input$N == 8) {
            input_mask <- c(1, 1, 1, 1, 1, 1, 1, 1)
        }
        
        sim <- scri(
            times = seq(0, max_t, by = dt),
            strength_loc = input$strength_loc * input_mask,
            strength_id = min(input$strength_id) * input_mask + (max(input$strength_id) - min(input$strength_id)) * c(1, rep(0, 7)),
            loc_peak = input$loc_peak, loc_spread = input$loc_spread,
            leak_vis = input$leak_vis, leak_id = input$leak_id,
            ff_loc = input$ff_loc, ff_id = input$ff_id,
            lat_vis = input$lat_vis, lat_id = input$lat_id,
            id_delay = input$id_delay, baseline = input$baseline,
            recurrent_gating = input$recurrent_gating
        )
        
        array_plot <- tibble(
            array_position = as.character(1:8),
            x = sin(seq(0, 2 * pi, length.out = 9)[1:8]),
            y = cos(seq(0, 2 * pi, length.out = 9)[1:8]),
            is_on = input_mask > 0,
            is_target = c(TRUE, rep(FALSE, 7))
        ) %>%
            ggplot(aes(x = x, y = y, label = array_position, color = array_position, fill = is_target, size = is_on)) +
            geom_label() +
            coord_equal(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)) +
            scale_color_okabeito(guide = "none") +
            scale_fill_manual(values = c("TRUE" = "#cccccc", "FALSE" = "white"), guide = "none") +
            scale_size_manual(values = c("TRUE" = 6, "FALSE" = 3), guide = "none") +
            labs(Title = "Search array") +
            theme_void()
        
        dyn_plot <- sim %>%
            ggplot(aes(x = t, y = act, color = array_position)) +
            geom_line() +
            expand_limits(y = 0) +
            scale_color_okabeito(guide = "none") +
            facet_wrap("unit_type", ncol = 1, scales = "free_y") +
            labs(x = "Time from array onset (ms)", y = "Activation", color = "Array\nposition")
        
        dyn_plot + array_plot + plot_layout(nrow = 1, widths = c(2, 1))
    }, res = 100)
}

# Run the application 
shinyApp(ui = ui, server = server)
