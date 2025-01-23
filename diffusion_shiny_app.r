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
library(WienR)
library(mvtnorm)
library(shiny)

tau <- 0.01
maxT <- 4
nsim <- 20

rt_p <- c(0.1, 0.3, 0.5, 0.7, 0.9)

source("wienr_fit_utils.r")

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Diffusion model"),

    # Sidebar with a slider input for number of bins 
    sidebarPanel(
        sliderInput("v", "Drift rate (v)", min = -5, max = 5, value = 0, step = 0.1),
        sliderInput("a", "Response caution (a)", min = 0, max = 5, value = 2, step = 0.1),
        sliderInput("w", "Response bias (w)", min = 0, max = 1, value = 0.5, step = 0.01),
        sliderInput("t0", "Residual time (t0)", min = 0, max = 1, value = 0.2, step = 0.05),
        hr(),
        sliderInput("sv", "Drift rate variability (sv)", min = 0, max = 5, value = 0, step = 0.1),
        sliderInput("sw", "Boundary variability (sw)", min = 0, max = 1, value = 0, step = 0.05),
        sliderInput("st0", "Residual time variability (st0)", min = 0, max = 1, value = 0, step = 0.05),
        hr(),
        sliderInput("n", "Num. simulated trials", min = 10, max = 1000, value = 100, step = 10),
        numericInput("sim_seed", "Seed for simulation", value = 1, step = 1)
    ),
    mainPanel(
        tabsetPanel(
            tabPanel("Choice and RT predictions", fluidPage(fluidRow(column(width = 8, plotOutput("main_plot", height = "600px")), column(width = 4, plotOutput("qp_plot", height = "600px"))))),
            tabPanel("Manual parameter fitting",
                     fluidPage(fluidRow(
                         column(width = 8, plotOutput("manual_fit_plot", height = "600px")),
                         column(width = 4,
                             textOutput("manual_param_label"),
                             sliderInput("v_fit", "v", min = -5, max = 5, value = 0, step = 0.1),
                             sliderInput("a_fit", "a", min = 0, max = 5, value = 2, step = 0.1),
                             sliderInput("w_fit", "w", min = 0, max = 1, value = 0.5, step = 0.01),
                             sliderInput("t0_fit", "t0", min = 0, max = 1, value = 0.2, step = 0.05),
                             sliderInput("sv_fit", "sv", min = 0, max = 5, value = 0, step = 0.1),
                             sliderInput("sw_fit", "sw", min = 0, max = 1, value = 0, step = 0.05),
                             sliderInput("st0_fit", "st0", min = 0, max = 1, value = 0, step = 0.05)
                         )))
            ),
            tabPanel("Parameter recovery", fluidPage(fluidRow(column(width = 8, plotOutput("recov_plot", height = "600px")), column(width = 4, textOutput("recov_param_label"), tableOutput("recov_pars")))))
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    sim <- reactive({
        set.seed(input$sim_seed)
        
        sim_data <- sampWiener(N = as.integer(input$n), a = as.numeric(input$a), v = as.numeric(input$v + 1E-6), w = as.numeric(input$w), t0 = as.numeric(input$t0), sv = as.numeric(input$sv), sw = as.numeric(input$sw), st0 = as.numeric(input$st0))
        
        sim_rt_quantiles <- tibble(rt = sim_data$q, response = sim_data$response) %>%
            group_by(response) %>%
            reframe(rt_q = quantile(rt, probs = rt_p)) %>%
            mutate(rt_p = rep(rt_p, n() / length(rt_p))) %>%
            complete(response, rt_p, fill = list(rt_q = NA))
        
        sim_p_resp <- tibble(rt = sim_data$q, response = sim_data$response) %>%
            group_by(response) %>%
            summarize(n_resp = n(), .groups = "keep") %>%
            ungroup() %>%
            complete(response, fill = list(n_resp = 0)) %>%
            mutate(p_resp = n_resp / sum(n_resp))
        
        list(sim_data = sim_data, sim_rt_quantiles = sim_rt_quantiles, sim_p_resp = sim_p_resp)
    })
    
    manual_fit <- reactive({
        manual_par <- c("v[1]" = unname(input$v_fit), "a[1]" = unname(input$a_fit), "w[1]" = unname(input$w_fit), "t0[1]" = unname(input$t0_fit), "sv[1]" = unname(input$sv_fit), "sw[1]" = unname(input$sw_fit), "st0[1]" = unname(input$st0_fit))
        
        nll <- -sum(WienerPDF(t = sim()$sim_data$q, response = (sim()$sim_data$response == "upper") + 1, a = input$a_fit, v = input$v_fit, w = input$w_fit, t0 = input$t0_fit, sv = input$sv_fit, sw = input$sw_fit, st0 = input$st0_fit)$logvalue)
        
        fitDF <- expand_grid(response = c("upper", "lower"), rt_p = rt_p) %>% mutate(rt_q = NA, p_resp = NA)
        
        for (i in 1:nrow(fitDF)) {
            fitDF$rt_q[i] <- qWDM(p = fitDF$rt_p[i], response = fitDF$response[i], a = manual_par["a[1]"], v = manual_par["v[1]"], w = manual_par["w[1]"], t0 = manual_par["t0[1]"], sv = manual_par["sv[1]"], sw = manual_par["sw[1]"], st0 = manual_par["st0[1]"])
        }
        
        fitDF$p_resp[fitDF$response == "upper"] <- WienerCDF(t = Inf, response = "upper", a = manual_par["a[1]"], v = manual_par["v[1]"], w = manual_par["w[1]"], t0 = manual_par["t0[1]"], sv = manual_par["sv[1]"], sw = manual_par["sw[1]"], st0 = manual_par["st0[1]"])$value
        fitDF$p_resp[fitDF$response == "lower"] <- WienerCDF(t = Inf, response = "lower", a = manual_par["a[1]"], v = manual_par["v[1]"], w = manual_par["w[1]"], t0 = manual_par["t0[1]"], sv = manual_par["sv[1]"], sw = manual_par["sw[1]"], st0 = manual_par["st0[1]"])$value
        
        sim_fit_data <- full_join(
            full_join(sim()$sim_p_resp, sim()$sim_rt_quantiles) %>% mutate(source = "Simulated"),
            fitDF %>% mutate(source = "Fitted")
        )
        
        list(sim_fit_data = sim_fit_data, nll = nll)
    })
    
    recov_fit <- reactive({
        init_par <- c("v[1]" = unname(input$v), "a[1]" = unname(input$a), "w[1]" = unname(input$w), "t0[1]" = unname(input$t0))
        if (input$sv > 0) init_par <- c(init_par, "sv[1]" = unname(input$sv))
        if (input$sw > 0) init_par <- c(init_par, "sw[1]" = unname(input$sw))
        if (input$st0 > 0) init_par <- c(init_par, "st0[1]" = unname(input$st0))
        
        fit <- fit_wienr(rt = sim()$sim_data$q, response = (sim()$sim_data$response == "upper") + 1, fit_sv = (input$sv > 0), fit_sw = (input$sw > 0), fit_st0 = (input$st0 > 0), init_par = init_par)
        
        nll <- fit$value
        fit_pars <- fit$par
        
        if (input$sv <= 0) fit_pars <- c(fit_pars, "sv[1]" = 0)
        if (input$sw <= 0) fit_pars <- c(fit_pars, "sw[1]" = 0)
        if (input$st0 <= 0) fit_pars <- c(fit_pars, "st0[1]" = 0)
        
        fitDF <- expand_grid(response = c("upper", "lower"), rt_p = rt_p) %>% mutate(rt_q = NA, p_resp = NA)
        
        for (i in 1:nrow(fitDF)) {
            fitDF$rt_q[i] <- qWDM(p = fitDF$rt_p[i], response = fitDF$response[i], a = fit_pars["a[1]"], v = fit_pars["v[1]"], w = fit_pars["w[1]"], t0 = fit_pars["t0[1]"], sv = fit_pars["sv[1]"], sw = fit_pars["sw[1]"], st0 = fit_pars["st0[1]"])
        }
        
        fitDF$p_resp[fitDF$response == "upper"] <- WienerCDF(t = Inf, response = "upper", a = fit_pars["a[1]"], v = fit_pars["v[1]"], w = fit_pars["w[1]"], t0 = fit_pars["t0[1]"], sv = fit_pars["sv[1]"], sw = fit_pars["sw[1]"], st0 = fit_pars["st0[1]"])$value
        fitDF$p_resp[fitDF$response == "lower"] <- WienerCDF(t = Inf, response = "lower", a = fit_pars["a[1]"], v = fit_pars["v[1]"], w = fit_pars["w[1]"], t0 = fit_pars["t0[1]"], sv = fit_pars["sv[1]"], sw = fit_pars["sw[1]"], st0 = fit_pars["st0[1]"])$value
        
        sim_fit_data <- full_join(
            full_join(sim()$sim_p_resp, sim()$sim_rt_quantiles) %>% mutate(source = "Simulated"),
            fitDF %>% mutate(source = "Fitted")
        )
        
        par_mat <- cbind(
            c(input$v, input$a, input$w, input$t0, input$sv, input$sw, input$st0),
            fit_pars[c("v[1]", "a[1]", "w[1]", "t0[1]", "sv[1]", "sw[1]", "st0[1]")]
        )
        
        rownames(par_mat) <- c("v", "a", "w", "t0", "sv", "sw", "st0")
        colnames(par_mat) <- c("Simulation", "Best fit")
        
        list(par_mat = par_mat, sim_fit_data = sim_fit_data, nll = nll)
    })
    
    output$main_plot <- renderPlot({
        set.seed(12222)
        
        t <- seq(0, maxT, by=tau)
        
        mu <- input$v * t
        Sigma <- outer(t, t, FUN=pmin) * (1 + input$sv^2)
        
        b_upper <- input$a * (1 - input$w)
        b_lower <- -input$a * input$w
        
        x <- rmvnorm(n=nsim, mean=mu, sigma=Sigma, method='chol')
        
        plotDF <- c()
        
        for (i in 1:nrow(x)) {
            x[i,] <- x[i,] + runif(n = 1, min = -0.5 * input$sw * input$a, max = 0.5 * input$sw * input$a)
            u <- min(c(which(x[i,] > b_upper), ncol(x)))
            l <- min(c(which(x[i,] < b_lower), ncol(x)))
            stopTime <- max(1, min(u, l) - 1)
            if (stopTime + 1 <= ncol(x)) {
                x[i, (stopTime + 1):ncol(x)] <- NA
            }
            
            plotDF <- rbind(
                plotDF,
                tibble(t = t + runif(n = 1, min = 0, max = input$st0), x = x[i,], index = i)
            )
        }
        
        predDF <- tibble(t = t) %>%
            filter(t > 0) %>%
            mutate(
                d_upper = WienerPDF(t = t, response = "upper", a = input$a, v = input$v, w = input$w, t0 = 0, sv = input$sv, sw = input$sw, st0 = input$st0)$value,
                d_lower = WienerPDF(t = t, response = "lower", a = input$a, v = input$v, w = input$w, t0 = 0, sv = input$sv, sw = input$sw, st0 = input$st0)$value
            )
        
        upperPlot <- predDF %>%
            ggplot(aes(x = t + input$t0, y = d_upper)) +
            geom_area(alpha = 0.5, color = "darkblue", fill = "darkblue") +
            coord_cartesian(xlim = c(0, max(t)), ylim = c(0, max(c(predDF$d_upper, predDF$d_lower)))) +
            labs(x = NULL, y = "P(upper at time t)")
        
        lowerPlot <- predDF %>%
            ggplot(aes(x = t + input$t0, y = d_lower)) +
            geom_area(alpha = 0.5, color = "darkred", fill = "darkred") +
            coord_cartesian(xlim = c(0, max(t)), ylim = c(0, max(c(predDF$d_upper, predDF$d_lower)))) +
            labs(x = NULL, y = "P(lower at time t)")
        
        trajPlot <- plotDF %>%
            filter(!is.na(x)) %>%
            ggplot(aes(x = t + input$t0, y = x, color = index, group = index)) +
            geom_line(alpha = 0.5) +
            geom_hline(yintercept = c(b_lower, b_upper), color = "black", linetype = "solid") +
            annotate(geom = "segment", x = input$t0, y = 0, xend = input$t0 + 1, yend = input$v, color = "black", linewidth = 1, arrow = arrow(length = unit(0.03, "npc"))) +
            coord_cartesian(xlim = c(0, max(t)), ylim = c(-input$a, input$a)) +
            scale_color_viridis_c(guide = "none") +
            labs(x = "Time", y = "Accumulated evidence")
        
        upperPlot / trajPlot / lowerPlot + plot_layout(heights = c(0.6, 1, 0.6))
    }, res = 100)
    
    output$qp_plot <- renderPlot({
        plotDF <- expand_grid(response = c("upper", "lower"), rt_p = rt_p) %>% mutate(rt_q = NA, p_resp = NA)
        
        for (i in 1:nrow(plotDF)) {
            plotDF$rt_q[i] <- qWDM(p = plotDF$rt_p[i], response = plotDF$response[i], a = input$a, v = input$v, w = input$w, t0 = input$t0, sv = input$sv, sw = input$sw, st0 = input$st0)
        }
        
        plotDF$p_resp[plotDF$response == "upper"] <- WienerCDF(t = Inf, response = "upper", a = input$a, v = input$v, w = input$w, t0 = input$t0, sv = input$sv, sw = input$sw, st0 = input$st0)$value
        plotDF$p_resp[plotDF$response == "lower"] <- WienerCDF(t = Inf, response = "lower", a = input$a, v = input$v, w = input$w, t0 = input$t0, sv = input$sv, sw = input$sw, st0 = input$st0)$value
        
        plotDF %>%
            ggplot(aes(x = p_resp, y = rt_q, shape = response)) +
            geom_line(aes(group = rt_p), alpha = 0.5) +
            geom_point() +
            expand_limits(x = c(0, 1)) +
            labs(x = "Response probability", y = "RT Quantile (s)", shape = "Response")
    }, res = 100)
    
    output$manual_fit_plot <- renderPlot({
        manual_fit()$sim_fit_data %>%
            ggplot(aes(x = p_resp, y = rt_q, color = source, shape = response)) +
            geom_line(aes(group = interaction(rt_p, source)), alpha = 0.5) +
            geom_point() +
            annotate(geom = "label", x = 0, y = Inf, label = paste0("NLL = ", round(manual_fit()$nll, 3)), hjust = "inward", vjust = "inward") +
            expand_limits(x = c(0, 1)) +
            labs(x = "Response probability", y = "RT Quantile (s)", color = NULL, shape = "Response")
    }, res = 100)
    
    output$recov_plot <- renderPlot({
        recov_fit()$sim_fit_data %>%
            ggplot(aes(x = p_resp, y = rt_q, color = source, shape = response)) +
            geom_line(aes(group = interaction(rt_p, source)), alpha = 0.5) +
            geom_point() +
            annotate(geom = "label", x = 0, y = Inf, label = paste0("NLL = ", round(recov_fit()$nll, 3)), hjust = "inward", vjust = "inward") +
            expand_limits(x = c(0, 1)) +
            labs(x = "Response probability", y = "RT Quantile (s)", color = NULL, shape = "Response")
    }, res = 100)
    
    output$recov_pars <- renderTable(recov_fit()$par_mat, striped = TRUE, hover = TRUE, bordered = TRUE, rownames = TRUE)
    
    output$recov_param_label <- renderText({"These are the \"best-fitting\" parameters that minimize the negative log-likelihood (NLL)."})
    output$manual_param_label <- renderText({"Adjust the parameter values below to try to reproduce the simulated data."})
}

# Run the application 
shinyApp(ui = ui, server = server)
