library(WienR)
library(ucminf)

ezddm <- function(propCorrect, rtCorrectVariance_seconds, rtCorrectMean_seconds, nTrials, s = 1) {
    # This function is adapted from the function of the same name in the `fddm` package
    s2 <- s^2
    v <- as.numeric(NA)
    a <- as.numeric(NA)
    Ter <- as.numeric(NA)
    propCorrect <- (2 * nTrials * propCorrect) / (2 * nTrials + 1) + 1 / (4 * nTrials * nTrials)
    L <- stats::qlogis(propCorrect)
    x <- L * (L * propCorrect^2 - L * propCorrect + propCorrect - 0.5) / rtCorrectVariance_seconds
    v <- sign(propCorrect - 0.5) * s * x^(1/4)
    a <- s2 * stats::qlogis(propCorrect)/v
    y <- -v * a/s2
    MDT <- (a/(2 * v)) * (1 - exp(y))/(1 + exp(y))
    Ter <- rtCorrectMean_seconds - MDT
    return(c("a" = a, "v" = v, "Ter" = Ter))
}

qWDM <- function(p, response, ...) {
    p_resp <- WienerCDF(Inf, response = response, ...)$value
    
    res <- try(uniroot(f = function(t) WienerCDF(t, response = response, ...)$value / p_resp - p, interval = c(0, 5), f.lower = -p, extendInt = "upX"))
    
    if (class(res) == "try-error") return(NA)
    
    return(res$root)
}

fit_wienr <- function(rt, response, fit_sv = FALSE, fit_sw = FALSE, fit_st0 = FALSE, optim_control = list(), init_par = NULL, drift_index = NULL, bound_index = NULL, resid_index = NULL, sv_index = NULL, sw_index = NULL, st0_index = NULL, ...) {
    if (!is.factor(response)) response <- as.factor(response)
    response_numeric <- as.numeric(response)
    
    par_names <- c()
    
    if (is.null(drift_index)) {
        drift_index <- rep(1, length(rt))
        n_drift <- 1
    } else {
        n_drift <- max(drift_index)
    }
    
    if (is.null(bound_index)) {
        bound_index <- rep(1, length(rt))
        n_bound <- 1
    } else {
        n_bound <- max(bound_index)
    }
    
    if (is.null(resid_index)) {
        resid_index <- rep(1, length(rt))
        n_resid <- 1
    } else {
        n_resid <- max(resid_index)
    }
    
    if (is.null(sv_index)) {
        sv_index <- rep(1, length(rt))
        n_sv <- 1
    } else {
        n_sv <- max(sv_index)
    }
    
    if (is.null(sw_index)) {
        sw_index <- rep(1, length(rt))
        n_sw <- 1
    } else {
        n_sw <- max(sw_index)
    }
    
    if (is.null(st0_index)) {
        st0_index <- rep(1, length(rt))
        n_st0 <- 1
    } else {
        n_st0 <- max(st0_index)
    }
    
    par_names <- c(paste0("a[", 1:n_bound, "]"), paste0("v[", 1:n_drift, "]"), paste0("w[", 1:n_bound, "]"), paste0("t0[", 1:n_resid, "]"))
    
    if (fit_sv) par_names <- c(par_names, paste0("sv[", 1:n_sv, "]"))
    if (fit_sw) par_names <- c(par_names, paste0("sw[", 1:n_sw, "]"))
    if (fit_st0) par_names <- c(par_names, paste0("st0[", 1:n_st0, "]"))
    
    init_to_use <- rep(NA, length(par_names))
    names(init_to_use) <- par_names
    
    ez_init <- ezddm(mean(response_numeric == 2), var(rt[response_numeric == 2]), mean(rt[response_numeric == 2]), length(response_numeric))
    
    init_to_use[startsWith(par_names, "a[")] <- unname(ez_init["a"])
    init_to_use[startsWith(par_names, "v[")] <- unname(ez_init["v"])
    init_to_use[startsWith(par_names, "w[")] <- 0.5
    init_to_use[startsWith(par_names, "t0[")] <- unname(min(0.99 * min(rt), ez_init["Ter"]))
    
    if (fit_sv) init_to_use[startsWith(par_names, "sv[")] <- 0
    if (fit_sw) init_to_use[startsWith(par_names, "sw[")] <- 0
    if (fit_st0) init_to_use[startsWith(par_names, "st0[")] <- 0
    
    if (!is.null(init_par)) {
        overlap <- intersect(names(init_to_use), names(init_par))
        init_to_use[overlap] <- init_par[overlap]
    }
    
    init_par <- init_to_use
    
    lower <- rep(-Inf, length(par_names))
    upper <- rep(Inf, length(par_names))
    
    names(lower) <- names(upper) <- par_names
    
    lower[startsWith(par_names, "a[")] <- 0
    lower[startsWith(par_names, "w[")] <- 0
    lower[startsWith(par_names, "t0[")] <- 0
    lower[startsWith(par_names, "sv[")] <- 0
    lower[startsWith(par_names, "sw[")] <- 0
    lower[startsWith(par_names, "st0[")] <- 0
    
    upper[startsWith(par_names, "w[")] <- 1
    upper[startsWith(par_names, "t0[")] <- unname(tapply(rt, INDEX = resid_index, FUN = min))
    upper[startsWith(par_names, "sw[")] <- 1
    
    neg_log_likelihood <- function(par, rt, response, bound_index, drift_index, resid_index, sv_index = NULL, sw_index = NULL, st0_index = NULL, ...) {
        a <- par[paste0("a[", bound_index, "]")]
        v <- par[paste0("v[", drift_index, "]")]
        w <- par[paste0("w[", bound_index, "]")]
        t0 <- par[paste0("t0[", resid_index, "]")]
        
        if (is.na(par["sv[1]"])) {
            sv <- 0
        } else {
            sv <- par[paste0("sv[", sv_index, "]")]
        }
        
        if (is.na(par["sw[1]"])) {
            sw <- 0
        } else {
            sw <- par[paste0("sw[", sw_index, "]")]
        }
        
        if (is.na(par["st0[1]"])) {
            st0 <- 0
        } else {
            st0 <- par[paste0("st0[", st0_index, "]")]
        }
        
        eval_pdf <- try(WienerPDF(t = rt, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, ...), silent = TRUE)
        
        if (any(class(eval_pdf) == "try-error")) return(Inf)
        
        return(-sum(eval_pdf$logvalue))
    }
    
    gradient <- function(par, rt, response, bound_index, drift_index, resid_index, sv_index = NULL, sw_index = NULL, st0_index = NULL, ...) {
        a <- par[paste0("a[", bound_index, "]")]
        v <- par[paste0("v[", drift_index, "]")]
        w <- par[paste0("w[", bound_index, "]")]
        t0 <- par[paste0("t0[", resid_index, "]")]
        
        if (is.na(par["sv[1]"])) {
            sv <- 0
            use_sv <- FALSE
        } else {
            sv <- par[paste0("sv[", sv_index, "]")]
            use_sv <- TRUE
        }
        
        if (is.na(par["sw[1]"])) {
            sw <- 0
            use_sw <- FALSE
        } else {
            sw <- par[paste0("sw[", sw_index, "]")]
            use_sw <- TRUE
        }
        
        if (is.na(par["st0[1]"])) {
            st0 <- 0
            use_st0 <- FALSE
        } else {
            st0 <- par[paste0("st0[", st0_index, "]")]
            use_st0 <- TRUE
        }
        
        eval_grad <- try(gradWienerPDF(t = rt, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, ...), silent = TRUE)
        
        if (any(class(eval_grad) == "try-error")) return(rep(NaN, length(par)))
        
        eval_pdf <- try(WienerPDF(t = rt, response = response, a = a, v = v, w = w, t0 = t0, sv = sv, sw = sw, st0 = st0, ...), silent = TRUE)
        if (any(class(eval_pdf) == "try-error")) return(rep(NaN, length(par)))
        
        # Derivative of log(f(x)) is f'(x) / f(x)
        
        grad <- rep(NaN, length(par))
        names(grad) <- names(par)
        
        for (i in 1:max(bound_index)) {
            grad[paste0("a[", i, "]")] <- sum(eval_grad$deriv[bound_index == i, "da"] / eval_pdf$value[bound_index == i])
            grad[paste0("w[", i, "]")] <- sum(eval_grad$deriv[bound_index == i, "dw"] / eval_pdf$value[bound_index == i])
        }
        
        if (use_sw) {
            for (i in 1:max(sw_index)) {
                grad[paste0("sw[", i, "]")] <- sum(eval_grad$deriv[sw_index == i, "dsw"] / eval_pdf$value[sw_index == i])
            }
        }
        
        for (i in 1:max(drift_index)) {
            grad[paste0("v[", i, "]")] <- sum(eval_grad$deriv[drift_index == i, "dv"] / eval_pdf$value[drift_index == i])
        }
        
        if (use_sv) {
            for (i in 1:max(sv_index)) {
                grad[paste0("sv[", i, "]")] <- sum(eval_grad$deriv[sv_index == i, "dsv"] / eval_pdf$value[sv_index == i])
            }
        }
        
        for (i in 1:max(resid_index)) {
            grad[paste0("t0[", i, "]")] <- sum(eval_grad$deriv[resid_index == i, "dt0"] / eval_pdf$value[resid_index == i])
        }
        
        if (use_st0) {
            for (i in 1:max(st0_index)) {
                grad[paste0("st0[", i, "]")] <- sum(eval_grad$deriv[st0_index == i, "dst0"] / eval_pdf$value[st0_index == i])
            }
        }
        
        return(-grad)
    }
    
    if (!is.na(init_par["sw[1]"]) | !is.na(init_par["st0[1]"])) {
        init_par0 <- init_par[!(startsWith(names(init_par), "sw") | startsWith(names(init_par), "st0"))]
        # lower0 <- lower[!(startsWith(names(lower), "sw") | startsWith(names(lower), "st0"))]
        # upper0 <- upper[!(startsWith(names(upper), "sw") | startsWith(names(upper), "st0"))]
        
        fit0 <- try(optim(
            par = init_par0,
            fn = neg_log_likelihood,
            gr = gradient,
            # lower = lower0,
            # upper = upper0,
            method = "Nelder-Mead",
            control = c(optim_control, maxit = 10000),
            rt = rt,
            response = response_numeric,
            bound_index = bound_index,
            drift_index = drift_index,
            resid_index = resid_index,
            sv_index = sv_index,
            sw_index = sw_index,
            st0_index = st0_index
        ), silent = TRUE)
        
        if (class(fit0) != "try-error") {
            overlap <- intersect(names(init_par), names(init_par0))
            init_par[overlap] <- fit0$par[overlap]
        }
    }
    
    fit1 <- optim(
        par = init_par,
        fn = neg_log_likelihood,
        gr = gradient,
        method = "Nelder-Mead",
        control = c(optim_control, maxit = 10000),
        rt = rt,
        response = response_numeric,
        bound_index = bound_index,
        drift_index = drift_index,
        resid_index = resid_index,
        sv_index = sv_index,
        sw_index = sw_index,
        st0_index = st0_index
    )
    
    fit2 <- ucminf(
        par = fit1$par,
        fn = neg_log_likelihood,
        gr = gradient,
        control = optim_control,
        rt = rt,
        response = response_numeric,
        bound_index = bound_index,
        drift_index = drift_index,
        resid_index = resid_index,
        sv_index = sv_index,
        sw_index = sw_index,
        st0_index = st0_index
    )
    
    return(fit2)
}

qp_fit <- function(rt, response, par = NULL, rt_p = c(0.1, 0.3, 0.5, 0.7, 0.9), drift_index = NULL, bound_index = NULL, resid_index = NULL, sv_index = NULL, sw_index = NULL, st0_index = NULL) {
    if (is.null(drift_index)) {
        drift_index <- rep(1, length(rt))
        n_drift <- 1
    } else {
        n_drift <- max(drift_index)
    }
    
    if (is.null(bound_index)) {
        bound_index <- rep(1, length(rt))
        n_bound <- 1
    } else {
        n_bound <- max(bound_index)
    }
    
    if (is.null(resid_index)) {
        resid_index <- rep(1, length(rt))
        n_resid <- 1
    } else {
        n_resid <- max(resid_index)
    }
    
    if (is.null(sv_index)) {
        sv_index <- rep(1, length(rt))
        n_sv <- 1
    } else {
        n_sv <- max(sv_index)
    }
    
    if (is.null(sw_index)) {
        sw_index <- rep(1, length(rt))
        n_sw <- 1
    } else {
        n_sw <- max(sw_index)
    }
    
    if (is.null(st0_index)) {
        st0_index <- rep(1, length(rt))
        n_st0 <- 1
    } else {
        n_st0 <- max(st0_index)
    }
    
    obs_rt_quantiles <- tibble(rt = rt, response = response, drift_index = drift_index, bound_index = bound_index, resid_index = resid_index, sv_index = sv_index, sw_index = sw_index, st0_index = st0_index) %>%
        group_by(drift_index, bound_index, resid_index, sv_index, sw_index, st0_index, response) %>%
        reframe(rt_q = quantile(rt, probs = rt_p)) %>%
        mutate(rt_p = rep(rt_p, n() / length(rt_p))) %>%
        complete(nesting(drift_index, bound_index, resid_index, sv_index, sw_index, st0_index), response, rt_p, fill = list(rt_q = NA))
    
    obs_p_resp <- tibble(rt = rt, response = response, drift_index = drift_index, bound_index = bound_index, resid_index = resid_index, sv_index = sv_index, sw_index = sw_index, st0_index = st0_index) %>%
        group_by(drift_index, bound_index, resid_index, sv_index, sw_index, st0_index, response) %>%
        summarize(n_resp = n(), .groups = "keep") %>%
        ungroup() %>%
        complete(nesting(drift_index, bound_index, resid_index, sv_index, sw_index, st0_index), response, fill = list(n_resp = 0)) %>%
        group_by(drift_index, bound_index, resid_index, sv_index, sw_index, st0_index) %>%
        mutate(p_resp = n_resp / sum(n_resp))
    
    if (!is.null(par)) {
        par_names <- c(paste0("a[", 1:n_bound, "]"), paste0("v[", 1:n_drift, "]"), paste0("w[", 1:n_bound, "]"), paste0("t0[", 1:n_resid, "]"), paste0("sv[", 1:n_sv, "]"), paste0("sw[", 1:n_sw, "]"), paste0("st0[", 1:n_st0, "]"))
        
        par_to_use <- rep(0, length(par_names))
        names(par_to_use) <- par_names
        overlap <- intersect(names(par_to_use), names(par))
        par_to_use[overlap] <- par[overlap]
        
        fitDF <- expand_grid(nesting(drift_index, bound_index, resid_index, sv_index, sw_index, st0_index), response = c("upper", "lower"), rt_p = rt_p) %>% mutate(rt_q = NA, p_resp = NA)
        
        for (i in 1:nrow(fitDF)) {
            fitDF$rt_q[i] <- qWDM(p = fitDF$rt_p[i], response = fitDF$response[i], a = par_to_use[paste0("a[", fitDF$bound_index[i], "]")], v = par_to_use[paste0("v[", fitDF$drift_index[i], "]")], w = par_to_use[paste0("w[", fitDF$bound_index[i], "]")], t0 = par_to_use[paste0("t0[", fitDF$resid_index[i], "]")], sv = par_to_use[paste0("sv[", fitDF$sv_index[i], "]")], sw = par_to_use[paste0("sw[", fitDF$sw_index[i], "]")], st0 = par_to_use[paste0("st0[", fitDF$st0_index[i], "]")])
            
            fitDF$p_resp[i] <- WienerCDF(t = Inf, response = fitDF$response[i], a = par_to_use[paste0("a[", fitDF$bound_index[i], "]")], v = par_to_use[paste0("v[", fitDF$drift_index[i], "]")], w = par_to_use[paste0("w[", fitDF$bound_index[i], "]")], t0 = par_to_use[paste0("t0[", fitDF$resid_index[i], "]")], sv = par_to_use[paste0("sv[", fitDF$sv_index[i], "]")], sw = par_to_use[paste0("sw[", fitDF$sw_index[i], "]")], st0 = par_to_use[paste0("st0[", fitDF$st0_index[i], "]")])$value
        }
        
        if (is.numeric(response)) {
            fitDF <- fitDF %>%
                mutate(response = as.numeric(factor(response, levels = c("lower", "upper"))))
        }
        
        obs_fit_data <- full_join(
            full_join(obs_p_resp, obs_rt_quantiles) %>% mutate(source = "Observed"),
            fitDF %>% mutate(source = "Fitted")
        )
        
        return(obs_fit_data)
    } else {
        return(full_join(obs_p_resp, obs_rt_quantiles))
    }
}