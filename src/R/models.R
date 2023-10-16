
# Core CAG Modeling classes

CAGBase <- R6Class("CAGBase", list(
    state_fields = NULL,

    export_state = function(field_names=self$state_fields) {
        class_name = class(self)[[1]]
        state_names = c("$CLASS$")
        state_list = list(class_name)
        for (name in field_names) {
            value = CAGBase$encode_value(get(name, self))
            state_names = c(state_names, name)
            state_list = append(state_list, list(value))
        }
        names(state_list) = state_names
        return(state_list)
    },

    import_state = function(data) {
        for (name in names(data)) {
            if (name != "$CLASS$") {
                self[[name]] = CAGBase$decode_value(data[[name]], name)
            }
        }
        return(invisible(self))
    }
))

CAGBase$encode_value = function(value) {
    if (is.environment(value) && is.function(value$export_state)) {
        return(value$export_state())
    } else if (is.list(value)) {
        return(lapply(value, CAGBase$encode_value))
    } else {
        return(value)
    }
}

CAGBase$is_object_state = function(value) {
    return(is.list(value) && !is.null(names(value)) && names(value)[[1]] == "$CLASS$")
}

debug.decode.value = FALSE
debug.decode.value.indent = 0
CAGBase$decode_value = function(value, tag=NULL) {
    if (debug.decode.value) {
        spacer = paste0(rep(".",debug.decode.value.indent),collapse="")
        context = if (is.null(tag)) "" else tag
        debug.decode.value.indent <<- debug.decode.value.indent + 1
        cat(sprintf("#DBG: %s In decode_value %s ...\n", spacer, context))
    }
    if (CAGBase$is_object_state(value)) {
        if (debug.decode.value) {
            cat(sprintf("#DBG: %s In decode_value %s [obj] ...\n", spacer, context))
        }
        classref = get(value[[1]])
        #cat(sprintf("#DBG: %s About to instantiate ...\n", spacer))
        if (!is.null(classref$instantiate)) {
            result = classref$instantiate()
        } else {
            result = classref$new()
        }
        result$import_state(value)
        if (debug.decode.value) {
            cat(sprintf("#DBG: %s Out decode_value ...\n", spacer))
            cat(sprintf("#DBG: %s object:\n", spacer))
            print(result)
            debug.decode.value.indent <<- debug.decode.value.indent - 1
        }
        return(result)
    }
    if (is.list(value)) {
        if (debug.decode.value) {
            cat(sprintf("#DBG: %s In decode_value %s [list] ...\n", spacer, context))
        }
        result = lapply(value, CAGBase$decode_value)
        if (debug.decode.value) {
            cat(sprintf("#DBG: %s Out decode_value ...\n", spacer))
            debug.decode.value.indent <<- debug.decode.value.indent - 1
        }
        return(result)
    }
    if (debug.decode.value) {
        cat(sprintf("#DBG: %s In decode_value %s [default] ...\n", spacer, context))
        cat(sprintf("#DBG: %s Out decode_value ...\n", spacer))
        debug.decode.value.indent <<- debug.decode.value.indent - 1
    }
    return(value)
}

CAGBase$load = function(inFile) {
    data = readRDS(inFile)
    #cat(sprintf("#DBG: read data:\n"))
    #print(data)
    return(CAGBase$decode_value(data, "root"))
}

CAGDonor <- R6Class("CAGDonor", inherit=CAGBase, public=list(
    id = NULL,
    age = NULL,
    inherited_cag = NULL,

    initialize = function(id, age=NULL, inherited_cag=NULL) {
        self$id = id
        self$age = age
        self$inherited_cag = inherited_cag
        self$state_fields = c("id", "age", "inherited_cag")
    }
))

CAGDonor$instantiate = function() {
    CAGDonor$new(NULL)
}

CAGDataSet <- R6Class("CAGDataSet", inherit=CAGBase, public=list(
    name = NULL,
    donor = NULL,
    region = NULL,
    celltype = NULL,
    observed = NULL,
    loss_fraction = NULL,

    initialize = function(donor, region=NULL, celltype=NULL, observed=NULL, name=NULL) {
        self$name = name
        self$donor = donor
        self$region = region
        self$celltype = celltype
        self$observed = observed
        self$state_fields = c("name", "donor", "region", "celltype", "observed", "loss_fraction")
    }
))

CAGDataSet$instantiate = function() {
    CAGDataSet$new(NULL)
}

optim.iteration = 0
optim.progress.interval = 10

CAGExpansionModel <- R6Class("CAGExpansionModel", inherit=CAGBase, public=list(
    name = NULL,
    data = NULL,
    fit_score = NULL,
    fitted_parameters = NULL,
    fixed_parameters = NULL,
    start_parameters_list = NULL,
    parameter_bounds = NULL,
    parameter_steps = NULL,
    debug = FALSE,
    verbose = TRUE,
    single_start = FALSE,

    initialize = function(name,
                          data=NULL,
                          fixed_parameters=NULL,
                          start_parameters_list=NULL,
                          parameter_bounds=NULL,
                          parameter_steps=NULL) {
       self$name = name
       self$data = data
       self$fixed_parameters = fixed_parameters
       self$start_parameters_list = start_parameters_list
       self$parameter_bounds = parameter_bounds
       self$parameter_steps = parameter_steps
       self$state_fields = c("name", "data", "fit_score", "fitted_parameters", "fixed_parameters", "start_parameters_list", "parameter_bounds", "parameter_steps")
    },

    fit = function(data=NULL, method="BFGS") {
        # TBD: Dispatch here?
        if (method == "BFGS") {
            self$fit_BFGS(data)
        } else {
            warn(sprintf("Unrecognized model fit method: %s", method))
        }
        return(invisible(self))
    },

    fit_BFGS = function(data=NULL, starts=self$start_parameters_list) {
        if (!is.null(data)) {
            self$data = data
        }
        best_score = -Inf
        best_params = NULL
        n_starts = length(starts)
        if (self$single_start) {
            n_starts = 1
        }
        for (i in 1:n_starts) {
            optim.iteration <<- 0
            opt = optim(par = self$scale_to_bounds(starts[[i]]),
                        fn = self$run_optim_fit_model,
                        method = "L-BFGS-B",
                        lower = 0,
                        upper = 1)
            score = -opt$value
            if (score > best_score) {
                best_score = score
                best_params = self$revert_scale(opt$par)
            }
            if (self$verbose) {
                params = self$revert_scale(opt$par)
                if (length(self$data) == 1) {
                    loss = if (is.null(self$data[[1]]$loss_fraction)) { "NA" } else { sprintf("%1.3f", self$data[[1]]$loss_fraction) }
                    cat(sprintf("#INFO %s: BFGS %s (age=%s, inh=%s, reg=%s, ct=%s, loss=%s, run=%d): %s -> %s\n", date(), 
                                self$data[[1]]$donor$id, self$data[[1]]$donor$age, self$data[[1]]$donor$inherited_cag,
                                self$data[[1]]$region, self$data[[1]]$celltype, loss, i,
                                paste0(params, collapse=" "), score))
                } else {
                    cat(sprintf("#INFO %s: BFGS FitToAll[n=%d] (run=%d): %s -> %s\n", date(), length(self$data),
                                paste0(params, collapse=" "), score))
                }
            }
        }
        self$fit_score = best_score
        if (!is.null(best_params)) {
            names(best_params) = names(self$parameter_bounds)
        }
        self$fitted_parameters = best_params
        return(invisible(self))
    },

    run_fit_model = function(params) {
        score = 0
        for (i in 1:length(self$data)) {
            age = self$data[[i]]$donor$age
            inherited = self$data[[i]]$donor$inherited_cag
            loss_factor = self$compute_loss_factor(self$data[[i]])
            #cat(sprintf("#DBG: generate_Qmatrix\n"))
            Qmatrix = self$generate_Qmatrix(params)
            #cat(sprintf("#DBG: generate_Qmatrix done\n"))
            #cat(sprintf("#DBG: calling expm\n"))
            cag_dist = expm::expm(Qmatrix * age)[inherited,]
            #cat(sprintf("#DBG: calling score_model\n"))
            data_set_score = self$score_model(self$data[[i]]$observed, loss_factor, cag_dist)
            #cat(sprintf("#DBG: adding score\n"))
            score = score + data_set_score
            if (self$debug) {
                cat(sprintf("#DBG: %s (age=%s, inh=%s, reg=%s, ct=%s, lossfact=%s): %s -> %s\n",
                            self$data[[i]]$donor$id, age, inherited,
                            self$data[[i]]$region, self$data[[i]]$celltype, loss_factor,
                            paste0(params, collapse=" "), data_set_score))
            }
        }
        return(score)
    },

    run_optim_fit_model = function(params) {
        # This wrapper is designed to be called by optim()
        # The params to are scaled and since optim() runs minimization by default, negate the returned score.
        optim.iteration <<- optim.iteration + 1
        if ((optim.iteration %% optim.progress.interval) == 0) {
            cat(sprintf("#DBG: %s Running optim iteration %d ...\n", date(), optim.iteration))
        }
        return(-self$run_fit_model(self$revert_scale(params)))
    },

    generate_Qmatrix = function(params, maxcag=NULL) {
        # This method needs to be overridden by subclasses
        stop(sprintf("Model %s failed to override generate_Qmatrix", self$name))
    },

    score_model = function(cag_observed, loss_factor, model_cag_dist) {
        log_likelihood = 0
        for (val in cag_observed) {
            log_likelihood = log_likelihood + log10(model_cag_dist[min(val, length(model_cag_dist))])
        }
        # Include cell loss factor adjustment
        log_likelihood = log_likelihood + loss_factor * length(cag_observed) * log10(model_cag_dist[length(model_cag_dist)])
        if (is.nan(log_likelihood) | is.infinite(log_likelihood)) {
            log_likelihood = -1e9
        }
        return(log_likelihood)
    },

    compute_loss_factor = function(data) {
        loss_factor = 0
        loss_fraction = self$get_loss_fraction(data)
        if (loss_fraction > 0) {
            # Cap for computing the loss fraction, by default limits the loss factor to at most 10.0
            if (!is.null(self$fixed_parameters$loss_cap)) {
                loss_cap = self$fixed_parameters$loss_cap
            } else {
                # Default loss cap and cache in fixed_parameters for posterity
                loss_cap = 0.90
                self$fixed_parameters = append(self$fixed_parameters, list(loss_cap))
            }
            adj_loss = min(loss_fraction, loss_cap)
            loss_factor = adj_loss/(1-adj_loss)
        }
        return(loss_factor)
    },

    get_loss_fraction = function(data) {
        if (is.null(data$loss_fraction) || is.na(data$loss_fraction) || length(data$loss_fraction) == 0) {
            return(0)
        }
        return(data$loss_fraction)
    },

    # Scales parameter values to be between 0 and 1 given bounds
    scale_to_bounds = function(params, bounds=self$parameter_bounds) {
        # TBD: Do this by name? Will it be too slow?
        #cat(sprintf("#DBG: scale_to_bounds:\n"))
        #print(params)
        #print(bounds)
        vapply(1:length(params),
               function(i) {
                   (params[[i]] - bounds[[i]][[1]]) / (bounds[[i]][[2]] - bounds[[i]][[1]]) },
               numeric(1))
    },

    save = function(outFile) {
        data = self$export_state()
        saveRDS(data, outFile)
        return(invisible(self))
    },

    # Reverses the scaling transformation of scale_to_bounds()
    revert_scale = function(params, bounds=self$parameter_bounds) {
        # TBD: Do this by name? Will it be too slow?
        vapply(1:length(params),
               function(i) {
                   params[[i]] * (bounds[[i]][[2]] - bounds[[i]][[1]]) + bounds[[i]][[1]] },
               numeric(1))
    },

    plot_fit = function(title=NULL, id=FALSE) {
        cat(sprintf("%s: plot_fit ...\n", date()))
        plots = lapply(self$data, self$plot_fit1, id=id)
        if (length(self$data) == 1) {
            result = plots[[1]]
        } else {
            result = plot_grid(plots, ncol=1)
        }
        return(result)
    },

    plot_fit1 = function(data, id=FALSE, plot.last.bin=FALSE) {
        pdf_plot = self$plot_fit_pdf(data, id, plot.last.bin=plot.last.bin)
        cdf_plot = self$plot_fit_cdf(data, "", plot.last.bin=plot.last.bin)
        qq_plot = self$plot_fit_qq(data, "")
        loss_plot = self$plot_fit_loss_bars(data, id="")
        plot_grid(pdf_plot, cdf_plot, qq_plot, loss_plot,
                  ncol=4, rel_widths=c(2,2,2,1.5))
    },

    plot_fit_pdf = function(data, id=FALSE, plot.last.bin=FALSE) {
        xlim1 = 1
        if (plot.last.bin) {
            xlim2 = self$fixed_parameters$max_cag
        } else {
            xlim2 = self$fixed_parameters$max_cag - 1
        }
        xseq = seq(xlim1, xlim2)
        cags = pmin(data$observed, rep(self$fixed_parameters$max_cag, length(data$observed)))
        data_hist = sapply(xseq, function(x) { sum(cags == x) })
        model_hist = self$get_model_dist()[xseq]
        loss_factor = self$compute_loss_factor(data)
        model_hist = model_hist * (1+loss_factor)
        plot_data = data.frame(CAG=xseq, data=data_hist, model=model_hist)
        plot =
            ggplot(plot_data) +
            geom_col(aes(x=CAG, y=model), color="orange", fill="orange", show.legend=FALSE) +
            geom_col(aes(x=CAG, y=data), color=NA, fill="black", width=0.5, show.legend=FALSE) +
            labs(x="CAG length", y="Cells")
        if (isTRUE(id)) {
            plot = plot + ggtitle(data$donor$id)
        } else if (typeof(id) == "character") {
            plot = plot + ggtitle(id)
        }
        return(plot)
    },

    plot_fit_cdf = function(data, id=FALSE, plot.last.bin=FALSE) {
        xlim1 = 1
        if (plot.last.bin) {
            xlim2 = self$fixed_parameters$max_cag
        } else {
            xlim2 = self$fixed_parameters$max_cag - 1
        }
        xseq = seq(xlim1, xlim2)
        cags = pmin(data$observed, rep(self$fixed_parameters$max_cag, length(data$observed)))
        data_hist = sapply(xseq, function(x) { sum(cags == x) })
        model_hist = self$get_model_dist()[xseq]
        data_cdf = cumsum(data_hist)/sum(data_hist)
        model_cdf = cumsum(model_hist)/sum(model_hist)
        plot_data = data.frame(CAG=xseq, data=data_cdf, model=model_cdf)
        plot = 
            ggplot(plot_data) +
            geom_line(aes(x=CAG, y=model), color="orange", show.legend=FALSE) +
            geom_line(aes(x=CAG, y=data), color="black", show.legend=FALSE) +
            labs(x="CAG length", y="CDF")
        if (isTRUE(id)) {
            plot = plot + ggtitle(data$donor$id)
        } else if (typeof(id) == "character") {
            plot = plot + ggtitle(id)
        }
        return(plot)
    },

    plot_fit_qq = function(data, id=FALSE) {
        age = data$donor$age
        inherited = data$donor$inherited_cag
        cags = data$observed
        ##maxcag = max(cags)
        maxcag = self$fixed_parameters$max_cag
        probs = tail(seq(0,maxcag)/maxcag, -1)
        data_quantiles = quantile(data$observed[data$observed <= maxcag], probs)

        Qmatrix = self$generate_Qmatrix(self$fitted_parameters, maxcag+1)
        model_dist = expm::expm(Qmatrix * age)[inherited,]
        model_dist = model_dist[seq(1,maxcag)]
        model_dist = model_dist/sum(model_dist)
        model_cdf = cumsum(model_dist)
        model_probs = (probs - 0.5/maxcag)
        model_quantiles = sapply(model_probs, function(p) { match(TRUE, model_cdf >= p) })


        plot_data = data.table(probs=probs, data=data_quantiles, model=model_quantiles)
        #cat(sprintf("#DBG:\n"))
        #write.table(plot_data)
        minplot = min(plot_data$data, plot_data$model)
        maxplot = max(plot_data$data, plot_data$model)
        plot_width = maxplot-minplot

        kstest = suppressWarnings(ks.test(plot_data$data, plot_data$model))
        kstext = sprintf("KS p=%1.3f", kstest$p.value)
        #cat(sprintf("#DBG: kstest:\n"))
        #print(kstext)

        plot =
            ggplot(plot_data) +
            geom_point(aes(x=model, y=data)) +
            geom_abline(slope=1, color="gray") +
            geom_text(data=data.frame(label=kstext, x=minplot + 0.02*plot_width, y=maxplot - 0.02*plot_width),
                      aes(x=x, y=y, label=label), size=3.0, hjust=0, vjust=1) +
            xlim(c(minplot,maxplot)) + ylim(c(minplot,maxplot)) +
            labs(x="CAG distribution (model)", y="CAG distribution (observed)")
        if (isTRUE(id)) {
            plot = plot + ggtitle(data$donor$id)
        } else if (typeof(id) == "character") {
            plot = plot + ggtitle(id)
        }
        return(plot)
    },

    plot_fit_loss_bars = function(data, id=FALSE) {
        maxcag = self$fixed_parameters$max_cag
        data_value = sum(data$observed >= maxcag)/length(data$observed)
        loss_factor = self$compute_loss_factor(data)
        # Adjust data value for loss estimate
        data_value = (data_value + loss_factor) / (1 + loss_factor)
        data_value = min(data_value, 1.0)
        model_value = self$get_model_pdf(data$donor)[maxcag]
        categories = c("Estimated", "Model")
        plot_data = data.table(type=factor(categories, levels=categories), value=c(data_value, model_value))
        plot = 
            ggplot(plot_data) +
            geom_col(aes(x=type, y=value, fill=type), width=0.6, show.legend=FALSE) +
            scale_fill_manual(values=c(Estimated="darkgray", Model="orange")) +
            ylim(c(0,1)) +
            labs(y="Cell loss", x="")
        if (isTRUE(id)) {
            plot = plot + ggtitle(data$donor$id)
        } else if (typeof(id) == "character") {
            plot = plot + ggtitle(id)
        }
        return(plot)
    },

    # private
    get_model_dist = function(donor=NULL) {
        data = self$get_data_for_donor(donor)
        if (is.null(data)) {
            return(NULL)
        }
        n = max(c(1,length(data$observed)))
        return(self$get_model_pdf(data$donor) * n)
    },

    # private
    get_model_pdf = function(donor, area.scale=1) {
        age = donor$age
        inherited = donor$inherited_cag
        Qmatrix = self$generate_Qmatrix(self$fitted_parameters)
        model_dist = expm::expm(Qmatrix * age)[inherited,]
        return(area.scale * model_dist)
    },

    get_model_rate_curve = function() {
        Qmatrix = self$generate_Qmatrix(self$fitted_parameters)
        rates = sapply(1:nrow(Qmatrix), function(ridx) { sum((1:nrow(Qmatrix)-ridx)*Qmatrix[ridx,]) })
        # Do not return the sync state
        return(head(rates,-1))
    },

    get_model_CAG_curves = function(donor=NULL, probs=c(0.05,0.5,0.95), times=seq(0,100,5), maxcag=NULL) {
        data = self$get_data_for_donor(donor)
        if (is.null(data)) {
            return(NULL)
        }
        inherited = data$donor$inherited_cag
        Qmatrix = self$generate_Qmatrix(self$fitted_parameters, maxcag)
        cags = sapply(times, function(t) { if (t == 0) { return(rep(inherited, length(probs))) }
                                           pdf = expm::expm(Qmatrix*t)[inherited,];
                                           af = approxfun(c(0,cumsum(pdf)),0:length(pdf),ties="ordered");
                                           return(af(probs)) })
        if (length(probs) == 1) {
            cags = matrix(cags, nrow=1)
        }
        colnames(cags) = times
        rownames(cags) = probs
#        cags = sapply(times, function(t) { pdf = expm::expm(Qmatrix*t)[inherited,];
#                                           #cat(sprintf("#DBG: pdf(%d):\n", t));
#                                           #print(pdf);
#                                           sum(seq(1:length(pdf))*pdf) })
        return(cags)
    },

    # private
    get_data_for_donor = function(donor=NULL) {
        if (is.null(donor)) {
            if (length(self$data) != 1) {
                stop("Must specify donor when model is fit to multiple data sets")
            } else {
                return(self$data[[1]])
            }
        } else {
            donors = lapply(self$data, function(d) { d$donor })
            inds = which(donors == donor)
            if (length(inds) == 0) {
                return(NULL)
            } else if (length(inds) > 0) {
                stop(sprintf("Multiple models for donor %s", donor$id))
            } else {
                return(self$data[[inds]])
            }
        }
    },

    compute_confidence_intervals = function(interval=0.95) {
        result = lapply(names(self$fitted_parameters), self$compute_confidence_interval1, interval)
        names(result) = names(self$fitted_parameters)
        return(result)
    },

    compute_confidence_interval1 = function(param, interval) {
        value = self$fitted_parameters[[param]]
        if (is.null(value)) {
            return(NULL)
        }
        logThreshold = -log10(1 - interval)
        low = self$compute_confidence_interval_search(param, -1, logThreshold)
        high = self$compute_confidence_interval_search(param, 1, logThreshold)
        return(c(low, high))
    },

    compute_confidence_interval_search = function(param, sign, logThreshold) {
        # Default method, binary search based on bounds and fitted value
        fit = self$fit_score
        value = self$fitted_parameters[[param]]
        bounds = self$parameter_bounds[[param]]
        debug.count = 0
        if (sign > 0) {
            bottom = value
            top = bounds[[2]]
        } else {
            bottom = bounds[[1]]
            top = value
        }
        value.epsilon = abs(value*0.0001)
        score.epsilon = 0.001
        target = fit - logThreshold
        while (TRUE) {
            mid = (bottom+top) / 2
            test_params = self$fitted_parameters
            test_params[[param]] = mid
            debug.count = debug.count + 1
            score = self$run_fit_model(test_params)
            #cat(sprintf("#DBG: %s confint loop %d %s %s/%s/%s -> %s\n", date(), debug.count, param, bottom, mid, top, target - score))
            if (abs(target-score) < score.epsilon || abs(top-bottom) < value.epsilon) {
                break
            } else if ((sign > 0 && score >= target) || (sign < 0 && score <= target)) {
                bottom = mid
            } else {
                top = mid
            }
        }
        cat(sprintf("#DBG: %s confint search iterations for %s: %d\n", date(), param, debug.count))
        return(mid)
    },

    compute_confidence_interval_search_v1 = function(param, sign, logThreshold) {
        # Older version, linear search based on step size followed by binary search
        # Does not use as fancy epsilon calculation
        # Saved for posterity in case we find we need this due to non-smooth distributions
        fit = self$fit_score
        value = self$fitted_parameters[[param]]
        step = self$parameter_steps[[param]]
        value.curr = value
        debug.count = 0
        while (TRUE) {
            value.next = value.curr + sign * step
            test_params = self$fitted_parameters
            test_params[[param]] = value.next
            debug.count = debug.count + 1
            score = self$run_fit_model(test_params)
            cat(sprintf("#DBG: %s confint loop 1/%d %s %s %s -> %s\n", date(), debug.count, param, value.curr, value.next, fit - logThreshold - score))
            if (score <= fit - logThreshold) {
                break
            }
            value.curr = value.next
        }
        if (sign > 0) {
            bottom = value.curr
            top = value.next
        } else {
            bottom = value.next
            top = value.curr
        }
        epsilon = step * 0.01
        while (top-bottom > epsilon) {
            mid = (bottom+top) / 2
            test_params = self$fitted_parameters
            test_params[[param]] = mid
            debug.count = debug.count + 1
            score = self$run_fit_model(test_params)
            cat(sprintf("#DBG: %s confint loop 2/%d %s %s -> %s\n", date(), debug.count, param, mid, fit - logThreshold - score))
            if (score <= fit - logThreshold) {
                bottom = mid
            } else {
                top = mid
            }
        }
        cat(sprintf("#DBG: %s confint search iterations for %s: %d\n", date(), param, debug.count))
        return((bottom + top)/2)
    },

    grid_search_internal = function(fits, params_to_search, param_bounds, param_steps, param_values) {
        ranges = NULL
        for (p in params_to_search) {
            low = param_bounds[[p]][[1]]
            high = param_bounds[[p]][[2]]
            step = param_steps[[p]]
            ranges = append(ranges, list(seq(low, high, step)))
        }
        names(ranges) = names(params_to_search)
        search_matrix = do.call(expand.grid, ranges)
        for (i in 1:nrow(search_matrix)) {
            test_params = append(as.list(c(search_matrix[i,])), param_values)
            fit = find_fit(fits, test_params)
            if (is.null(fit)) {
                fit = fit_model(test_params)
                fits = rbind(fits, fit)
            }
        }
        return(fits)
    },

    generate_Q_matrix_internal = function(pexp, rates, max_cag, contraction_limit=NULL) {
        Q = matrix(0, nrow=max_cag, ncol=max_cag)
        for (i in 2:(max_cag-1)) {
            if (is.vector(rates)) {
                rate = rates[[i]]
            } else if (is.function(rates)) {
                rate = rates(i)
            }
            expansion_rate = rate * pexp
            contraction_rate = rate * (1-pexp)
            if (!is.null(contraction_limit) && (i - 1 < contraction_limit) && rate > 0) {
                ## Random experiment (double expansion rate when limiting contraction)
                ##cat(sprintf("#DBG: override expansion rate for Q[%d] -> %1.3f\n", i, expansion_rate + contraction_rate))
                ##expansion_rate = expansion_rate + contraction_rate
                #cat(sprintf("#DBG: override contraction rate for Q[%d] -> 0\n", i))
                contraction_rate = 0
            }
            Q[i, i-1] = contraction_rate
            Q[i, i+1] = expansion_rate
            Q[i, i] = -(expansion_rate + contraction_rate)
        }
        return(Q)
    }
))

CAGExpansionModel$load = function(inFile) {
    CAGBase$load(inFile)
}

TwoPhaseModel = R6Class("TwoPhaseModel", inherit=CAGExpansionModel, public=list(
    initialize = function(data=NULL) {
        super$initialize(name="TwoPhase",
                         data=data,
                         fixed_parameters=
                             list(max_cag = 180, threshold1 = 33.5),
                         start_parameters_list=
                             list(list(rate1 = 0.1, rate2 = 1, pexp = 0.7, threshold2 = 65),
                                  list(rate1 = 0.125, rate2 = 1.5, pexp = 0.7, threshold2 = 70),
                                  list(rate1 = 0.15, rate2 = 2, pexp = 0.7, threshold2 = 75)),
                         parameter_bounds=
                             list(rate1 = c(0, 0.3),
                                  rate2 = c(0, 3),
                                  pexp =  c(0.5, 1),
                                  threshold2 = c(40, 100)),
                         parameter_steps =
                             list(rate1 = 0.001,
                                  rate2 = 0.01,
                                  pexp = 0.001,
                                  threshold2 = 0.1))
    },
    generate_Qmatrix = function(params, maxcag=NULL) {
        # Note: params here are by position and are not named
        r1 = params[[1]]
        r2 = params[[2]]
        pexp = params[[3]]
        thresh1 = self$fixed_parameters$threshold1
        thresh2 = params[[4]]
        rates = function(x) {
            max(r1 * (x - thresh1), 0) + max(r2 * (x - thresh2), 0)
        }
        if (is.null(maxcag)) {
            maxcag = self$fixed_parameters$max_cag
        }
        return(self$generate_Q_matrix_internal(pexp, rates, maxcag, self$fixed_parameters$contraction_limit))
    }
))

TwoPhaseModelSmoothed = R6Class("TwoPhaseModelSmoothed", inherit=CAGExpansionModel, public=list(
    initialize = function(data=NULL) {
        super$initialize(name="TwoPhaseSmoothed",
                         data=data,
                         fixed_parameters=
                             list(max_cag = 180, threshold1 = 33.5),
                         start_parameters_list=
                             list(list(rate1 = 0.1, rate2 = 1, pexp = 0.7, threshold2 = 65, smooth.radius=3),
                                  list(rate1 = 0.125, rate2 = 1.5, pexp = 0.7, threshold2 = 70, smooth.radius=3),
                                  list(rate1 = 0.15, rate2 = 2, pexp = 0.7, threshold2 = 75, smooth.radius=3)),
                         parameter_bounds=
                             list(rate1 = c(0.05, 0.3),
                                  rate2 = c(0.05, 3),
                                  pexp =  c(0.5, 1),
                                  threshold2 = c(40, 100),
                                  smooth.radius = c(0, 10)),
                         parameter_steps =
                             list(rate1 = 0.001,
                                  rate2 = 0.01,
                                  pexp = 0.001,
                                  threshold2 = 0.1,
                                  smooth.radius = 0.1))
    },
    generate_Qmatrix = function(params) {
        # Note: params here are by position and are not named
        r1 = params[[1]]
        r2 = params[[2]]
        pexp = params[[3]]
        thresh1 = self$fixed_parameters$threshold1
        thresh2 = params[[4]]
        smooth.radius = params[[5]]
        rates = function(x) {
            if (x < thresh1) {
                return(0)
            }
            if (x < thresh2) {
                return(r1 * (x-thresh1))
            }
            return(r1 * (x-thresh1) + smooth_between(x,0,r2*(x-thresh2),thresh2+smooth.radius/2,smooth.radius))

            ## This didn't really pan out ...
            ##return(smooth_between(x,r1,r2,thresh2,smooth.radius)*(x-thresh1))
        }
        return(self$generate_Q_matrix_internal(pexp, rates, self$fixed_parameters$max_cag, self$fixed_parameters$contraction_limit))
    }
))

TwoPhaseModelExp1 = R6Class("TwoPhaseModelExp1", inherit=CAGExpansionModel, public=list(
    initialize = function(data=NULL) {
        super$initialize(name="TwoPhaseExp1",
                         data=data,
                         fixed_parameters=
                             list(max_cag = 180, threshold1 = 33.5),
                         start_parameters_list=
                             list(list(rate1 = 0.1, rate2 = 1, pexp = 0.7, threshold2 = 65, exp2=1)),
                         parameter_bounds=
                             list(rate1 = c(0.01, 0.3),
                                  rate2 = c(0.01, 1),
                                  pexp =  c(0.5, 1),
                                  threshold2 = c(40, 100),
                                  exp2 = c(0.5,3)),
                         parameter_steps =
                             list(rate1 = 0.001,
                                  rate2 = 0.01,
                                  pexp = 0.001,
                                  threshold2 = 0.1,
                                  exp2 = 0.1))
    },
    generate_Qmatrix = function(params) {
        # Linear to thresh2, then power law to see what value of the exponent is suggested:
        # r1 * (x-thresh1) + r2 * (x-thresh2)^a
 
        # Note: params here are by position and are not named
        r1 = params[[1]]
        r2 = params[[2]]
        pexp = params[[3]]
        thresh1 = self$fixed_parameters$threshold1
        thresh2 = params[[4]]
        exp2 = params[[5]]
        rates = function(x) {
            if (x < thresh1) {
                return(0)
            }
            if (x < thresh2) {
                return(r1 * (x-thresh1))
            }
            return(r1 * (x-thresh1) + r2 * (x-thresh2)^exp2)
        }
        return(self$generate_Q_matrix_internal(pexp, rates, self$fixed_parameters$max_cag))
    }
))

TwoPhaseModelExp2 = R6Class("TwoPhaseModelExp2", inherit=CAGExpansionModel, public=list(
    initialize = function(data=NULL) {
        super$initialize(name="TwoPhaseExp2",
                         data=data,
                         fixed_parameters=
                             list(max_cag = 180, threshold1 = 33.5),
                         start_parameters_list=
                             list(list(rate1 = 0.1, rate2 = 1, pexp = 0.7, threshold2 = 65, exp1=1, exp2=1)),
                         parameter_bounds=
                             list(rate1 = c(0.001, 0.3),
                                  rate2 = c(0.001, 3),
                                  pexp =  c(0.5, 1),
                                  threshold2 = c(40, 100),
                                  exp1 = c(0.5,3),
                                  exp2 = c(0.1,3)),
                         parameter_steps =
                             list(rate1 = 0.001,
                                  rate2 = 0.01,
                                  pexp = 0.001,
                                  threshold2 = 0.1,
                                  exp1 = 0.1,
                                  exp2 = 0.1))
    },
    generate_Qmatrix = function(params, maxcag=NULL) {
        # Arbitrary power in both phases:
        # r1 * (x-thresh1)^a1 + r2 * (x-thresh2)^a2
 
        # Note: params here are by position and are not named
        r1 = params[[1]]
        r2 = params[[2]]
        pexp = params[[3]]
        thresh1 = self$fixed_parameters$threshold1
        thresh2 = params[[4]]
        exp1 = params[[5]]
        exp2 = params[[6]]
        rates = function(x) {
            if (x < thresh1) {
                return(0)
            }
            if (x < thresh2) {
                return(r1 * (x-thresh1)^exp1)
            }
            return(r1 * (x-thresh1)^exp1 + r2 * (x-thresh2)^exp2)
        }
        if (is.null(maxcag)) {
            maxcag = self$fixed_parameters$max_cag
        }
        return(self$generate_Q_matrix_internal(pexp, rates, maxcag))
    }
))

smooth_between <- function(x, v1, v2, mid, radius) {
    if (x <= mid - radius) {
        return(v1)
    } else if (x >= mid + radius) {
        return(v2)
    } else if (radius == 0 || x == mid) {
        # Prevent overflow
        return(0.5 * (v1+v2))
    } else {
        t1 = exp(-2*radius/((x-mid)-radius))
        t2 = exp(2*radius/((x-mid)+radius))
        f = t1/(t1+t2)
        if (is.nan(f) || is.infinite(f)) {
            if (x < mid) {
                return(v1)
            } else {
                return(v2)
            }
        }
        result = v1 + f*v2;
        if (is.nan(result) || is.infinite(result)) {
            cat(sprintf("#DBG: smooth_between(%s, %s, %s, %s, %s) -> %s\n", x, v1, v2, mid, radius, result))
            cat(sprintf("#DBG: t1 = %s t2 = %s f = %s\n", t1, t2, f))
        }
        return(v1 + f*v2);
    }
}

TwoPhasePowerModel = R6Class("TwoPhasePowerModel", inherit=CAGExpansionModel, public=list(
    initialize = function(data=NULL) {
        super$initialize(name="TwoPhasePower",
                         data=data,
                         fixed_parameters=
                             list(max_cag = 180, threshold1 = 33.5),
                         start_parameters_list=
                             list(list(exp1 = 1.0, exp2 = 2.0, pexp = 0.7, threshold2 = 70),
                                  list(exp1 = 1.5, exp2 = 2.0, pexp = 0.65, threshold2 = 70)),
                         parameter_bounds=
                             list(exp1 = c(0.5,4),
                                  exp2 = c(0.5,4),
                                  pexp =  c(0.5, 1),
                                  threshold2 = c(40,100)),
                         parameter_steps =
                             list(exp1 = 0.01,
                                  exp2 = 0.01,
                                  pexp = 0.001,
                                  threshold2 = 0.1))
    },
    generate_Qmatrix = function(params) {
        # Note: params here are by position and are not named
        exp1 = params[[1]]
        exp2 = params[[2]]
        pexp = params[[3]]
        thresh1 = self$fixed_parameters$threshold1
        thresh2 = params[[4]]
        rates = function(x) {
            max(x - thresh1, 0)^exp1 + max(x - thresh2, 0)^exp2
        }
        return(self$generate_Q_matrix_internal(pexp, rates, self$fixed_parameters$max_cag))
    }
))

OnePhasePowerModel = R6Class("OnePhasePowerModel", inherit=CAGExpansionModel, public=list(
    initialize = function(data=NULL) {
        super$initialize(name="OnePhasePower",
                         data=data,
                         fixed_parameters=
                             list(max_cag = 180, threshold1 = 30),
                         start_parameters_list=
                             list(list(rate1 = 1, exp1 = 1, pexp = 0.7),
                                  list(rate1 = 1.5, exp1 = 1.5, pexp = 0.65)),
                         parameter_bounds=
                             list(rate1 = c(0.1,2),
                                  exp1 = c(1.0,4.0),
                                  pexp =  c(0.5, 1)),
                         parameter_steps =
                             list(rate1 = 0.01,
                                  exp1 = 0.01,
                                  pexp = 0.001))
    },
    generate_Qmatrix = function(params) {
        # Note: params here are by position and are not named
        rate1 = params[[1]]
        exp1 = params[[2]]
        pexp = params[[3]]
        thresh1 = self$fixed_parameters$threshold1
        rates = function(x) {
            rate1 * max(x - thresh1, 0)^exp1
        }
        return(self$generate_Q_matrix_internal(pexp, rates, self$fixed_parameters$max_cag))
    }
))

BioModel5a = R6Class("BioModel5a", inherit=CAGExpansionModel, public=list(
    initialize = function(data=NULL) {
        super$initialize(name="BioModel5a",
                         data=data,
                         fixed_parameters=
                             list(max_cag = 180, threshold1 = 33.5),
                         start_parameters_list=
                             list(list(rate1 = 0.003, rate2 = 0.03, pexp = 0.7, threshold2 = 65),
                                  list(rate1 = 0.005, rate2 = 0.05, pexp = 0.7, threshold2 = 70)),
                         parameter_bounds=
                             list(rate1 = c(0, 0.01),
                                  rate2 = c(0, 0.1),
                                  pexp =  c(0.5, 1),
                                  threshold2 = c(40, 100)),
                         parameter_steps =
                             list(rate1 = 0.001,
                                  rate2 = 0.01,
                                  pexp = 0.001,
                                  threshold2 = 0.1))
    },
    generate_Qmatrix = function(params) {
        # Note: params here are by position and are not named
        c1 = params[[1]]
        c2 = params[[2]]
        pexp = params[[3]]
        thresh1 = self$fixed_parameters$threshold1
        thresh2 = params[[4]]
        rates = function(x) {
            rate = 0
            if (x > thresh1) {
                rate = 0.5 * c1 * (min(x,thresh2) - thresh1)^2
            }
            if (x > thresh2) {
                rate = rate + 0.5 * c2 * (x-thresh2)^2 + c1 * (thresh2-thresh1) * (x - thresh2)
                ##rate = rate + 0.5 * c2 * (x-thresh2)^2 # creates bimodal distribution on AN09336
            }
            return(rate)
        }
        return(self$generate_Q_matrix_internal(pexp, rates, self$fixed_parameters$max_cag))
    }
))

BioModel5x = R6Class("BioModel5x", inherit=CAGExpansionModel, public=list(
    initialize = function(data=NULL) {
        super$initialize(name="BioModel5x",
                         data=data,
                         fixed_parameters=
                             list(max_cag = 180, threshold1 = 33.5),
                         start_parameters_list=
                             list(list(rate1 = 0.003, rate2 = 0.03, pexp = 0.7, threshold2 = 65),
                                  list(rate1 = 0.005, rate2 = 0.05, pexp = 0.7, threshold2 = 70)),
                         parameter_bounds=
                             list(rate1 = c(0, 0.1),
                                  rate2 = c(0, 0.2),
                                  pexp =  c(0.5, 1),
                                  threshold2 = c(40, 100)),
                         parameter_steps =
                             list(rate1 = 0.001,
                                  rate2 = 0.01,
                                  pexp = 0.001,
                                  threshold2 = 0.1))
    },
    generate_Qmatrix = function(params) {
        # Note: params here are by position and are not named
        r1 = params[[1]]
        r2 = params[[2]]
        pexp = params[[3]]
        thresh1 = self$fixed_parameters$threshold1
        thresh2 = params[[4]]
        rates = function(x) {
            max(r1 * (x - thresh1)^2, 0) + max(r2 * (x - thresh2)^2, 0)
        }
        return(self$generate_Q_matrix_internal(pexp, rates, self$fixed_parameters$max_cag))
    }
))

BioModel5b = R6Class("BioModel5b", inherit=CAGExpansionModel, public=list(
    initialize = function(data=NULL) {
        super$initialize(name="BioModel5b",
                         data=data,
                         fixed_parameters=
                             list(max_cag = 180, threshold1 = 33.5),
                         start_parameters_list=
                             list(list(rate1 = 0.5, rate2 = 1, pexp = 0.7, threshold2 = 65),
                                  list(rate1 = 1, rate2 = 3, pexp = 0.7, threshold2 = 70)),
                         parameter_bounds=
                             list(rate1 = c(0, 10),
                                  rate2 = c(0, 100),
                                  pexp =  c(0.5, 1),
                                  threshold2 = c(40, 100)),
                         parameter_steps =
                             list(rate1 = 0.001,
                                  rate2 = 0.01,
                                  pexp = 0.001,
                                  threshold2 = 0.1))
    },
    generate_Qmatrix = function(params) {
        # Note: params here are by position and are not named
        r1 = params[[1]]
        r2 = params[[2]]
        pexp = params[[3]]
        thresh1 = self$fixed_parameters$threshold1
        thresh2 = params[[4]]
        rates = function(x) {
            rate = 0
            if (x > thresh1) {
                rate = rate + r1*(0.5*(x-thresh1)^2 - thresh1*(x-thresh1) + 0.5*thresh1^2)
            }
            if (x > thresh2) {
                rate = rate + r2*(0.5*(x-thresh2)^2 - thresh2*(x-thresh2) + 0.5*thresh2^2)
            }
            return(rate)
        }
        return(self$generate_Q_matrix_internal(pexp, rates, self$fixed_parameters$max_cag))
    }
))

BioModel5c = R6Class("BioModel5c", inherit=CAGExpansionModel, public=list(
    initialize = function(data=NULL) {
        super$initialize(name="BioModel5c",
                         data=data,
                         fixed_parameters=
                             list(max_cag = 180, threshold1 = 33.5),
                         start_parameters_list=
                             list(list(rate1 = 0.5, rate2 = 1, pexp = 0.7, threshold2 = 65),
                                  list(rate1 = 1, rate2 = 3, pexp = 0.7, threshold2 = 70)),
                         parameter_bounds=
                             list(rate1 = c(0, 10),
                                  rate2 = c(0, 100),
                                  pexp =  c(0.5, 1),
                                  threshold2 = c(40, 100)),
                         parameter_steps =
                             list(rate1 = 0.001,
                                  rate2 = 0.01,
                                  pexp = 0.001,
                                  threshold2 = 0.1))
    },
    generate_Qmatrix = function(params) {
        # Note: params here are by position and are not named
        r1 = params[[1]]
        r2 = params[[2]]
        pexp = params[[3]]
        thresh1 = self$fixed_parameters$threshold1
        thresh2 = params[[4]]
        rates = function(x) {
            rate = 0
            if (x > thresh1) {
                rate = rate + r1*(min(x,thresh2) - thresh1 + thresh1*log(min(x,thresh2)/thresh1))
            }
            if (x > thresh2) {
                rate = rate + r2*(x - thresh2 + thresh2*log(x/thresh2))
            }
            return(rate)
        }
        return(self$generate_Q_matrix_internal(pexp, rates, self$fixed_parameters$max_cag))
    }
))

                        