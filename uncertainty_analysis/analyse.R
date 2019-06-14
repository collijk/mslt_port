#!/usr/bin/Rscript --vanilla

SCRIPT <- 'analyse.R'


file_pattern <- function(population, delay, tob_prev, interv) {
    delay_str <- paste0(delay, '-years')
    file_base <- paste('mslt_tobacco', population, delay_str, tob_prev,
                       interv, sep = '_')
    mslt_base <- paste0(file_base, '_mm_')
    mslt_pattern <- paste0(mslt_base, '.*\\.csv')
    return(mslt_pattern)
}


mean_file_pattern <- function(population, delay, tob_prev, interv) {
    delay_str <- paste0(delay, '-years')
    file_base <- paste('mslt_tobacco', population, delay_str, tob_prev,
                       interv, sep = '_')
    mslt_base <- paste0(file_base, '_mm')
    mslt_pattern <- paste0(mslt_base, '\\.csv')
    return(mslt_pattern)
}


summarise_files <- function(delay, tob_prev, interv, verbose = FALSE) {
    file_list <- list.files(
        pattern = file_pattern('.*', delay, tob_prev, interv))

    num_files <- length(file_list)

    if (num_files < 500) {
        return(NULL)
    }

    ## NOTE: for LYs and HALYs, report data for:
    ## 1. All cohorts combined (Maori and non-Maori).
    ## 2. Maori females.
    ## 3. Non-Maori males.

    ## For ACMR and YLDR, report data for:
    ## 1. Maori females aged 62 in 2041 and 2061; and
    ## 2. Non-Maori males aged 62 in 2041 and 2061.
    rate_years <- c(2041, 2061)
    rate_ages <- c(62)
    rate_cols <- c('year', 'age', 'sex')

    df_ly <- NULL
    df_acmr <- NULL
    df_yldr <- NULL

    for (ix in seq_along(file_list)) {
        file <- file_list[ix]
        if (verbose) {
            cat('Loading', file, '...')
        }
        is_non_maori <- grepl('_non-maori_', file)
        if (is_non_maori) {
            popn <- 'non-maori'
            rate_sex <- c('male')
        } else {
            popn <- 'maori'
            rate_sex <- c('female')
        }

        draw_number <- as.numeric(gsub('^.*_([0-9]+)\\.csv$', '\\1', file))
        ## TODO: load data for draw number 0 as the expected value.
        ## NOTE: this file will end in '_mm.csv'.

        df_in <- read.csv(file, header = TRUE)
        df_in$LY <- df_in$person_years
        df_in$bau_LY <- df_in$bau_person_years

        df_acmr_in <- df_in[df_in$year %in% rate_years
                            & df_in$age %in% rate_ages
                            & df_in$sex %in% rate_sex,
                            c(rate_cols, 'bau_acmr', 'acmr')]
        df_acmr_in$popn <- popn
        df_acmr <- rbind(df_acmr, df_acmr_in)
        df_yldr_in <- df_in[df_in$year %in% rate_years
                            & df_in$age %in% rate_ages
                            & df_in$sex %in% rate_sex,
                            c(rate_cols, 'bau_yld_rate', 'yld_rate')]
        df_yldr_in$popn <- popn
        df_yldr <- rbind(df_yldr, df_yldr_in)

        totals <- colSums(df_in[, c('bau_LY', 'LY', 'bau_HALY', 'HALY')])
        male_t <- colSums(df_in[df_in$sex == 'male',
                                c('bau_LY', 'LY', 'bau_HALY', 'HALY')])
        female_t <- colSums(df_in[df_in$sex == 'female',
                                  c('bau_LY', 'LY', 'bau_HALY', 'HALY')])
        df_ly_in <- as.data.frame(t(totals))
        df_ly_in$sex <- 'All'
        df_ly_in_m <- as.data.frame(t(male_t))
        df_ly_in_m$sex <- 'male'
        df_ly_in_f <- as.data.frame(t(female_t))
        df_ly_in_f$sex <- 'female'
        df_ly_in <- rbind(df_ly_in, df_ly_in_m, df_ly_in_f)
        df_ly_in$popn <- popn
        df_ly_in$draw_number <- draw_number
        df_ly <- rbind(df_ly, df_ly_in)
        if (verbose) {
            cat(' done\n')
        }
    }

    ## Calculate the net LYs and HALYs for the total population --- both Maori
    ## and non-Maori --- by pairing simulations based on 'draw_number'.
    df_all_popn <- df_ly[df_ly$sex == 'All', ]
    if (length(unique(df_all_popn$popn)) > 1) {
        ## NOTE: until all simulations are complete, only retain simulations
        ## where we have results for both populations.
        max_draw_nm <- max(df_all_popn$draw_number[
                                           df_all_popn$popn == 'non-maori'])
        max_draw_m <- max(df_all_popn$draw_number[
                                          df_all_popn$popn == 'maori'])
        max_draw <- min(max_draw_nm, max_draw_m)
        df_all_popn <- df_all_popn[df_all_popn$draw_number <= max_draw, ]

        df_all_popn <- aggregate(cbind(bau_LY, LY, bau_HALY, HALY)
                                 ~ draw_number + sex,
                                 data = df_all_popn, FUN = sum)
        df_all_popn$popn <- 'All'
        df_ly <- rbind(df_ly, df_all_popn)
    }

    ## Calculate relative gains.
    df_ly$LY_gain <- df_ly$LY - df_ly$bau_LY
    df_ly$LY_pcnt <- 100 * df_ly$LY_gain / df_ly$bau_LY

    df_ly$HALY_gain <- df_ly$HALY - df_ly$bau_HALY
    df_ly$HALY_pcnt <- 100 * df_ly$HALY_gain / df_ly$bau_HALY

    ## NOTE: multiply ACMR by 1e5
    df_acmr$acmr <- 1e5 * df_acmr$acmr
    df_acmr$bau_acmr <- 1e5 * df_acmr$bau_acmr

    df_acmr$acmr_gain <- df_acmr$acmr - df_acmr$bau_acmr
    df_acmr$acmr_pcnt <- 100 * df_acmr$acmr_gain / df_acmr$acmr

    df_yldr$yld_rate_gain <- df_yldr$yld_rate - df_yldr$bau_yld_rate
    df_yldr$yld_rate_pcnt <- 100 * df_yldr$yld_rate_gain / df_yldr$yld_rate

    df_ly <- calculate_ci(
        df_ly,
        cbind(bau_LY, LY_gain, LY_pcnt, bau_HALY, HALY_gain, HALY_pcnt)
        ~ popn + sex)
    df_acmr <- calculate_ci(
        df_acmr,
        cbind(bau_acmr, acmr_gain, acmr_pcnt) ~ popn + year + age + sex)
    df_yldr <- calculate_ci(
        df_yldr,
        cbind(bau_yld_rate, yld_rate_gain, yld_rate_pcnt)
        ~ popn + year + age + sex)

    ## Round numbers
    df_ly$bau_LY <- round(df_ly$bau_LY)
    df_ly$LY_gain <- round(df_ly$LY_gain)
    df_ly$LY_pcnt <- round(df_ly$LY_pcnt, digits = 2)

    df_ly$bau_HALY <- round(df_ly$bau_HALY)
    df_ly$HALY_gain <- round(df_ly$HALY_gain)
    df_ly$HALY_pcnt <- round(df_ly$HALY_pcnt, digits = 2)

    df_acmr$bau_acmr <- round(df_acmr$bau_acmr, digits = 1)
    df_acmr$acmr_gain <- round(df_acmr$acmr_gain, digits = 1)
    df_acmr$acmr_pcnt <- round(df_acmr$acmr_pcnt, digits = 2)

    df_yldr$bau_yld_rate <- round(df_yldr$bau_yld_rate, digits = 4)
    df_yldr$yld_rate_gain <- round(df_yldr$yld_rate_gain, digits = 4)
    df_yldr$yld_rate_pcnt <- round(df_yldr$yld_rate_pcnt, digits = 2)

    ## Add scenario metadata.
    df_ly$delay <- delay
    df_ly$tob_prev <- tob_prev
    df_ly$interv <- interv
    df_acmr$delay <- delay
    df_acmr$tob_prev <- tob_prev
    df_acmr$interv <- interv
    df_yldr$delay <- delay
    df_yldr$tob_prev <- tob_prev
    df_yldr$interv <- interv

    return(list(df_ly = df_ly, df_acmr = df_acmr, df_yldr = df_yldr))
}


calculate_ci <- function(df, formula, width=0.95) {
    pr_lower <- 0.5 * (1 - width)
    pr_upper <- 0.5 * (1 + width)
    pr_median <- 0.5

    df_lower <- do.call(
        data.frame,
        aggregate(formula, data = df,
                  FUN = function(x) quantile(x, probs = pr_lower)))
    df_upper <- do.call(
        data.frame,
        aggregate(formula, data = df,
                  FUN = function(x) quantile(x, probs = pr_upper)))
    df_median <- do.call(
        data.frame,
        aggregate(formula, data = df,
                  FUN = function(x) quantile(x, probs = pr_median)))
    df_mean <- do.call(
        data.frame,
        aggregate(formula, data = df,
                  FUN = mean))

    df_lower$bound <- 'lower'
    df_upper$bound <- 'upper'
    df_median$bound <- 'median'
    df_mean$bound <- 'mean'

    rbind(df_lower, df_upper, df_median, df_mean)
}


run <- function() {
    populations <- c('maori', 'non-maori')
    interventions <- c('erad', 'tfg', 'tax')
    bau_delay <- c(20, 0, 20)
    bau_tob_prev <- c('decreasing', 'decreasing', 'constant')

    df_all_ly <- NULL
    df_all_acmr <- NULL
    df_all_yldr <- NULL

    for (bau_ix in seq_along(bau_delay)) {
        delay <- bau_delay[bau_ix]
        tob_prev <- bau_tob_prev[bau_ix]
        for (interv in interventions) {
            dfs <- summarise_files(delay, tob_prev, interv)
            if (is.null(dfs)) {
                next
            }
            df_all_ly <- rbind(df_all_ly, dfs$df_ly)
            df_all_acmr <- rbind(df_all_acmr, dfs$df_acmr)
            df_all_yldr <- rbind(df_all_yldr, dfs$df_yldr)
        }
    }

    df_all_ly <- df_all_ly[order(df_all_ly$delay, df_all_ly$tob_prev,
                                 df_all_ly$interv, df_all_ly$sex,
                                 df_all_ly$popn, df_all_ly$bound), ]

    df_all_acmr <- df_all_acmr[order(df_all_acmr$delay, df_all_acmr$tob_prev,
                                     df_all_acmr$interv, df_all_acmr$sex,
                                     df_all_acmr$age, df_all_acmr$year,
                                     df_all_acmr$popn, df_all_acmr$bound), ]

    df_all_yldr <- df_all_yldr[order(df_all_yldr$delay, df_all_yldr$tob_prev,
                                     df_all_yldr$interv, df_all_yldr$sex,
                                     df_all_yldr$age, df_all_yldr$year,
                                     df_all_yldr$popn, df_all_yldr$bound), ]

    return(list(df_all_ly = df_all_ly, df_all_acmr = df_all_acmr,
                df_all_yldr = df_all_yldr))
}


to_wide <- function(df) {
    ## NOTE: list ID columns in the order in which we want to sort the rows.
    all_id_cols <- c('delay', 'tob_prev', 'interv',
                     'sex', 'age', 'year', 'popn')

    id_cols <- all_id_cols[all_id_cols %in% names(df)]
    time_col <- 'bound'
    other_cols <- names(df)[! (names(df) %in% c(time_col, id_cols))]

    df_wide <- reshape(df,
                       idvar = id_cols,
                       timevar = time_col,
                       direction = 'wide')

    ## Sort the rows by each of the ID columns in turn.
    df_wide <- df_wide[
        do.call(order, lapply(id_cols, function(col) df_wide[[col]])), ]

    ## Retain the same column order as the original data frame.
    out_cols <- c()
    for (col in names(df)) {
        if (col == 'bound') {
            next
        } else if (col %in% other_cols) {
            out_cols <- c(out_cols,
                          paste0(col, '.lower'),
                          paste0(col, '.upper'),
                          paste0(col, '.median'),
                          paste0(col, '.mean'))
        } else {
            out_cols <- c(out_cols, col)
        }
    }

    return(df_wide[, out_cols])
}


main <- function(args = NULL) {
    dfs <- run()

    df_ly <- to_wide(
        dfs$df_all_ly[, c('popn', 'sex',
                          'bau_LY', 'LY_gain', 'LY_pcnt',
                          'bound', 'delay', 'tob_prev', 'interv')])
    df_haly <- to_wide(
        dfs$df_all_ly[, c('popn', 'sex',
                          'bau_HALY', 'HALY_gain', 'HALY_pcnt',
                          'bound', 'delay', 'tob_prev', 'interv')])
    df_acmr <- to_wide(dfs$df_all_acmr)
    df_yldr <- to_wide(dfs$df_all_yldr)

    ## Save each table to disk.
    write.table(df_ly, 'uncertainty-ly.ssv',
                row.names = FALSE, quote = FALSE)
    write.table(df_haly, 'uncertainty-haly.ssv',
                row.names = FALSE, quote = FALSE)
    write.table(df_acmr, 'uncertainty-acmr.ssv',
                row.names = FALSE, quote = FALSE)
    write.table(df_yldr, 'uncertainty-yldr.ssv',
                row.names = FALSE, quote = FALSE)

    options(width = 120)
    print(df_ly)
    cat('\n')
    print(df_haly)
    cat('\n')
    print(df_acmr)
    cat('\n')
    print(df_yldr)
    cat('\n')

    return(invisible(0))
}


##
## Only call main() if this script is being run from the command-line.
##
if (! interactive()) {
    file.arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
    file.path <- substring(file.arg, 8)
    file.regex <- paste0('^(|.*/)', SCRIPT)
    if (grepl(file.regex, file.path)) {
        quit(status = main(args = commandArgs(TRUE)))
    }
}
