#!/usr/bin/Rscript --vanilla

SCRIPT <- 'tabulate.R'


##
## Return only rows pertaining to a single BAU scenario.
##
## - bau:   "const" or "decr"
## - delay: 0 or 20
##
bau_only <- function(df, bau = c('const', 'decr'), delay) {
    bau <- match.arg(bau)
    bau <- c('const' = 'constant', 'decr' = 'decreasing')[bau]
    df[df$delay == delay
      & df$tob_prev == bau,]
}


##
## Calculate the mean LYs for each population.
##
## - bau:   "const" or "decr"
## - delay: 0 or 20
##
mean_values <- function(bau, delay) {
    nm_base_name <- paste0('../exp_non-maori_bau-', bau, '_delay-',
                           delay, 'yrs_')
    m_base_name <- paste0('../exp_maori_bau-', bau, '_delay-',
                          delay, 'yrs_')

    nm_erad <- read.csv(paste0(nm_base_name, 'erad_mm.csv'))
    nm_tax <- read.csv(paste0(nm_base_name, 'tax_mm.csv'))
    nm_tfg <- read.csv(paste0(nm_base_name, 'tfg_mm.csv'))
    m_erad <- read.csv(paste0(m_base_name, 'erad_mm.csv'))
    m_tax <- read.csv(paste0(m_base_name, 'tax_mm.csv'))
    m_tfg <- read.csv(paste0(m_base_name, 'tfg_mm.csv'))

    ## Extract the ACMR and YLDR for cohorts of interest.
    age <- 63
    years <- c(2041, 2061)
    nm_erad_rates <- nm_erad[nm_erad$age == age & nm_erad$year %in% years, ]
    nm_tax_rates <- nm_tax[nm_tax$age == age & nm_tax$year %in% years, ]
    nm_tfg_rates <- nm_tfg[nm_tfg$age == age & nm_tfg$year %in% years, ]
    m_erad_rates <- m_erad[m_erad$age == age & m_erad$year %in% years, ]
    m_tax_rates <- m_tax[m_tax$age == age & m_tax$year %in% years, ]
    m_tfg_rates <- m_tfg[m_tfg$age == age & m_tfg$year %in% years, ]

    nm_erad_rates$popn <- 'non-maori'
    nm_erad_rates$interv <- 'erad'
    nm_tax_rates$popn <- 'non-maori'
    nm_tax_rates$interv <- 'tax'
    nm_tfg_rates$popn <- 'non-maori'
    nm_tfg_rates$interv <- 'tfg'
    m_erad_rates$popn <- 'maori'
    m_erad_rates$interv <- 'erad'
    m_tax_rates$popn <- 'maori'
    m_tax_rates$interv <- 'tax'
    m_tfg_rates$popn <- 'maori'
    m_tfg_rates$interv <- 'tfg'
    df_rates <- rbind(nm_erad_rates, m_erad_rates,
                      nm_tax_rates, m_tax_rates,
                      nm_tfg_rates, m_tfg_rates)
    df_rates <- df_rates[, c('popn', 'year', 'age', 'sex',
                             'bau_acmr', 'acmr',
                             'bau_yld_rate', 'yld_rate',
                             'interv')]
    df_rates$delay <- delay
    df_rates$bau_acmr <- 1e5 * df_rates$bau_acmr
    df_rates$acmr <- 1e5 * df_rates$acmr

    ## Extract the LYs and HALYs.
    nm_erad <- aggregate(
        cbind(person_years, bau_person_years, HALY, bau_HALY) ~ sex,
        data = nm_erad, FUN = sum)
    nm_erad$popn <- 'non-maori'
    nm_erad$interv <- 'erad'
    nm_tax <- aggregate(
        cbind(person_years, bau_person_years, HALY, bau_HALY) ~ sex,
        data = nm_tax, FUN = sum)
    nm_tax$popn <- 'non-maori'
    nm_tax$interv <- 'tax'
    nm_tfg <- aggregate(
        cbind(person_years, bau_person_years, HALY, bau_HALY) ~ sex,
        data = nm_tfg, FUN = sum)
    nm_tfg$popn <- 'non-maori'
    nm_tfg$interv <- 'tfg'

    m_erad <- aggregate(
        cbind(person_years, bau_person_years, HALY, bau_HALY) ~ sex,
        data = m_erad, FUN = sum)
    m_erad$popn <- 'maori'
    m_erad$interv <- 'erad'
    m_tax <- aggregate(
        cbind(person_years, bau_person_years, HALY, bau_HALY) ~ sex,
        data = m_tax, FUN = sum)
    m_tax$popn <- 'maori'
    m_tax$interv <- 'tax'
    m_tfg <- aggregate(
        cbind(person_years, bau_person_years, HALY, bau_HALY) ~ sex,
        data = m_tfg, FUN = sum)
    m_tfg$popn <- 'maori'
    m_tfg$interv <- 'tfg'

    df <- rbind(nm_erad, m_erad, nm_tax, m_tax, nm_tfg, m_tfg)
    df$delay <- delay
    df$LY <- df$person_years
    df$bau_LY <- df$bau_person_years
    df$person_years <- NULL
    df$bau_person_years <- NULL

    df_net <- aggregate(
        cbind(bau_LY, LY, bau_HALY, HALY) ~ interv + delay,
        data = df, FUN = sum)
    df_net$sex <- 'All'
    df_net$popn <- 'All'

    df <- rbind(df, df_net)

    list(df_mean = df, df_rates = df_rates)
}


ly_gain <- function(df, interv) {
    df <- df[df$interv == interv, ]
    paste(format(round(df$LY_gain), big.mark = ","),
          " (", format(df$LY_pcnt, nsmall = 2), "%)", sep = "")
}


haly_gain <- function(df, interv) {
    df <- df[df$interv == interv, ]
    paste(format(round(df$HALY_gain), big.mark = ","),
          " (", format(df$HALY_pcnt, nsmall = 2), "%)", sep = "")
}

acmr_gain <- function(df, interv) {
    df <- df[df$interv == interv, ]
    paste(format(round(df$acmr_gain), big.mark = ","),
          " (", format(df$acmr_pcnt, nsmall = 2), "%)", sep = "")
}


yldr_gain <- function(df, interv) {
    df <- df[df$interv == interv, ]
    paste(format(round(df$yldr_gain, 6), nsmall = 6),
          " (", format(df$yldr_pcnt, nsmall = 2), "%)", sep = "")
}

ly_gain_ci <- function(df, interv, sep='-') {
    df <- df[df$interv == interv, ]
    paste(format(round(df$LY_gain.lower), big.mark = ","),
          sep,
          format(round(df$LY_gain.upper), big.mark = ","),
          " (",
          format(df$LY_pcnt.lower, nsmall = 2),
          sep,
          format(df$LY_pcnt.upper, nsmall = 2),
          "%)", sep = "")
}


haly_gain_ci <- function(df, interv, sep='-') {
    df <- df[df$interv == interv, ]
    paste(format(round(df$HALY_gain.lower), big.mark = ","),
          sep,
          format(round(df$HALY_gain.upper), big.mark = ","),
          " (",
          format(df$HALY_pcnt.lower, nsmall = 2),
          sep,
          format(df$HALY_pcnt.upper, nsmall = 2),
          "%)", sep = "")
}


acmr_gain_ci <- function(df, interv, sep=' - ') {
    df <- df[df$interv == interv, ]
    paste(format(round(df$acmr_gain.lower), big.mark = ","),
          sep,
          format(round(df$acmr_gain.upper), big.mark = ","),
          " (",
          format(df$acmr_pcnt.lower, nsmall = 2),
          sep,
          format(df$acmr_pcnt.upper, nsmall = 2),
          "%)", sep = "")
}


yldr_gain_ci <- function(df, interv, sep=' - ') {
    df <- df[df$interv == interv, ]
    paste(format(round(df$yld_rate_gain.lower, 6), nsmall = 6),
          sep,
          format(round(df$yld_rate_gain.upper, 6), nsmall = 6),
          " (",
          format(df$yld_rate_pcnt.lower, nsmall = 2),
          sep,
          format(df$yld_rate_pcnt.upper, nsmall = 2),
          "%)", sep = "")
}


##
## Compare each tobacco intervention to the BAU, showing the difference in:
## LYs, HALYs, ACMR, and YLDR.
##
build_table <- function(bau = 'decr', delay = 20) {
    df_ly <- read.table('uncertainty-ly.ssv', header = TRUE, as.is = TRUE)
    df_haly <- read.table('uncertainty-haly.ssv', header = TRUE, as.is = TRUE)
    df_acmr <- read.table('uncertainty-acmr.ssv', header = TRUE, as.is = TRUE)
    df_yldr <- read.table('uncertainty-yldr.ssv', header = TRUE, as.is = TRUE)

    ## Only retain results for the standard BAU scenario.
    df_ly <- bau_only(df_ly, bau, delay)
    df_haly <- bau_only(df_haly, bau, delay)
    df_acmr <- bau_only(df_acmr, bau, delay)
    df_yldr <- bau_only(df_yldr, bau, delay)

    ## Only retain results for the cohorts of interest.
    ly_mask <- ((df_ly$popn == 'All' & df_ly$sex == 'All')
        | (df_ly$popn == 'maori' & df_ly$sex == 'female')
        | (df_ly$popn == 'non-maori' & df_ly$sex == 'male'))
    df_ly <- df_ly[ly_mask, ]
    haly_mask <- ((df_haly$popn == 'All' & df_haly$sex == 'All')
        | (df_haly$popn == 'maori' & df_haly$sex == 'female')
        | (df_haly$popn == 'non-maori' & df_haly$sex == 'male'))
    df_haly <- df_haly[haly_mask, ]

    ## Calculate the mean LY / HALY values and gains.
    dfs <- mean_values(bau, delay)
    df_mean <- dfs$df_mean
    df_ly <- merge(df_ly, df_mean, sort = FALSE)
    df_haly <- merge(df_haly, df_mean, sort = FALSE)
    df_ly$LY_gain <- df_ly$LY - df_ly$bau_LY
    df_ly$LY_pcnt <- round(100 * df_ly$LY_gain / df_ly$bau_LY, 2)
    df_haly$HALY_gain <- df_haly$HALY - df_haly$bau_HALY
    df_haly$HALY_pcnt <- round(100 * df_haly$HALY_gain / df_haly$bau_HALY, 2)

    ## Calculate the ACMR / YLDR values and gains.
    df_rates <- dfs$df_rates
    df_acmr <- merge(df_acmr, df_rates, sort = FALSE)
    df_acmr$acmr_gain <- df_acmr$acmr - df_acmr$bau_acmr
    df_acmr$acmr_pcnt <- round(100 * df_acmr$acmr_gain / df_acmr$bau_acmr, 2)
    df_yldr <- merge(df_yldr, df_rates, sort = FALSE)
    df_yldr$yldr_gain <- df_yldr$yld_rate - df_yldr$bau_yld_rate
    df_yldr$yldr_pcnt <- round(100 * df_yldr$yldr_gain / df_yldr$bau_yld_rate, 2)

    ## Retain cohorts of interest and sort rows.
    df_acmr <- df_acmr[(df_acmr$popn == 'maori' & df_acmr$sex == 'female')
                       | (df_acmr$popn == 'non-maori' & df_acmr$sex == 'male'), ]
    df_acmr <- df_acmr[order(df_acmr$interv, df_acmr$popn, df_acmr$year), ]
    df_yldr <- df_yldr[(df_yldr$popn == 'maori' & df_yldr$sex == 'female')
                       | (df_yldr$popn == 'non-maori' & df_yldr$sex == 'male'), ]
    df_yldr <- df_yldr[order(df_yldr$interv, df_yldr$popn, df_yldr$year), ]

    ## Assemble the table content, one column at a time.
    output <- c('LYs', '', '',
                'HALYs', '', '',
                'ACMR', '', '', '',
                'YLDR', '', '', '')
    demographic <- c('Total', 'Maori female', 'Non-Maori male',
                     'Total', 'Maori female', 'Non-Maori male',
                     'Maori female', 'Maori female',
                     'Non-Maori male', 'Non-Maori male',
                     'Maori female', 'Maori female',
                     'Non-Maori male', 'Non-Maori male')
    age <- c(rep('', 6),
             df_acmr$age[1:4],
             df_yldr$age[1:4])
    calendar_year <- c(rep('', 6),
                       df_acmr$year[1:4],
                       df_yldr$year[1:4])
    bau <- c(format(df_ly$bau_LY[1:3], nsmall = 0, big.mark = ","),
             format(df_haly$bau_HALY[1:3], nsmall = 0, big.mark = ","),
             format(round(df_acmr$bau_acmr[1:4])),
             format(round(df_yldr$bau_yld_rate[1:4], 6), nsmall = 6))
    tobacco_eradication <- c(
        ly_gain(df_ly, 'erad'),
        haly_gain(df_haly, 'erad'),
        acmr_gain(df_acmr, 'erad'),
        yldr_gain(df_yldr, 'erad'))
    tobacco_eradication_ci <- c(
        ly_gain_ci(df_ly, 'erad'),
        haly_gain_ci(df_haly, 'erad'),
        acmr_gain_ci(df_acmr, 'erad'),
        yldr_gain_ci(df_yldr, 'erad'))
    tobacco_tax <- c(
        ly_gain(df_ly, 'tax'),
        haly_gain(df_haly, 'tax'),
        acmr_gain(df_acmr, 'tax'),
        yldr_gain(df_yldr, 'tax'))
    tobacco_tax_ci <- c(
        ly_gain_ci(df_ly, 'tax'),
        haly_gain_ci(df_haly, 'tax'),
        acmr_gain_ci(df_acmr, 'tax'),
        yldr_gain_ci(df_yldr, 'tax'))
    tobacco_free_generation <- c(
        ly_gain(df_ly, 'tfg'),
        haly_gain(df_haly, 'tfg'),
        acmr_gain(df_acmr, 'tfg'),
        yldr_gain(df_yldr, 'tfg'))
    tobacco_free_generation_ci <- c(
        ly_gain_ci(df_ly, 'tfg'),
        haly_gain_ci(df_haly, 'tfg'),
        acmr_gain_ci(df_acmr, 'tfg'),
        yldr_gain_ci(df_yldr, 'tfg'))

    df_out <- data.frame(
        output = output,
        demographic = demographic,
        age = age,
        calendar_year = calendar_year,
        bau = bau,
        tobacco_eradication = tobacco_eradication,
        tobacco_eradication_ci = tobacco_eradication_ci,
        tobacco_tax = tobacco_tax,
        tobacco_tax_ci = tobacco_tax_ci,
        tobacco_free_generation = tobacco_free_generation,
        tobacco_free_generation_ci = tobacco_free_generation_ci
    )

    names(df_out) <- c(
        'Output',
        'Demographic',
        'Age',
        'Calendar Year',
        'BAU',
        'Tobacco Eradication',
        'Tobacco Eradication CIs',
        'Tobacco Tax',
        'Tobacco Tax CIs',
        'Tobacco-Free Generation',
        'Tobacco-Free Generation CIs')

    df_out
}


##
## The script entry-point.
##
main <- function(args = NULL) {
    table_3 <- build_table(bau = 'decr', delay = 20)
    table_4a <- build_table(bau = 'const', delay = 20)
    table_4b <- build_table(bau = 'decr', delay = 0)
    write.csv(table_3, 'table_3.csv', row.names = FALSE)
    write.csv(table_4a, 'table_4a_constant_prevalence.csv', row.names = FALSE)
    write.csv(table_4b, 'table_4b_immediate_recovery.csv', row.names = FALSE)
    invisible(0)
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
