#!/usr/bin/Rscript --vanilla

SCRIPT <- 'print-tables.R'

exp_tags_from_filename <- function(filename) {
    popn <- c('maori' = 'Maori', 'non-maori' = 'non-Maori')
    bau <- c('decreasing' = 'Decreasing', 'constant' = 'Constant')
    delay <- c('0-years' = 'Immediate', '20-years' = 'Delayed')
    interv <- c('erad' = 'Tobacco Eradication',
                'tax' = 'Tobacco Tax',
                'tfg' = 'Tobacco-Free Generation')
    reduce <- c('acmr' = 'Reduce ACMR by 5%',
                'chd' = 'Reduce CHD by 5%')

    filename <- basename(filename)
    name_parts <- strsplit(filename, '_')[[1]]
    num_parts <- length(name_parts)

    result <- NULL
    valid <- TRUE
    if (num_parts < 1) {
        valid <- FALSE
    }
    valid <- valid && endsWith(name_parts[num_parts], '.csv')
    if (name_parts[1] == 'mslt' && name_parts[2] == 'tobacco') {
        valid <- valid && num_parts == 7
        valid <- valid && name_parts[3] %in% names(popn)
        valid <- valid && name_parts[4] %in% names(delay)
        valid <- valid && name_parts[5] %in% names(bau)
        valid <- valid && name_parts[6] %in% names(interv)
        if (valid) {
            result <- c('popn' = popn[[name_parts[3]]],
                        'delay' = delay[[name_parts[4]]],
                        'bau' = bau[[name_parts[5]]],
                        'interv' = interv[[name_parts[6]]])
        }
    } else if (name_parts[1] == 'mslt' && name_parts[2] == 'reduce') {
        valid <- valid && num_parts == 4
        valid <- valid && name_parts[3] %in% names(reduce)
        if (valid) {
            result <- c('popn' = 'non-Maori',
                        'reduce' = reduce[[name_parts[3]]])
        }
    }

    result
}

load_csv <- function(filename) {
    df <- read.csv(filename, header = TRUE, as.is = TRUE)
    tags <- exp_tags_from_filename(filename)
    tag_names <- names(tags)
    for (i in 1:length(tags)) {
        df[[tag_names[i]]] <- tags[[i]]
    }
    df
}

load_csv_files <- function(filenames) {
    df <- NULL
    for (filename in filenames) {
        df_in <- load_csv(filename)
        df <- rbind(df, df_in)
    }
    df
}

load_data <- function(data_dir = 'results') {
    pattern <- 'mslt_tobacco_.*_(erad|tax|tfg)_mm\\.csv'
    filenames <- list.files(path = data_dir, pattern = pattern,
                            full.names = TRUE)
    df <- load_csv_files(filenames)
}


comparison_table <- function(df, column, ix_columns, prec = 6) {
    df_bau <- df[df$interv == 'BAU', ]
    df_int <- df[df$interv != 'BAU', ]
    df_out <- df_bau
    for (interv in unique(df_int$interv)) {
        mask <- df_int$interv == interv
        df_tmp <- df_int[mask, ]
        df_tmp[[column]] <- df_tmp[[column]] - df_bau[[column]]
        pcnt <- 100 * df_tmp[[column]] / df_bau[[column]]
        df_tmp[[column]] <- sprintf('%.*f (%0.2f%%)', prec, df_tmp[[column]],
                                    pcnt)
        names(df_tmp)[names(df_tmp) == column] <- interv
        df_out <- merge(df_out,
                        df_tmp[, c(ix_columns, interv)])
    }
    df_out
}


print_LY_HALY_ACMR_YLD_table <- function(df) {
    prev_width <- getOption('width')
    options(width = 1e3)

    common_cols <- c('age', 'sex', 'year_of_birth',
                     'popn', 'bau', 'delay', 'interv')

    bau_cols <- c(common_cols,
                  'bau_acmr', 'bau_population',
                  'bau_person_years',
                  'bau_yld_rate', 'bau_HALY',
                  'bau_pr_death', 'bau_deaths')

    int_cols <- c(common_cols,
                  'acmr', 'population',
                  'person_years',
                  'yld_rate', 'HALY',
                  'pr_death', 'deaths')

    df_bau <- df[, bau_cols]
    names(df_bau) <- gsub('^bau_', '', names(df_bau))
    df_bau$interv <- 'BAU'
    df_bau <- df_bau[! duplicated(df_bau), ]

    df_int <- df[, int_cols]
    names(df_int) <- gsub('^int_', '', names(df_int))

    df0 <- rbind(df_bau, df_int)
    df <- df0

    ## Print LYs for the total population, Maori females, non-Maori males.
    df_ly_net <- aggregate(person_years ~ interv, data = df, FUN = sum)
    df_ly_net$popn <- '    Total'
    df_ly_net$sex <- '      '
    df_out <- comparison_table(df_ly_net, 'person_years',
                               c('popn', 'sex'),
                               prec = 0)
    df_out$interv <- NULL
    print(df_out, row.names = FALSE)
    df_ly <- aggregate(person_years ~ popn + sex + interv, data = df,
                       FUN = sum)
    df_ly <- df_ly[(df_ly$popn == 'non-Maori' & df_ly$sex == 'male')
                   | (df_ly$popn == 'Maori' & df_ly$sex == 'female'), ]
    df_out <- comparison_table(df_ly, 'person_years',
                               c('popn', 'sex'),
                               prec = 0)
    df_out$interv <- NULL
    print(df_out, row.names = FALSE)

    ## Print HALYs for the total population, Maori females, non-Maori males.
    df_haly_net <- aggregate(HALY ~ interv, data = df, FUN = sum)
    df_haly_net$popn <- '    Total'
    df_haly_net$sex <- '      '
    df_out <- comparison_table(df_haly_net, 'HALY',
                               c('popn', 'sex'),
                               prec = 0)
    df_out$interv <- NULL
    print(df_out, row.names = FALSE)
    df_haly <- aggregate(HALY ~ popn + sex + interv, data = df, FUN = sum)
    df_haly <- df_haly[(df_haly$popn == 'non-Maori' & df_haly$sex == 'male')
                   | (df_haly$popn == 'Maori' & df_haly$sex == 'female'), ]
    df_out <- comparison_table(df_haly, 'HALY',
                               c('popn', 'sex'),
                               prec = 0)
    df_out$interv <- NULL
    print(df_out, row.names = FALSE)

    ## Print ACMRs for:
    ## 62-yo non-Maori males and Maori females in 2041 and 2061
    df_acmr <- df0[df0$age == 62
                   & df0$year_of_birth %in% c(1979, 1999),
                   c('age', 'sex', 'popn', 'year_of_birth',
                     'interv', 'acmr')]
    df_out <- comparison_table(df_acmr, 'acmr',
                               c('popn', 'sex', 'year_of_birth'))
    df_out$acmr <- sprintf('%.6f', df_out$acmr)
    df_out <- df_out[(df_out$sex == 'male' & df_out$popn == 'non-Maori')
                     | (df_out$sex == 'female' & df_out$popn == 'Maori'), ]
    df_out$year_of_birth <- df_out$year_of_birth + df_out$age
    names(df_out)[names(df_out) == 'year_of_birth'] <- 'year'
    df_out$interv <- NULL
    print(df_out, row.names = FALSE)

    ## Print YLDRs for:
    ## 62-yo non-Maori males and Maori females in 2041 and 2061
    df_yldr <- df0[df0$age == 62
                   & df0$year_of_birth %in% c(1979, 1999),
                   c('age', 'sex', 'popn', 'year_of_birth',
                     'interv', 'yld_rate')]
    df_out <- comparison_table(df_yldr, 'yld_rate',
                               c('popn', 'sex', 'year_of_birth'))
    df_out$yld_rate <- sprintf('%.6f', df_out$yld_rate)
    df_out <- df_out[(df_out$sex == 'male' & df_out$popn == 'non-Maori')
                     | (df_out$sex == 'female' & df_out$popn == 'Maori'), ]
    df_out$year_of_birth <- df_out$year_of_birth + df_out$age
    names(df_out)[names(df_out) == 'year_of_birth'] <- 'year'
    df_out$interv <- NULL
    print(df_out, row.names = FALSE)

    options(width = prev_width)
    invisible(NULL)
}


load_table2_files <- function(data_dir = 'results') {
    re_mm <- 'mslt_reduce_(acmr|chd)_mm\\.csv'

    mm_files <- list.files(path = data_dir, pattern = re_mm,
                           full.names = TRUE)
    df <- load_csv_files(mm_files)

    year_of_birth <- 1959
    popn <- 'non-Maori'
    sex <- 'male'
    ages <- c(52, 53, 109, 110)

    df <- df[df$year_of_birth == year_of_birth
             & df$popn == popn
             & df$sex == sex, ]


    bau_cols <- c('age', 'bau_acmr', 'bau_population',
                  'bau_person_years',
                  'bau_yld_rate', 'bau_HALY',
                  'bau_pr_death', 'bau_deaths',
                  'reduce')

    int_cols <- c('age', 'acmr', 'population',
                  'person_years',
                  'yld_rate', 'HALY',
                  'pr_death', 'deaths',
                  'reduce')

    df_bau <- df[, bau_cols]
    names(df_bau) <- gsub('^bau_', '', names(df_bau))
    df_bau$reduce <- 'BAU'

    df_int <- df[, int_cols]
    names(df_int) <- gsub('^int_', '', names(df_int))

    df <- rbind(df_bau, df_int)
    df <- df[order(df$reduce, df$age), ]
    df <- df[! duplicated(df), ]

    df_net <- aggregate(list(LY = df$person_years, HALY = df$HALY),
                        by = list(reduce = df$reduce),
                        FUN = sum)
    names(df_net) <- c('reduce', 'NetLY', 'NetHALY')

    ## Calculate life expectancy.
    df$LE <- ave(df$person_years, df$reduce, FUN = function(x) rev(cumsum(rev(x)))) / df$population
    df$HALE <- ave(df$HALY, df$reduce, FUN = function(x) rev(cumsum(rev(x)))) / df$population

    df <- df[df$age %in% ages, ]

    df$survive <- df$population + df$deaths

    df$acmr <- round(df$acmr, 4)
    df$pr_death <- round(df$pr_death, 4)
    df$population <- round(df$population)
    df$deaths <- round(df$deaths)
    df$survive <- round(df$survive)
    df$person_years <- round(df$person_years)
    df$LE <- round(df$LE, 2)
    df$yld_rate <- round(df$yld_rate, 4)
    df$HALY <- round(df$HALY)
    df$HALE <- round(df$HALE, 2)

    df <- df[, c('age', 'acmr', 'pr_death', 'survive', 'deaths', 'population',
                 'person_years', 'LE',
                 'yld_rate',
                 'HALY', 'HALE',
                 'reduce')]

    list(df = df, df_net = df_net)
}

main <- function(args) {
    prev_width <- getOption('width')
    options(width = 1e3)
    cat('\n')
    cat('=========\n')
    cat(' Table 2 \n')
    cat('=========\n')
    cat('\n')
    dfs <- load_table2_files()
    df_tbl2 <- dfs$df
    df_tbl2_net <- dfs$df_net
    for (reduce in unique(df_tbl2$reduce)) {
        print(df_tbl2[df_tbl2$reduce == reduce, ], row.names = FALSE)
        df_tmp <- df_tbl2_net[df_tbl2_net$reduce == reduce, ]
        cat(sprintf('%s%7.0f%s%7.0f\n',
                    paste0(rep(' ', 52), collapse = ''),
                    df_tmp$NetLY[1],
                    paste0(rep(' ', 15), collapse = ''),
                    df_tmp$NetHALY[1]))
    }
    options(width = prev_width)

    df <- load_data()
    df_tbl3 <- df[df$delay == 'Delayed'
                  & df$bau == 'Decreasing', ]
    df_tbl4a <- df[df$delay == 'Delayed'
                  & df$bau == 'Constant', ]
    df_tbl4b <- df[df$delay == 'Immediate'
                   & df$bau == 'Decreasing', ]

    cat('\n')
    cat('=========\n')
    cat(' Table 3 \n')
    cat('=========\n')
    cat('\n')
    print_LY_HALY_ACMR_YLD_table(df_tbl3)

    cat('\n')
    cat('=========\n')
    cat(' Table 4 \n')
    cat('=========\n')
    cat('\n')
    print_LY_HALY_ACMR_YLD_table(df_tbl4a)
    cat('\n')
    print_LY_HALY_ACMR_YLD_table(df_tbl4b)
    cat('\n')

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
