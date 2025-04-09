
quick_sig <- function(stats_df){
  stats_df %>%
    mutate(sig_text = ifelse(p_value < 0.05, 'sig', ''),
           sig_symbol = ifelse(p_value < 0.001, '***',
                               ifelse(p_value < 0.01, '**',
                                      ifelse(p_value < 0.05, '*',
                                             'n.s.'))))
}

save_aov <- function(dir_output, prefix, anova_summary, pairs_hsd, df){
  if(!is.na(dir_output)){
    fn <- str_c(dir_output, '/', prefix, '_aov.csv')
    anova_summary %>% write_csv(fn)
    print(str_c('saved file: ', fn))
    
    fn <- str_c(dir_output, '/', prefix, '_aov_hsd.csv')
    pairs_hsd %>% write_csv(fn)
    print(str_c('saved file: ', fn))
    
    fn <- str_c(dir_output, '/', prefix, '_data.csv')
    df %>% write_csv(fn)
    print(str_c('saved file: ', fn))
  }
}


save_ttest <- function(dir_output, prefix, stats_result){
  if(!is.na(dir_output)){
    fn <- str_c(dir_output, '/', prefix, '_ttest.csv')
    stats_result %>% write_csv(fn)
    print(str_c('saved file: ', fn))
  }
}


tidy_cor_output_reshape <- function(df, n_table, var_name){
  df_tidy <- df[[n_table]] %>%
    as_tibble()
  
  df_tidy <- df_tidy %>%
    mutate(var_y = colnames(df_tidy)) %>%
    gather('var_x', !!as.name(var_name), -var_y) 
  
  
  
  return(df_tidy)
}

tidy_cor_output <- function(df){
  tidy_cor_output_reshape(df, 1, 'r') %>%
    left_join(tidy_cor_output_reshape(df, 2, 'n'), by = c("var_y", "var_x")) %>%
    left_join(tidy_cor_output_reshape(df, 3, 'p'), by = c("var_y", "var_x"))
}


get_cormatrix <- function(df){
  df_corr = rcorr(as.matrix(df))
  tidy_cor_output(df_corr)
}

tidy_cormatrix <- function(df, var_cors, var_groups = NA){
  require(Hmisc)
  
  df <- df %>%
    ungroup()
  
  
  if(sum(is.na(var_groups)) > 0){
    df <- df %>% select(var_cors)
    
    df %>%
      group_modify(~ get_cormatrix(.))
    
  } else {
    df <- df %>% select(var_cors, var_groups)
    
    df %>%
      group_by_at(var_groups) %>%
      group_modify(~ get_cormatrix(.))
    
  }
}


# ttest paired
do_ttest_paired <- function(df, var1, var2){
  t.test(
    df %>% pull(var1),
    df %>% pull(var2),
    paired = TRUE,
    alternative = 'two.sided'
  ) %>%
    broom::tidy() %>%
    mutate(ag_stat = 'do_ttest_paired')
}

tidy_ttest_paired <- function(df, dir_output, prefix, var1, var2, ...){
  stats_result <- df %>%
    group_by_(...) %>%
    do(do_ttest_paired(., var1, var2)) %>%
    mutate(var1 = var1,
           var2 = var2) %>%
    select(var1, var2, everything()) %>%
    clean_names()
  
  save_ttest(dir_output, prefix, stats_result)
  
  return(stats_result)
}




# wilcox test (independent samples, non-parameteric substitute for t test)
do_wilcox_unpaired <- function(df, var1, var2){
  wilcox.test(
    df %>% pull(var1),
    df %>% pull(var2),
    paired = FALSE,
    alternative = 'two.sided'
  ) %>%
    broom::tidy() %>%
    mutate(ag_stat = 'do_wilcox_unpaired')
}

tidy_wilcox_unpaired <- function(df, var1, var2, ...){
  df %>%
    group_by_(...) %>%
    do(do_wilcox_unpaired(., var1, var2))
}



## anova one between --------------------------------------------------------------------------------------------------
aov_one_between <-  function(df, dir_output, prefix, var_id, var_dv, var_between){
  
  anova_model <- aov_ez(
    id = var_id,
    dv = var_dv,
    data = df,
    between = var_between,
    anova_table=list(correction="none", es = "none"))
  
  
  anova_summary <-  tibble(
    var_dv = var_dv,
    var_between = var_between,
    df_num = summary(anova_model)[[1]],
    df_den = summary(anova_model)[[2]],
    mse = summary(anova_model)[[3]],
    f = summary(anova_model)[[4]],
    p_value = summary(anova_model)[[5]]
  )
  
  interaction <- paste("~", var_between, sep = "")
  
  emm <- emmeans(anova_model, formula(interaction))
  
  
  pairs_hsd <- pairs(emm) %>%
    as_tibble() %>%
    separate(contrast, into = c("cont1", "cont2"), "-") %>%
    mutate(cont1 = cont1%>%str_trim(),
           cont2 = cont2%>%str_trim()
    ) %>%
    clean_names()
  
  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)
  
  return(list(anova_model, anova_summary, pairs_hsd))
  
}


## rm anova one within ------------------------------------------------------------------------------------------------
aov_rm_one_within <- function(df, dir_output, prefix, var_id, var_dependent, var_within){
  require(afex)
  
  anova_model <- aov_ez(
    id = var_id, 
    dv = var_dependent, 
    data = df,
    within = var_within,
    anova_table=list(correction="none", es = "none"))
  
  
  label_anova <- c("intercept", var_within) %>% as_tibble()
  
  anova_summary <-  summary(anova_model)[[4]][] %>% 
    as_tibble() %>%
    bind_cols(label_anova,.) %>%
    rename("effect" = "value",
           "ss" = "Sum Sq",
           "df_num" = "num Df",
           "ss_error" = "Error SS",
           "df_den" = "den Df",
           "f" = "F value",
           "p_value" = "Pr(>F)") 
  
  anova_summary <- anova_summary %>%
    mutate(stat = 'aov_rm_one_within',
           var_dv      = var_dependent,
           var_id      = var_id,
           var_within  = var_within) %>%
    select(stat, var_dv, var_id, var_within, everything()) %>%
    quick_sig()
  
  
  interaction <- paste("~", var_within, sep = "")
  
  emm <- emmeans(anova_model, formula(interaction))
  
  pairs_hsd <- pairs(emm) %>%
    as_tibble() %>%
    separate(contrast, into = c('contrast1_within', 'contrast2_within'), sep = '-') %>%
    mutate_if(is.character, str_trim) %>%
    rename('p_value' = 'p.value',
           't_ratio' = 't.ratio')  %>%
    quick_sig()
  
  # calculate Cohen's d effect size
  pairs_effect_size <- eff_size(emm, sigma = sqrt(mean(sigma(anova_model$lm)^2)), edf = df.residual(anova_model$lm))%>%
    as_tibble() %>%
    separate(contrast, into = c('contrast1_within', 'contrast2_within'), sep = '-') %>%
    mutate_if(is.character, str_trim) %>%
    clean_names() %>%
    rename(es_se = se,
           es_lower_cl = lower_cl,
           es_upper_cl = upper_cl) %>%
    select(-df)
  
  
  pairs_hsd <- pairs_hsd %>%
    left_join(pairs_effect_size, by = join_by(contrast1_within, contrast2_within))
  
  
  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)
  
  return(list(anova_model, anova_summary, pairs_hsd))
}

## rm anova two between ------------------------------------------------------------------------------------------------
aov_rm_two_between <- function(df, dir_output, prefix, var_id, var_dependent, var_between1, var_between2){
  require(afex)
  
  anova_model <- aov_ez(
    id = var_id, 
    dv = var_dependent, 
    data = df,
    between = c(var_between1, var_between2),
    anova_table=list(correction="none", es = "none"))
  
  label_anova <- c(var_between1,
                   var_between2,
                   str_c(var_between1, '*', var_between2)
  ) %>%
    as_tibble()
  
  
  
  anova_summary <-  summary(anova_model)[] %>%
    as_tibble() %>%
    bind_cols(label_anova,.) %>%
    rename("effect" = "value",
           "df_num" = "num Df",
           "df_den" = "den Df",
           'mse' = 'MSE',
           "f" = "F",
           "p_value" = "Pr(>F)")
  
  
  anova_summary <- anova_summary %>%
    mutate(stat = 'aov_rm_two_between',
           var_dv   = var_dependent,
           var_id   = var_id,
           var_between1 = var_between1,
           var_between2 = var_between2
    ) %>%
    select(stat, var_dv, everything()) %>%
    quick_sig()
  
  interaction <- paste("~", var_between1, '*', var_between2, sep = "")
  
  emm <- emmeans(anova_model, formula(interaction))
  
  pairs_hsd <- pairs(emm) %>%
    as_tibble() %>%
    separate(contrast, into = c('contrast1', 'contrast2'), sep = '-') %>%
    mutate_if(is.character, str_trim) %>%
    separate(contrast1, into = c('contrast1_between1', 'contrast1_between2'), sep = ' ') %>%
    separate(contrast2, into = c('contrast2_between1', 'contrast2_between2'), sep = ' ') %>%
    mutate_if(is.character, str_trim) %>%
    rename('p_value' = 'p.value',
           't_ratio' = 't.ratio')  %>%
    quick_sig()
  
  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)
  
  return(list(anova_model, anova_summary, pairs_hsd))
  
}

## rm anova two within ------------------------------------------------------------------------------------------------
aov_rm_three_within <- function(df, dir_output, prefix, var_id, var_dependent, var_within1, var_within2, var_within3){
  require(afex)
  
  anova_model <- aov_ez(
    id = var_id,
    dv = var_dependent,
    data = df,
    within = c(var_within1, var_within2, var_within3),
    anova_table=list(correction="none", es = "none"))
  
  label_anova <- c("intercept",
                   var_within1,
                   var_within2,
                   var_within3,
                   str_c(var_within1, '*', var_within2),
                   str_c(var_within1, '*', var_within3),
                   str_c(var_within2, '*', var_within3),
                   str_c(var_within1, '*', var_within2, '*', var_within3)
  ) %>%
    as_tibble()
  
  anova_summary <-  summary(anova_model)[[4]][] %>%
    as_tibble() %>%
    bind_cols(label_anova,.) %>%
    rename("effect" = "value",
           "ss" = "Sum Sq",
           "df_num" = "num Df",
           "ss_error" = "Error SS",
           "df_den" = "den Df",
           "f" = "F value",
           "p_value" = "Pr(>F)")
  
  
  anova_summary <- anova_summary %>%
    mutate(stat = 'aov_rm_three_within',
           var_dv   = var_dependent,
           var_id   = var_id,
           var_within1 = var_within1,
           var_within2 = var_within2
    ) %>%
    select(stat, var_dv, everything()) %>%
    quick_sig()
  
  interaction <- paste("~", var_within1, '*', var_within2, '*', var_within3, sep = "")
  
  emm <- emmeans(anova_model, formula(interaction))
  
  pairs_hsd <- pairs(emm) %>%
    as_tibble() %>%
    separate(contrast, into = c('contrast1', 'contrast2'), sep = '-') %>%
    mutate_if(is.character, str_trim) %>%
    separate(contrast1, into = c('contrast1_within1', 'contrast1_within2', 'contrast1_within3'), sep = ' ') %>%
    separate(contrast2, into = c('contrast2_within1', 'contrast2_within2', 'contrast2_within3'), sep = ' ') %>%
    mutate_if(is.character, str_trim) %>%
    rename('p_value' = 'p.value',
           't_ratio' = 't.ratio')  %>%
    quick_sig()
  
  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)
  
  return(list(anova_model, anova_summary, pairs_hsd))
  
}

## rm anova one between one within-------------------------------------------------------------------------------------
aov_rm_one_between_one_within <- function(df, dir_output, prefix, var_id, var_dependent, var_between, var_within){
  require(afex)
  
  anova_model <- aov_ez(
    id = var_id, 
    dv = var_dependent, 
    data = df,
    within = var_within,
    between = var_between,
    anova_table=list(correction="none", es = "none"))
  
  
  label_anova <- c("intercept", var_between, var_within, str_c(var_between, '*', var_within)) %>% as_tibble()
  
  anova_summary <-  summary(anova_model)[[4]][] %>% 
    as_tibble() %>%
    bind_cols(label_anova,.) %>%
    rename("effect" = "value",
           "ss" = "Sum Sq",
           "df_num" = "num Df",
           "ss_error" = "Error SS",
           "df_den" = "den Df",
           "f" = "F value",
           "p_value" = "Pr(>F)") 
  
  anova_summary <- anova_summary %>%
    mutate(stat = 'aov_rm_one_between_one_within',
           var_dv   = var_dependent,
           var_id   = var_id,
           var_between = var_between,
           var_within = var_within
    ) %>%
    select(stat, var_dv, var_id, var_between, var_within, everything()) %>%
    quick_sig()
  
  interaction <- paste("~", var_within, '*', var_between, sep = "")
  
  emm <- emmeans(anova_model, formula(interaction))
  
  pairs_hsd <- pairs(emm) %>%
    as_tibble() %>%
    separate(contrast, into = c('contrast1', 'contrast2'), sep = '-') %>%
    mutate_if(is.character, str_trim) %>%
    separate(contrast1, into = c('contrast1_within', 'contrast1_between'), sep = ' ') %>%
    separate(contrast2, into = c('contrast2_within', 'contrast2_between'), sep = ' ') %>%
    mutate_if(is.character, str_trim) %>%
    rename('p_value' = 'p.value',
           't_ratio' = 't.ratio')  %>%
    quick_sig()
  
  save_aov(dir_output, prefix, anova_summary, pairs_hsd, df)
  
  return(list(anova_model, anova_summary, pairs_hsd))
}

