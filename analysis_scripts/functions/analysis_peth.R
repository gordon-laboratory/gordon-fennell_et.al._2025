combine_session_files <- function(process_blocknames, dir_localdata_sessions, fn_suffix){
  
  for(process_blockname in process_blocknames){
    dir_session_process <- str_c(dir_localdata_sessions, process_blockname, '/')
    
    # read in behavior data
    if(fn_suffix %>% str_detect('.csv')){
      df <- read.csv(str_c(dir_session_process, process_blockname, fn_suffix)) 
    } else if(fn_suffix %>% str_detect('.feather')){
      df <- read_feather(str_c(dir_session_process, process_blockname, fn_suffix)) 
    }
    
    if(process_blockname == process_blocknames[1]){ # if first loop
      df_combined <- df
    } else {
      df_combined <- df %>% bind_rows(df_combined,.)
    }
  }
  
  return(df_combined)
}



combine_session_fp <- function(process_blocknames, dir_localdata_sessions){
  
  for(process_blockname in process_blocknames){
    dir_session_process <- str_c(dir_localdata_sessions, process_blockname, '/')
    
    ## read in fp data
    df <- read_feather(str_c(dir_session_process, process_blockname, '_streams_peth_preprocessed.feather')) %>%
      filter(event_id_char %in% c('access_period', 'spout_extended')) %>%
      rename(trial_num = event_number) %>%
      select(any_of(c("blockname", "region", "trial_num", "time_rel", "time", "signal", "delta_signal_poly", "zscore", "zscore_blsub")))
    
    df_trial_info <- read.csv(str_c(dir_session_process, process_blockname, '_data_trial_summary.csv')) %>%
      select(blockname, trial_num, spout, solution, trial_lick, lick_count, lick_ts_first, lick_ts_last) %>%
      mutate(solution_previous = lag(solution))
    
    df <- df %>%
      left_join(df_trial_info, by = join_by(blockname, trial_num))
    
    if(process_blockname == process_blocknames[1]){ # if first loop
      df_combined <- df
    } else {
      df_combined <- df %>% bind_rows(df_combined,.)
    }
  }
  
  return(df_combined)
}

combine_session_fp_bl <- function(process_blocknames, dir_localdata_sessions){
  
  for(process_blockname in process_blocknames){
    dir_session_process <- str_c(dir_localdata_sessions, process_blockname, '/')
    
    ## read in fp data
    df <- read_feather(str_c(dir_session_process, process_blockname, '_streams_baseline_preprocessed.feather')) %>%
      select(any_of(c("blockname", "region", "time", "signal", "delta_signal_poly", "zscore", "zscore_blsub")))
    
    if(process_blockname == process_blocknames[1]){ # if first loop
      df_combined <- df
    } else {
      df_combined <- df %>% bind_rows(df_combined,.)
    }
  }
  
  return(df_combined)
}


get_mean_signals <- function(df){
  df %>%
    summarise(
      delta_signal_poly_sem = sd(delta_signal_poly)/sqrt(n()),
      delta_signal_poly_mean = delta_signal_poly %>% mean(),
      zscore_sem = sd(zscore)/sqrt(n()),
      zscore_mean = zscore %>% mean(),
      zscore_blsub_sem = sd(zscore_blsub)/sqrt(n()),
      zscore_blsub_mean = zscore_blsub %>% mean(), 
      .groups = 'drop'
    )
}


get_peth_binned_summary <- function(df_peth, peth_tm_bins, grouping_vars){
  
  # for each bin defined in peth_tm_bins
  for(bin_id in seq(1,nrow(peth_tm_bins))){
    
    peth_tm_bin <- peth_tm_bins[bin_id,] # get tm_bin info for loop
    
    # filter time_rel to tm_bins, define time_bin,  summarise and combine
    if(bin_id == 1){
      
      df_fp_peth_summary <- df_peth %>%
        filter(time_rel > peth_tm_bin$tm_bin_start, time_rel < peth_tm_bin$tm_bin_end) %>%
        mutate(time_bin = peth_tm_bin$tm_bin_id) %>%
        group_by(across(all_of(grouping_vars))) %>%
        get_mean_signals()
      
    } else {
      df_fp_peth_summary <- df_peth %>%
        filter(time_rel > peth_tm_bin$tm_bin_start, time_rel < peth_tm_bin$tm_bin_end) %>%
        mutate(time_bin = peth_tm_bin$tm_bin_id) %>%
        group_by(across(all_of(grouping_vars))) %>%
        get_mean_signals() %>%
        bind_rows(df_fp_peth_summary,.)
    }
  }
  
  # compute mean signal at baseline of subsequent trial
  if(sum('solution' %in% grouping_vars)){
    grouping_vars_previous <- grouping_vars
    grouping_vars_previous[grouping_vars_previous=='solution'] <- 'solution_previous'
    
    
    if(sum('solution_previous' %in% names(df_fp_peth)) == 1){
      df_fp_peth_summary <- df_peth %>%
        filter(time_rel > -3, time_rel < 0) %>%
        mutate(time_bin = 'baseline_subsequent') %>%
        group_by(across(all_of(grouping_vars_previous))) %>%
        get_mean_signals() %>%
        rename(solution = solution_previous) %>%
        filter(!is.na(solution)) %>%
        bind_rows(df_fp_peth_summary,.)
    }
  }
  
  
  
  
  return(df_fp_peth_summary)
}

