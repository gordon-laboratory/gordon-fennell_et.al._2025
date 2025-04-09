
get_peaks_moving_zscore <- function(streams_baseline, window_size, minpeakheight, minpeakdistance) {
  # window_size, minpeakheight, and minpeakdistance are all in samples
  
  suppressWarnings({
    # Find peaks in moving_zscore
    peaks_info <- findpeaks(streams_baseline$moving_zscore, 
                            minpeakheight = minpeakheight, 
                            minpeakdistance = minpeakdistance)
    
    # If no peaks found, return an empty tibble
    if (is.null(peaks_info)) {
      return(streams_baseline %>% head(0))
    }
    
    # Extract the peak times (indices) and values
    peak_times <- streams_baseline[peaks_info[,2],]  # 2nd column gives peak indices
    
    # Return a tibble with region, peak times, and values
    return(peak_times)
  })
}


get_peaks_zscore_from_peaks_moving_zscore <- function(df_peaks_moving_zscore, streams_baseline, sample_range){
  
  df_peaks_temp <- df_peaks_moving_zscore[1,]  %>%
    select(blockname, region, time) 
  
  df_peaks_zscore <- streams_baseline %>%
    filter(region == df_peaks_temp$region) %>%
    filter(time >= df_peaks_temp$time) %>%
    arrange(time) %>%
    mutate(sample_number = row_number()) %>%
    filter(sample_number <= sample_range) %>%
    filter(zscore == max(zscore)) %>%
    filter(time == min(time)) %>%
    select(-sample_number)
  
  return(df_peaks_zscore)
}

get_peak_zscore_peth <- function(df_peaks_zscore, streams_baseline, time_window){
  
  df_peaks_temp <- df_peaks_zscore[1,]  %>%
    select(blockname, region, time, peak_number) 
  
  tm_bins <- seq(time_window[1], time_window[2], 0.05) # generate time bins and set labels to end of each window
  
  df_peak_zscore_peth <- streams_baseline %>%
    filter(region == df_peaks_temp$region) %>%
    mutate(time_rel = time - as.numeric(df_peaks_temp$time)) %>%
    mutate(time_rel = time_rel %>% round(2)) %>%
    filter(
      time_rel > time_window[1],
      time_rel < time_window[2]
    ) %>%
    mutate(time_rel = cut(time_rel, tm_bins, tm_bins[2:length(tm_bins)])) %>%
    mutate(time_rel = time_rel %>% as.character() %>% as.double()) %>%           # convert time to numeric
    filter(time_rel < max(time_rel)) %>%
    mutate(peak_number = df_peaks_temp$peak_number)
  
  return(df_peak_zscore_peth)
}

# Function to fit the exponential decay model and extract time constant (tau)
calculate_time_constant <- function(df) {
  # Fit an exponential decay model: zscore = y0 * exp(-time_rel / tau)
  fit <- tryCatch(
    nls(zscore ~ y0 * exp(-time_rel / tau), 
        data = df, 
        start = list(y0 = max(df$zscore), tau = 1),
        control = nls.control(maxiter = 500)),
    error = function(e) NULL  # Return NULL in case of a fitting error
  )
  
  # If the fit is successful, extract the time constant (tau), else return NA
  if (!is.null(fit)) {
    tau_estimate <- coef(fit)["tau"]
  } else {
    tau_estimate <- NA
  }
  
  return(tibble(tau_estimate = tau_estimate ))
}

compute_acf <- function(df) {
  # Extract the time series data for the region
  ts_data <- df$zscore
  
  # Compute autocorrelation without plotting
  acf_values <- acf(ts_data, lag.max = 30*10, plot = FALSE)
  
  # Extract autocorrelation values and corresponding lags
  acf_df <- data.frame(
    lag = acf_values$lag,
    acf = acf_values$acf
  )
  
  return(acf_df)
}

get_streams_baseline_peaks <- function(dir_localdata_sessions, window_size, minpeakheight, minpeakdistance, overwrite){
  require(zoo)
  require(pracma)
  
  
  print('obtaining baseline peaks')
  
  dir_sessions <- list.dirs(dir_localdata_sessions, recursive = F) 
  #dir_sessions <- dir_sessions[dir_sessions %>% str_detect('acd01')] 
  
  # for each session 
  for(dir_session in dir_sessions){
    
    
    loop_blockname <- sub(".*/", "", dir_session) 
    
    fn_output <- str_c(dir_session, '/', loop_blockname, '_streams_baseline_preprocessed_peaks.csv')
    
    if(file.exists(fn_output) & overwrite == 0){
      print(str_c(' ~ ', loop_blockname, '- skipped'))
      
    }else{
      print(str_c(' ~ ', loop_blockname))
      
      
      streams_baseline <- read_feather(str_c(dir_session, '/', loop_blockname, '_streams_baseline_preprocessed.feather')) %>%
        mutate(time = (time + 0.001) %>% round_down_to_nearest_05())
      
      
      # compute moving z-score used for detecting spontaneous peaks
      streams_baseline <- streams_baseline %>%
        mutate(
          rolling_mean = rollapply(signal, window_size, mean, align = "right", fill = NA),
          rolling_sd = rollapply(signal, window_size, sd, align = "right", fill = NA),
          moving_zscore = (signal - rolling_mean) / rolling_sd
        ) %>%
        filter(!is.na(moving_zscore)) %>%
        mutate(time = time %>% round(2))
      
      return(streams_baseline)
      
      # get peaks from moving z-score data
      df_peaks <- streams_baseline %>%
        group_by(region) %>%
        group_modify(~ get_peaks_moving_zscore(.x, window_size, minpeakheight, minpeakdistance)) %>%
        ungroup() %>%
        mutate(peak_number = row_number()) %>%
        mutate(time = time %>% round(2))
      
      # use peaks from moving z-score to find corresponding peaks in z-score
      df_peaks_zscore <- do.call(
        rbind, apply(df_peaks, 1, function(row) {
          get_peaks_zscore_from_peaks_moving_zscore(
            as.data.frame(t(row)),
            streams_baseline,
            minpeakdistance/2)
        }
        )
      ) %>%
        ungroup() %>%
        mutate(time = time %>% round(2)) %>%
        mutate(peak_number = row_number())
      
      
      # return peth relative to peak times
      df_peaks_zscore_peth <- do.call(
        rbind, apply(df_peaks_zscore, 1, function(row) {
          get_peak_zscore_peth(
            as.data.frame(t(row)), 
            streams_baseline %>% ungroup(), 
            c(-6, 9)
          )
        }
        )
      ) %>%
        mutate(time = time %>% round(2))%>%
        mutate(time_rel = time_rel %>% round(2))
      
      # get peak rate for each region
      df_peakcount_by_region <- df_peaks_zscore %>%
        group_by(blockname, region) %>%
        summarise(peak_count = n(), .groups = 'drop') %>%
        mutate(peak_rate_per_min = peak_count / 
                 ((streams_baseline$time %>% max() - streams_baseline$time %>% min())/60)
        ) 
      
      # return basleine correlation 
      streams_baseline_wide <- streams_baseline %>%
        select(blockname, time, region, zscore) %>%
        spread(region, zscore)
      
      df_corr_mat  <- tidy_cormatrix(streams_baseline_wide, unique(streams_baseline$region), 'blockname')
      
      # save data
      df_peaks_zscore %>% write.csv(str_c(dir_session, '/', loop_blockname, '_streams_baseline_preprocessed_peaks.csv'))
      df_peaks_zscore_peth %>% write.csv(str_c(dir_session, '/', loop_blockname, '_streams_baseline_preprocessed_peaks_peth.csv'))
      df_peakcount_by_region %>% write.csv(str_c(dir_session, '/', loop_blockname, '_streams_baseline_preprocessed_peakcount.csv'))
      df_corr_mat %>% write.csv(str_c(dir_session, '/', loop_blockname, '_streams_baseline_preprocessed_corr.csv'))
      
    }
  } 
}



combined_baseline_analysis <- function(dir_localdata_sessions, dir_baseline_summary_output){
  session_blocknames <- list.files(dir_localdata_sessions) 
  session_subjects <- str_sub(session_blocknames, 12, 16) %>% unique()
  
  # read in quality metrics for fiber x session exclusion
  fns_qc_metrics <- list.files(dir_localdata_sessions, pattern = 'quality_metrics.csv', recursive = T, full.names = T)
  combined_qc_metrics <- bind_rows(lapply(fns_qc_metrics, read.csv))
  
  for(loop_subject in session_subjects){
   
    
    # get acf for each subject x region across all sessions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    file_paths <- list.files(dir_localdata_sessions, pattern = str_c(loop_subject, '_streams_baseline_preprocessed.feather'), recursive = T, full.names = T)
   
    # skip subject if 
    if(length(file_paths) == 0){
      print(str_c(' ~ ', loop_subject, ' - skipped'))
      next
    }
    
    print(str_c(' ~ ', loop_subject))
    
    combind_streams_baseline <- file_paths %>%
      map_dfr(read_feather)  
    
    # filter out files with poor signal
    combind_streams_baseline <- combind_streams_baseline %>%
      left_join(combined_qc_metrics %>% select(blockname, region, exclude_poor_signal), by = c('blockname', 'region')) %>%
      filter(exclude_poor_signal == 0 | is.na(exclude_poor_signal)) %>% # is.na for lha ratio
      select(-exclude_poor_signal)
    
    df_acf_by_region <- combind_streams_baseline %>%
      group_by(region) %>%
      group_modify(~ compute_acf(.x)) %>%
      ungroup() 
    
    df_acf_by_region <- df_acf_by_region %>%
      mutate(session_count = length(unique(combind_streams_baseline$blockname)))  %>%
      mutate(subject = loop_subject) %>%
      select(subject, region, session_count, everything())
    
    df_acf_by_region %>% write.csv(str_c(dir_baseline_summary_output, loop_subject, '_allsessions_acf.csv'))
    
    
    # get mean spontaneous spike info for each subject x region across all sessions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # compute and summarise tau
    file_paths <- list.files(dir_localdata_sessions, pattern = str_c(loop_subject, '_streams_baseline_preprocessed_peaks_peth.csv'), recursive = T, full.names = T)
    
    combind_peaks_peth <- file_paths %>% map_dfr(read.csv)
    
    # filter out files with poor signal
    combind_peaks_peth <- combind_peaks_peth %>%
      left_join(combined_qc_metrics %>% select(blockname, region, exclude_poor_signal), by = c('blockname', 'region')) %>%
      filter(exclude_poor_signal == 0 | is.na(exclude_poor_signal)) %>% # is.na for lha ratio
      select(-exclude_poor_signal)
    
    df_tau_by_region <- combind_peaks_peth %>%
      filter(time_rel >= 0) %>% # filter to decay
      group_by(region, time_rel) %>%
      summarise(zscore = mean(zscore), .groups = 'drop') %>%
      group_by(region) %>%
      do(calculate_time_constant(.)) %>%
      mutate(session_count = length(unique(combind_peaks_peth$blockname))) %>%
      mutate(subject = loop_subject) %>%
      select(subject, region, session_count, everything())
    
    # combine and summarise peak count per min
    file_paths <- list.files(dir_localdata_sessions, pattern = str_c(loop_subject, '_streams_baseline_preprocessed_peakcount.csv'), recursive = T, full.names = T)
    
    combined_peakcount <- file_paths %>%
      map_dfr(read.csv)
    
    # filter out files with poor signal
    combined_peakcount <- combined_peakcount %>%
      left_join(combined_qc_metrics %>% select(blockname, region, exclude_poor_signal), by = c('blockname', 'region')) %>%
      filter(exclude_poor_signal == 0 | is.na(exclude_poor_signal)) %>% # is.na for lha ratio
      select(-exclude_poor_signal)
    
    combined_peakcount <- combined_peakcount %>%
      mutate(subject = loop_subject) %>%
      group_by(subject, region) %>%
      summarise(peak_rate_per_min = mean(peak_rate_per_min), .groups = 'drop')
    
    # write combined summary
    df_tau_by_region %>%
      left_join(combined_peakcount, by = join_by(subject, region)) %>%
      write.csv(str_c(dir_baseline_summary_output, loop_subject, '_allsessions_spon_summary.csv'))
    
  }
}
