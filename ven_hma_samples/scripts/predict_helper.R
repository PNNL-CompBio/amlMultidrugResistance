


## Function to compute signal value using the markers contained in 'full_model'. This is an element in the output produced by
## the 'logistic_model' function, in the script 's2n_model' which is loaded at the top.
## NOTE: The markers in the 'full_model' objects are intersected with the set of COMPLETE (NO NA) features from the new data.
predict_helper <- function(new_msnset, full_model, case_group, control_group){
   case_marker = full_model$features %>% filter(group_ == case_group) %>% pull(feature)
   control_marker = full_model$features %>% filter(group_ == control_group) %>% pull(feature)
   
   ## This ensures we don't have any NA when computing the signal in the new samples.
   case_marker = intersect(case_marker, complete_features) %>% intersect(., featureNames(new_msnset))
   control_marker = intersect(control_marker, complete_features) %>% intersect(., featureNames(new_msnset))
   
   case_mat = t(exprs(new_msnset)[case_marker, , drop = FALSE]) %>% as.data.frame()
   control_mat = t(exprs(new_msnset)[control_marker, , drop = FALSE]) %>% as.data.frame()
   if (ncol(control_mat) > 0){
      control_mat$control_avg = rowMeans(control_mat, na.rm = T)
      control_mat <- control_mat %>% select(control_avg)
   }
   if (ncol(case_mat) > 0){
      case_mat$case_avg = rowMeans(case_mat, na.rm = T)
      case_mat <- case_mat %>% select(case_avg)
   }
   
   logistic_df <- data.frame(sample = sampleNames(new_msnset)) %>% 
      cbind(case_mat, control_mat) 
   logistic_pred = predict.glm(full_model$model, logistic_df, type = "response")
   coefs = full_model$model$coefficients
   logistic_df <- logistic_df %>%
      mutate(logistic_signal = coefs[['(Intercept)']] + coefs[['case_avg']]*case_avg + coefs[['control_avg']]*control_avg,
             logistic_signal2 = case_avg - control_avg,
             prob = logistic_pred,
             prob2 = 1/(1+exp(-logistic_signal)))
   return(logistic_df)
}



## Function to compute signal value using two sets of markers. The signal will be mean(case_markers) - mean(control_markers)
## NOTE: The markers are intersected with the set of COMPLETE (NO NA) features from the new data.
predict_helper_ <- function(new_msnset, case_markers, control_markers){
   complete_features = featureNames(new_msnset)[(rowSums(is.na(exprs(new_msnset))) == 0)]
   case_markers = intersect(case_markers, complete_features)
   control_markers = intersect(control_markers, complete_features)
   case_mat = t(exprs(new_msnset)[case_markers, , drop = FALSE]) %>% as.data.frame()
   control_mat = t(exprs(new_msnset)[control_markers, , drop = FALSE]) %>% as.data.frame()
   if (ncol(control_mat) > 0){
      control_mat$control_avg = rowMeans(control_mat, na.rm = T)
      control_mat <- control_mat %>% select(control_avg)
   }
   if (ncol(case_mat) > 0){
      case_mat$case_avg = rowMeans(case_mat, na.rm = T)
      case_mat <- case_mat %>% select(case_avg)
   }
   logistic_df <- data.frame(sample = sampleNames(new_msnset)) %>% 
      cbind(case_mat, control_mat) %>%
      mutate(signal = case_avg - control_avg)
   return(logistic_df)
}


