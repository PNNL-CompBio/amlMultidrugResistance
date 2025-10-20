
library(dplyr)


## group_var is the group variable in the pData of the msnset. This is like 'response' in rf_modeling
## pred_group is like pred.cls. It is the level of the group variable that is 
## identified as '1', the other is '0'. 

## Notes: This function test the LOOCV performance of a logistic model based on a logratio of markers.
## To use, it requires: 1) An msnset (logscale) and 2) A variable of interest (group_var) with exactly 2 levels (case vs control).
## The output is a list with the following items:
## 1) The probabilities for each held out sample, 
## 2) The features selected in each LOO model
## 3) A count of frequency for each feature chosen in the LOO models
## 4) The AUC performance of the LOO model outputs
## 5) The LOO model prediction for each sample
## 6) Another list, containing 
##    6.1) The model trained on the FULL data (not LOO)
##    6.2) The features selected by the FULL model

## Note: pred_group should be the level which corresponds to the 'case' condition. 
## control_group should be the level corresponding to the 'control' condition. If NA then 'group_var' must be a variable with two levels
## To select markers, the function runs a t.test (same variance) and filters according to pval_cutoff. 
## If method is 'simple', then the most significant features are chosen (according to t.test pvalue). N_markers determines the number chosen.
## If method is 'kmeans', then markers are chosen after K-means clustering of the covariance matrix of the features below pval_cutoff. 
##    (This is meant to prioritize features with independent information)

logistic_model <- function(msnset, group_var, pred_group, N_markers, sample_pooling = NA, control_group = NA,
                           pval_cutoff = 0.01, method = "kmeans", var.equal = TRUE){
   m1 <- msnset
   meta <- pData(m1) %>% dplyr::select(!!group_var) %>% mutate(sample_name_unique = rownames(.)) %>%
      dplyr::select(group_ = !!group_var, sample = sample_name_unique) 
   group_lvls <- unique(meta$group_)
   
   if(length(group_lvls) != 2 & is.na(control_group)){
      print("Only implemented for comparisons of two classes.")
      break
   } else if (is.na(control_group)) {
      control_group = setdiff(unique(meta$group_), pred_group)
   }
   other_groups = setdiff(unique(meta$group_), c(pred_group, control_group))
   ## Ensures limma_gen computes LogFC as case_group - control_group
   meta <- meta %>%
      mutate(group_ = factor(group_, levels = c(control_group, pred_group, other_groups)))
   m1$group_ <- meta$group_
   
   all_samples <- meta %>% filter(group_ %in% c(control_group, pred_group)) %>% pull(sample)
   loocv_probs = data.frame(sample = meta$sample, group_ = meta$group_) %>%
      filter(sample %in% all_samples) %>%
      mutate(group_ = case_when(group_ == pred_group ~ 1,
                                group_ == control_group ~ 0))
   rownames(loocv_probs) <- loocv_probs$sample
   feature_df <- data.frame()
   m <- m1
   print(all_samples)
   for (loo_sample in all_samples){
      print(loo_sample)
      if (is.data.frame(sample_pooling)){
         ## Find the patient of the loo_sample
         patient <- sample_pooling %>%
            filter(sample_name == loo_sample) %>% pull(patient_name)
         ## Find the samples from that patient
         exclude_samples <- sample_pooling %>%
            filter(patient_name == patient) %>%
            pull(sample_name)
         ## Exclude the other samples from that patient, since they're not independent from loo_sample
         exclude_samples <- setdiff(exclude_samples, loo_sample)
         m1 <- m[, setdiff(sampleNames(m), exclude_samples)]
      }
      samples_keep = sampleNames(m1) != loo_sample
      m1_ <- m1[, samples_keep]
      
      ###############################################
      # choose features to maximize signal-to-noise #
      ###############################################
      
      ## The choice of features is based on all BUT the loo_sample
      model_prep <- logistic_model_helper1(m1_, pred_group, control_group, pval_cutoff, N_markers, loo_sample, method, var.equal)
      new_feature_df = model_prep[[3]]
      case_marker = model_prep[[1]]
      control_marker = model_prep[[2]]
      feature_df <- rbind(feature_df, new_feature_df)
      
      logistic_df <- logistic_model_helper2(m1, case_marker, control_marker, TRUE, pred_group)
      
      logistic_df <- logistic_df %>%
         ## excluding samples from other levels, different from pred_group, control_group
         filter(sample %in% all_samples)
      
      logistic_df_test <- logistic_df %>% filter(sample == loo_sample) %>% dplyr::select(-sample, -group_)
      logistic_df_train <- logistic_df %>% filter(sample != loo_sample) %>% dplyr::select(-sample)
      model_logistic = suppressWarnings(glm(group_ ~ ., data = logistic_df_train, family = "binomial"))
      loo_logistic_pred = predict.glm(model_logistic, logistic_df_test, type = "response")
      loocv_probs[loo_sample, "loocv_logistic_prob"] <- loo_logistic_pred[[1]]
   }
   
   predProb = loocv_probs[, "loocv_logistic_prob"]
   pred = ROCR::prediction(predProb, loocv_probs[, "group_"])
   
   feature_summary <- feature_df %>% group_by(loo_sample) %>% 
      summarize(features = list(feature))
   selected <- feature_summary$features
   names(selected) <- feature_summary$loo_sample
   
   model_prep <- logistic_model_helper1(m1, pred_group, control_group, pval_cutoff, N_markers, "full_model", method, var.equal)
   new_feature_df = model_prep[[3]]
   case_marker = model_prep[[1]]
   control_marker = model_prep[[2]]
   cov_mats = model_prep[[4]]
   
   logistic_df <- logistic_model_helper2(m1, case_marker, control_marker, TRUE, pred_group)
   
   logistic_df <- logistic_df %>%
      ## excluding samples from other levels, different from pred_group, control_group
      filter(sample %in% all_samples)
   logistic_df_train <- logistic_df %>% dplyr::select(-sample)
   full_model = suppressWarnings(glm(group_ ~ ., data = logistic_df_train, family = "binomial"))
   
   return(list(prob = predProb, features = selected, top = rev(sort(table(unlist(selected)))), 
               auc = ROCR::performance(pred, "auc")@y.values[[1]], pred = pred, 
               full_model = list('model' = full_model, 'features' = new_feature_df),
               #cov_mats = cov_mats
               ))
}




logistic_model_helper1 <- function(m1_, pred_group, control_group, pval_cutoff, N_markers, loo_sample, method, var.equal){
   if (var.equal){
      test <- diffexp_helper(m1_, "group_", paste(pred_group, control_group, sep = "-")) %>%
         select(feature, logFC, P.Value = t_test_pval)
   } else {
      test <- diffexp_helper(m1_, "group_", paste(pred_group, control_group, sep = "-")) %>%
         select(feature, logFC, P.Value = welch_pval)
   }
   
   if ((length(sampleNames(m1_)) < 2*N_markers + 1) && !combine_markers){
      print("Reducing number of markers to have full rank in model regression")
      N_markers = floor(length(sampleNames(m1_))/2) - 1
   }
   case_marker <- test %>% filter(logFC > 0) %>% arrange(P.Value) %>% filter(P.Value < pval_cutoff) 
   case_marker = logistic_model_helper3(m1_, case_marker, paste(pred_group, control_group, sep = "-"), N_markers, var.equal, method)
   case_cov <- case_marker[[2]]
   case_marker <- case_marker[[1]]
   control_marker <- test %>% filter(logFC < 0) %>% arrange(P.Value) %>% filter(P.Value < pval_cutoff) 
   control_marker = logistic_model_helper3(m1_, control_marker, paste(pred_group, control_group, sep = "-"), N_markers, var.equal, method)
   control_cov <- control_marker[[2]]
   control_marker <- control_marker[[1]]
   case_df = data.frame()
   control_df = data.frame()
   try(case_df <- data.frame(feature = case_marker$feature, group_ = pred_group, loo_sample = loo_sample, 
                             s2n = case_marker$group_s2n[[1]], silhouette_width = case_marker$silhouette_width))
   try(control_df <- data.frame(feature = control_marker$feature, group_ = control_group, loo_sample = loo_sample, 
                                s2n = control_marker$group_s2n[[1]], silhouette_width = control_marker$silhouette_width))
   feature_df <- rbind(case_df, control_df)
   
   return(list(case_marker$feature, control_marker$feature, feature_df, list(case_cov, control_cov)))
}


logistic_model_helper2 <- function(m1, case_marker, control_marker, combine_markers, pred_group){
   case_mat = t(exprs(m1)[case_marker, , drop = FALSE]) %>% as.data.frame()
   control_mat = t(exprs(m1)[control_marker, , drop = FALSE]) %>% as.data.frame()
   if (combine_markers){
      if (ncol(control_mat) > 0){
         control_mat$control_avg = rowSums(control_mat)/ncol(control_mat)
         # control_mat$control_avg = apply(control_mat, 1)
         control_mat <- control_mat %>% select(control_avg)
      }
      if (ncol(case_mat) > 0){
         case_mat$case_avg = rowSums(case_mat)/ncol(case_mat)
         # case_mat$case_avg = apply(case_mat, 1, mean)
         case_mat <- case_mat %>% select(case_avg)
      }
   }
   
   logistic_df <- data.frame(sample = sampleNames(m1), group_ = m1$group_ == pred_group) %>% 
      cbind(case_mat, control_mat) 
   return(logistic_df)
}


logistic_model_helper3 <- function(m1_, test, contrast, N_markers, var.equal, method){
   groups_ = strsplit(contrast, "-")[[1]]
   compute_s2n <- function(markers){
      marker_mat = exprs(m1_)[markers, , drop = FALSE]
      new_row = data.frame(value = colSums(marker_mat)/nrow(marker_mat), 
                           sample = colnames(marker_mat),
                           group_ = m1_$group_)

      t_test = t.test(value ~ group_, data = new_row %>% filter(group_ %in% groups_), 
                      alternative = "two.sided", var.equal = var.equal)
      return(abs(t_test[[1]])[[1]])
   }
   track_df = data.frame()
   N = min(nrow(test), 5)
   N_markers = min(N_markers, nrow(test))
   if (method == "simple"){
      markers = test$feature[1:N_markers]
      s2n = compute_s2n(markers)
      xx = data.frame(feature = markers, s2n = s2n, group_s2n = s2n, silhouette_width = NA,
                      N_pastcutoff = nrow(test))
      cov_mat = NA
   } else if (method == "s2n_max_v1"){
      for (i in 1:N){
         potential_markers = test$feature[[i]]
         track_df = rbind(track_df,
                          data.frame(feature = potential_markers, s2n = NA, start = i))
         for (j in 2:N_markers){
            other_features = setdiff(test$feature, potential_markers)
            next_df = data.frame()
            for (feature in other_features){
               next_df <- rbind(next_df,
                                data.frame(feature = feature, 
                                           s2n = compute_s2n(c(potential_markers, feature)), 
                                           start = i))
               
            }
            next_df <- next_df %>% arrange(-s2n) %>% slice(1)
            potential_markers = c(potential_markers, next_df$feature[[1]])
            track_df = rbind(track_df, next_df)
         }
      }
      xx = track_df %>% group_by(start) %>% mutate(group_s2n = s2n[length(s2n)]) %>%
         arrange(-group_s2n) %>% ungroup() %>% slice(1:N_markers)
   } else if (method == "kmeans"){
      N_markers = min(N_markers, nrow(test)-1)
      cov_mat = t(exprs(m1_)[test$feature, m1_$group_ %in% groups_]) %>% as.data.frame()
      cov_mat = sweep(cov_mat, 2, apply(cov_mat, 2, sd), FUN = '/')
      cov_mat = cov(cov_mat)
      
      set.seed(123)
      km.out = kmeans(cov_mat, centers = N_markers, nstart = 200)
      test$cluster = km.out$cluster
      x = cluster::silhouette(test$cluster, dist(cov_mat))
      sil_width = summary(x)[[4]]
      markers = test %>% group_by(cluster) %>%
         slice(1) %>% pull(feature)
      s2n = compute_s2n(markers)
      xx = data.frame(feature = markers, s2n = s2n, group_s2n = s2n, silhouette_width = sil_width,
                      N_pastcutoff = nrow(test))
   }
   return(list(xx, cov_mat))
}










