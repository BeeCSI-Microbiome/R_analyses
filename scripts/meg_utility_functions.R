###############
## Functions ##
###############
## Utility functions that can be optionally used later for producing
## graphs and performing reshaping operations

import("data.table")
import("metagenomeSeq")
import("utils")
import("reshape2")
import("Biobase")
import("stats")
import("limma")
# Call variables from parent scope
`..` <- function (..., .env = sys.parent(2)) {
  get(deparse(substitute(...)), env = .env)
}

# Misc reshape function for data table
melt_dt <- function(D, level_id) { # where is it pulling D, and level_id
  temp <- melt(D, variable.name = "Sample", value.name = "Normalized_Count") # error occurs here. 
  names(temp) <- c("Name", "ID", "Normalized_Count")
  temp <- data.table(cbind(rep(level_id, nrow(temp)), temp))
  names(temp)[1] <- "Level_ID"
  return(temp)
}


data_subset <- function(data_obj, subsets) {
  local_meta <- data.table(pData(data_obj))
  local_subset <- c()          
  for( c in 1:length(subsets) ) {
    
    conditional_terms <- unlist(strsplit(subsets[[c]], " "))
    conditional_string <- paste("local_meta[[\'", conditional_terms[1],
                                "\']] ", conditional_terms[2],
                                " \'", conditional_terms[3], "\'",
                                sep = "", collapse = "")
    if(length(local_subset) > 0) {
      local_subset <- intersect(local_subset, which(eval(parse(text = conditional_string))))
    }
    else {
      local_subset <- which(eval(parse(text = conditional_string)))
    }
  }
  return(data_obj[, local_subset])
}


meg_fitZig <- function(data_list,
                       data_names,
                       metadata,
                       zero_mod,
                       data_mod,
                       filter_min_threshold,
                       contrast_list,
                       random_effect_var,
                       outdir,
                       analysis_name,
                       analysis_subset,
                       data_type,
                       pval = 0.1,
                       top_hits = 200) {
  settings <- zigControl(maxit = 100, verbose = T)
  
  local_obj <- data_list
  res <- list()
  for( l in 1:length(local_obj) ) {
    
    filter_threshold <- quantile(rowSums(MRcounts(local_obj[[l]])), filter_min_threshold)
    if( filter_threshold > filter_min_threshold ) filter_threshold <- filter_min_threshold
    local_obj[[l]] <- local_obj[[l]][which(rowSums(MRcounts(local_obj[[l]])) >= filter_threshold ), ]
    
    if(length(analysis_subset) > 0) {
      local_obj[[l]] <- data_subset(local_obj[[l]], analysis_subset)
    }
    
    
    col_selection <- as.integer(which(colSums(MRcounts(local_obj[[l]]) > 0) >= 1))
    local_obj[[l]] <- local_obj[[l]][, col_selection]
    
    mod_select <- model.matrix(eval(parse(text = data_mod)), data = pData(local_obj[[l]]))
    zero_mod_select <- zero_mod[col_selection, ]
    
    #cumNorm(local_obj[[l]])  # This is a placeholder for metagenomeSeq; we don't actually use these values
    
    tryCatch(
      {
        if( is.na(random_effect_var) ) {
          res[[l]] <- fitZig(obj = local_obj[[l]],
                             mod = mod_select,
                             zeroMod = zero_mod_select,
                             control = settings,
                             useCSSoffset = F)
        }
        else {
          res[[l]] <- fitZig(obj = local_obj[[l]],
                             mod = mod_select,
                             zeroMod = zero_mod_select,
                             control = settings,
                             useCSSoffset = F,
                             useMixedModel = T,
                             block = pData(local_obj[[l]])[, random_effect_var])
        }
      },
      error = function(e) {
        print(paste("Encountered an error performing fitZig for", data_type, " ", data_names[l], " ", analysis_name,
                    sep = "", collapse = ""))
      },
      finally = {
        if( length(res) != l ) {
          next
        }
      }
    )
    
    local_contrasts <- contrast_list
    local_contrasts[[length(local_contrasts)+1]] <- res[[l]]@fit$design
    # $ S4 error caused here
    
    names(local_contrasts)[length(local_contrasts)] <- "levels"
    
    contrast_matrix <- do.call(makeContrasts, local_contrasts)
    colnames(contrast_matrix) <- make.names(contrast_list)
    
    contrast_fit <- contrasts.fit(res[[l]]@fit, contrast_matrix)        
    contrast_fit <- eBayes(contrast_fit)
    
    stats_results <- data.table(
      Node.Name = character(),
      Contrast = character(),
      logFC = numeric(),
      CI.L = numeric(),
      CI.R = numeric(),
      AveExpr = numeric(),
      t = numeric(),
      P.Value = numeric(),
      adj.P.Val = numeric(),
      B = numeric()
    )
    
    for( c in 1:ncol(contrast_fit$contrasts) ) {
      tophits <- topTable(contrast_fit, p.value = pval, confint = T,
                          number = top_hits, sort.by = "AveExpr", coef = c)
      
      if( nrow(tophits) > 0) {
        temp_res <- data.table(
          Node.Name = rownames(tophits),
          Contrast = rep(colnames(contrast_fit$contrasts)[c], nrow(tophits))
        )
        temp_res <- cbind(temp_res, tophits)
        stats_results <- rbind(stats_results, temp_res)
      }
      else {
        print(paste("No significant results for", data_type,
                    data_names[l], analysis_name,
                    colnames(contrast_fit$contrasts)[c],
                    sep = " ", collapse = ""))
      }
    }
    
    if( nrow(stats_results) > 0 ) {
      write.csv(stats_results,
                file = paste(outdir, "/", analysis_name, "_", data_type, "_",
                           data_names[l], "_",
                           contrast_list[1], "_Model_Contrasts.csv",
                           sep = "", collapse = ""),
                quote = F, row.names = F)
    }
  }
}
