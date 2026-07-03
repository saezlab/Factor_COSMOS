#'\code{make_df_and_targets_great_again}
#'
#'This function Check wether sample names are coherent between the measurment dataframe and the target dataframe. If they are coherent,
#' the target dataframe rows are reordered to match the column order of the measurment dataframe.
#'@param df the measurment n*m dataframe (n is number of omic features, m is number of samples) where columns are ordered by conditions.
#'@param targets A n*2 dataframe, where n is the number of samples. First column correspond to samples, second column correspond to conditions.
#'@return A list, first element is the cleaned measurment dataframe and second element is the cleaned target dataframe.
make_df_and_targets_great_again <- function(df, targets)
{
  names(targets)[c(1,2)] <- c("sample","condition")
  bad_samples <- c(NA)
  i <- 1
  for (sample in names(df))
  {
    if (!(sample %in% targets[,1]))
    {
      bad_samples <- c(bad_samples, sample)
    }
    else
    {
      if (i == 1)
      {
        clean_targets <- as.data.frame(matrix(NA, 1, length(targets[1,])))
        names(clean_targets) <- names(targets)
        clean_targets[1,] <- targets[targets$sample == sample,]
        i <- i+1
      }
      else
      {
        clean_targets <- as.data.frame(rbind(clean_targets,targets[targets$sample == sample,]))
        i <- i+1
      }
    }
  }
  targets <- clean_targets
  if (length(bad_samples) > 1)
  {
    bad_samples <- bad_samples[-1]
    print(paste("These samples were not found : ", bad_samples, sep = ""))
    print(bad_samples)
    df <- df[,!(names(df) %in% bad_samples)]
  }
  return(list(df,targets))
}


#'\code{makeContrastsAlt}
#'
#'This function create a contrast matrix to be used by limma.
#'
#'@param targets A n*2 dataframe, where n is the number of samples. First column correspond to samples, second column correspond to conditions.
#'@param comparisons a list of numeric vectors. Each vector represent which condition should be conpared. Example :
#'c(2,-1) means that the first condition should be substracted from second condition. Vectors can be more than two element for complex contrasts.
#'@return a contrast matrix
makeContrastsAlt <- function(targets, comparisons)
{
  cont.matrix <- matrix(0,nrow = length(unique(targets$condition)), ncol = length(comparisons))
  i <- 1
  for (comparison in comparisons)
  {
    for (j in 1:length(comparison))
    {
      cont.matrix[abs(comparison[j]),i] <- cont.matrix[abs(comparison[j]),i]+(comparison[j]/abs(comparison[j]))
    }
    i <- i + 1
  }
  return(cont.matrix)
}

#'\code{checkInputs}
#'
#'This function makes sure that the input of runLimma are properly formatted
#'
#'@param measurments the measurment n*m dataframe (n is number of omic features, m is number of samples) where columns are ordered by conditions.
#'@param targets A n*2 dataframe, where n is the number of samples. First column correspond to samples, second column correspond to conditions.
#'@return TRUE if all is good, FALSE otherwise
checkInputs <- function(measurments, targets)
{
  if(class(measurments) != "data.frame")
  {
    error_message <- paste("The measurments argument should be a data.frame. It's currently a", paste(class(measurments), ".",sep = ""))
    return(list(FALSE, error_message))
  }
  else
  {
    if(dim(measurments)[1] == 0)
    {
      error_message <- "The measurments dataframe doesn't seem to contain any measurments..."
      return(list(FALSE, error_message))
    }
    else
    {
      if(dim(measurments)[2] == 0)
      {
        error_message <- "The measurments dataframe doesn't seem to contain any samples..."
        return(list(FALSE, error_message))
      }
      else
      {
        if(class(as.matrix(measurments)[,1]) != "numeric")
        {
          return(list(FALSE, "The measurments dataframe should contain only numerical values (or NAs)."))
        }
        else
        {
          if(class(targets) != "data.frame")
          {
            error_message <- paste("The targets argument should be a data.frame. It's currently a", paste(class(targets), ".",sep = ""))
            return(list(FALSE, error_message))
          }
          else
          {
            if(dim(targets)[2] < 2)
            {
              return(list(FALSE,"The targets dataframe should have at least two columns, sample names and conditions."))
            }
            else
            {
              if(dim(targets)[1] != dim(measurments)[2])
              {
                error_message <- paste("The targets dataframe should have as many samples (targets rows) as the measurements (measurments columns). Currently, the targets dataframe has", paste(dim(targets)[1], "samples and the measurements have", paste(dim(measurments)[2],"samples.")))
                return(list(FALSE, error_message))
              }
              else
              {
                #placeholder in case i think of more stuff to check
              }
            }
          }
        }
      }
    }
  }
  return(list(TRUE, "All seems to be in order..."))
}

#'\code{runLimma}
#'
#'This function is a wrapper of limma made to facilitate the use of limma differential analysis
#'
#'@param measurments the measurment n*m dataframe (n is number of omic features, m is number of samples) where columns are ordered by conditions.
#'@param targets A n*2 dataframe, where n is the number of samples. First column correspond to samples, second column correspond to conditions.
#'@param comparisons a list of numeric vectors. Each vector represent which condition should be conpared. Example :
#'c(2,-1) means that the first condition should be substracted from second condition. Vectors can be more than two element for complex contrasts.
#'@param regress_out in case the user which to exclude possible confounding factors from the analysis, the user can provide additional columns in the targets dataframe.
#'then, the confounding factor can be regressed out by indicating the number of the column of the target dataframe describing it. Only one factor can be regressed out at the present time.
#'@return a list. First element is the limma model fitted with the contrast matrix, this is the usual output of limma. Second element is the contrast matrix that was used. third element is the fitted limma object without contrasts.
runLimma <- function(measurements, targets, comparisons = NULL, regress_out = NULL)
{
  input_check <- checkInputs(measurements, targets)
  if (input_check[[1]]) #input has correct format
  {
    measurements <- measurements[,targets[,1]]
    if (!is.null(comparisons))
    {
      if (!is.null(regress_out))
      {
        for (regressor in regress_out)
        {
          measurements <- removeBatchEffect(measurements, targets[,regressor])
        }
      }
      
      cont.matrix <- makeContrastsAlt(targets, comparisons)
      
      cont.matrix <- as.data.frame(cont.matrix)
      row.names(cont.matrix) <- unique(targets$condition)
      cont.matrix <- as.matrix(cont.matrix)
      
      fcond <- factor(targets$condition, levels = unique(targets$condition))
      
      design <- model.matrix(~0+fcond)
      design <- as.data.frame(design)
      names(design) <- unique(targets$condition)
      design <- as.matrix(design)
      View(design)
      print(cont.matrix)
      View(measurements)
      
      fit <- lmFit(measurements, design)
      fit2 <- contrasts.fit(fit, cont.matrix)
      fit2 <- eBayes(fit2)
      
      return(list(fit2, cont.matrix, fit))
    }
  }
  else
  {
    print(input_check[[2]])
    return(input_check[[1]])
  }
}

#'\code{ttopFormatter}
#'
#'This function is simply designed to format the toptable of limma with first column as gene identifiers instead of only row.names.
#'
#'@param ttop a toptable dataframe generated by the topTable function of limma
#'
#'@return a dataframe similar to the output of topTable function of limma but with first column as IDs instead of onyl row.names.
ttopFormatter <- function(ttop)
{
  ttop$ID <- row.names(ttop)
  ttop <- ttop[,c(7,1,2,3,4,5,6)]
  ttop <- ttop[complete.cases(ttop),]
  return(ttop)
}


#' ttop_list_to_t_table
#'
#' ipsum...
#'
#' @param ttop_list  ipsum...
#' @return ipsum...
#' @export
ttop_list_to_t_table <- function(ttop_list)
{
  if(length(ttop_list) > 1)
  {
    t_table <- merge(ttop_list[[1]][,c(1,4)], ttop_list[[2]][,c(1,4)], by = "ID", all = T)
    if(length(ttop_list) > 2)
    {
      for(i in 3:length(ttop_list))
      {
        t_table <- merge(t_table, ttop_list[[i]][,c(1,4)], by = "ID", all = T)
      }
    }
  }
  else
  {
    t_table <- ttop_list[[1]][,c(1,4)]
  }
  names(t_table) <- c("ID",names(ttop_list))
  return(t_table)
}

#' limma_res_to_ttop_list
#'
#' ipsum...
#'
#' @param limma_res  ipsum...
#' @param comp_names  ipsum...
#' @param number  ipsum...
#' @param adjust.method  ipsum...
#' @return ipsum...
#' @export
limma_res_to_ttop_list <- function(limma_res, comp_names, number, adjust.method = "fdr")
{
  ttop_list <- list()
  n_comp <- length(limma_res[[2]][1,])
  for(i in 1:n_comp)
  {
    ttop_list[[i]] <- ttopFormatter(topTable(limma_res[[1]], coef = i, number = number, adjust.method = adjust.method))
    ttop_list[[i]] <- ttop_list[[i]][complete.cases(ttop_list[[i]]),]
  }
  names(ttop_list) <- comp_names
  return(ttop_list)
}


# Define a custom function to find the most extreme value
most_extreme <- function(x) {
  x <- na.omit(x) # Remove NAs
  max_x <- max(x) # Find maximum value
  min_x <- min(x) # Find minimum value
  
  if (abs(min_x) > abs(max_x)) {
    return(min_x)
  } else {
    return(max_x)
  }
}

summarise_t_table_by_most_extrem <- function(t_table)
{
  row.names(t_table) <- t_table$ID
  t_table <- t_table[,-1]
  t_table[is.na(t_table)] <- 0
  
  t_table$ID <- row.names(t_table)
  # Group the z-score data frame by 'psite' and calculate the most extreme values for each group
  t_table$ID <- gsub("___.*","",t_table$ID)
  
  
  t_table <- t_table %>% group_by(ID) %>% summarise_each(funs(most_extreme))
  
  # Convert the grouped z-score data back to a data frame and set row names using the 'psite' column
  t_table <- as.data.frame(t_table)
  row.names(t_table) <- t_table$ID
  t_table <- t_table[,-1]
  
  t_table[t_table == 0] <- NA
  
  t_table$ID <- row.names(t_table)
  
  t_table <- t_table[,c(length(t_table[1,]),1:(length(t_table[1,]) - 1))]
  
  
  
  
  return(t_table)
}