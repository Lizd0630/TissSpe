# Function requires data frame with expression values
# Mean values between replicares are calculated
# namelist must be the unique word-start identifier to recognize sample replicates
# identifier is the colname of gene names or AS events
rep_mean <- function(df, name_vect, identifier) {
  mat <- matrix(NA, ncol = length(name_vect) + 1, nrow = nrow(df))
  colnames(mat) <- c(identifier, name_vect)
  for (i in name_vect) {
    mat[, i] <- rowMeans(df[,regexpr(paste0("^", i), colnames(x)) > 0],
                         na.rm = TRUE,
                         dim = 1)
  }
  rownames(mat) <- rownames(df)
  df_new <- as.data.frame(mat)
  df_new[, identifier] <- df[, identifier]
  return(df_new)
}


#Function requires data frame to be normalized
#1. All 0 are set to NA, to exclude them from quantile normalization
#2. Data are quantile normalized
#3. 0 values (the one set to NA) are set back to 0
quant_norm <- function(df) {
  df[df == 0] <- NA
  mat <- as.matrix(df)
  mat <- normalize.quantiles(mat)
  mat[is.na(mat)] <- 0
  return(data.frame(mat))
}

# Function to cut data.frame range into n equal-width interval points,
# max point and min point.
# this function is design for processing data, in which detecting differantial
# expression rely on diffenrence, like PSI(delta psi).
# require a data.frame as input, n medium intervals
seq_cut <- function(df,
                    n = 10) {
  region <- range(as.matrix(df), na.rm = T)
  bks <- quantile(region[1]:region[2], probs = seq(0, 1, 1/n))
  bks[1] <- bks[1] + 1e-6
  bks[n] <- bks[n] - 1e-6
  bks <- c(region[1], bks, region[2])

  df_rank <- t(apply(df, 1, function(i) {as.numeric(as.vector(cut(i, breaks = bks, labels = c(0:(n+1)), include.lowest = T)))}))
  colnames(df_rank) <- colnames(df)
  df_rank <- as.data.frame(df_rank)
  return(df_rank)
}

# Function to cut data.frame range into n equal-density interval points,
# over-represent interval point and under-represent interval point
# require a data.frame as input, n medium intervals
quant_cut <- function(df,
                      min = 0.2,
                      max = 2,
                      n = 10) {
  vect <- as.vector(t(df))
  region <- c(min(vect, na.rm = T), max(vect, na.rm = T))
  vect <- vect[(vect >= min & vect <= max)]
  bks <- as.vector(quantile(vect, probs = seq(0, 1, 1/n)))
  bks[1] <- bks[1] - 1e-6
  bks[n+1] <- bks[n+1] + 1e-6
  bks <- c(region[1], bks, region[2])

  df_rank <- t(apply(df, 1, function(i) {as.numeric(as.vector(cut(i, breaks = bks, labels = c(0:(n+1)), include.lowest=T)))}))
  colnames(df_rank) <- colnames(df)
  df_rank <- as.data.frame(df_rank)
  return(df_rank)
}


# Function require a vector with expression of one gene in different tissues.
# Mean is calculated taking in account tissues with 0 expression: 2+0+4=2
expr_mean <- function(x) {
  if(!all(is.na(x)))
  {
    res <- mean(x, na.rm = T)
  } else {
    res <- NA
  }
  return(res)
}



# Function require a vector with expression of one gene in different tissues.
# Max is calculated taking in account tissues with 0 expression: 2+0+4=4
expr_max <- function(x)
{
  if(!all(is.na(x)))
  {
    res <- max(x, na.rm = T)
  } else {
    res <- NA
  }
  return(res)
}
