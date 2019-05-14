#######################
# specificity methods #
#######################
#' Calculate specificity for a given numeric data.frame
#'
#' Function to calculate specificity for a given numeric data.frame and return
#' a data.frame with original value with specificity index of 9 methhods.
#'
#' @param df data.frame of numeric.
#' @param cutoff Numeric. Values under cutoff will set to 0(unexpressed) in
#' Specificity method Counts.
#' @return data.frame with original value with specificity index of 9 methhods.
ts_index <- function(df,
                     cutoff = 1) {
  len <- 1:ncol(df)
  df$Tau <- apply(df[, len], 1, ts_Tau)
  df$Gini <- apply(df[, len], 1, ts_Gini)
  df$Tsi <- apply(df[, len], 1, ts_Tsi)
  df$Counts <- apply(df[, len], 1, function(x) {x <- ts_Counts(x, cutoff)})
  df$Hg <- apply(df[, len], 1, ts_Hg)
  df$Zscore <- ts_Z(df[, len])
  df$Spm <- apply(df[, len], 1, ts_Spm)
  df$Ee <- ts_Ee(df[, len])
  df$Pem <- ts_Pem(df[, len])

  df$Mean <- apply(df[, len], 1, expr_mean)
  df$Max <- apply(df[, len], 1, expr_max)
  return(df)
}



# N is the number of tissues in the data set, n is the number of tissues
# expressed; X is the maximal specificity value for a certain gene among all
# tissues, maxX is maximal values.

#' Tau
#'
#' Function require a vector with expression of one gene in different tissues.
#' If expression for one tissue is not known, gene specificity for this gene is
#' NA. Minimum 2 tissues.
#'
#' @param x A numeric vector.
#' @examples
#' \dontrun{
#' ts_Tau(c(1,2,3,2,8))
#' }
ts_Tau <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(max(x)!=0)
      {
        x <- (1-(x/max(x)))
        res <- sum(x, na.rm=TRUE)
        res <- res/(length(x)-1)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    }
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}



#' Gini
#'
#' Function require a vector with expression of one gene in different tissues.
#' If expression for one tissue is not known, gene specificity for this gene is
#' NA. Transformation: x*(N/(N-1)).
#'
#' @param x A numeric vector.
#' @importFrom DescTools Gini
#' @examples
#' \dontrun{
#' ts_Gini(c(1,2,3,2,8))
#' }
ts_Gini <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(sum(x!=0))
      {
        res <- DescTools::Gini(x)*(length(x)/(length(x)-1))
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    }
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}



#' Tsi
#'
#' Function require a vector with expression of one gene in different tissues.
#' If expression for one tissue is not known, gene specificity for this gene is
#' NA.
#'
#' @param x A numeric vector.
#' @examples
#' \dontrun{
#' ts_Tsi(c(1,2,3,2,8))
#' }
ts_Tsi <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(sum(x!=0))
      {
        res <- max(x) / sum(x)
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    }
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}



#' Counts
#'
#' Function require a vector with expression of one gene in different tissues.
#' If expression for one tissue is not known, gene specificity for this gene is
#' NA. Function requires setting of a \code{cutoff} (rpkm/FPKM/TPM).
#' Transformation: (1-x/N)*(N/(N-1)).
#'
#' @param x A numeric vector.
#' @param cutoff numeric. values under \code{cutoff} will be set to
#' 0(unexpression).
#'
#' @examples
#' \dontrun{
#' ts_Counts(c(0,1,2,3,2,8), cutoff = 1)
#' }
ts_Counts <- function(x, cutoff)
{
  if(all(!is.na(x)))
  {
    res <- length(which(x > cutoff))
    if (res > 0)
    {
      res <- (1 - res/length(x))*(length(x)/(length(x)-1))  #Modification: To bring to normalized scale
    }
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}



#' Ee
#'
#' Function require a data frame with expression data, and give back a vector
#' with EEi values for each gene. If expression for one tissue is not known,
#' gene specificity for this gene is NA. Transformation: x/maxX.
#'
#' @param x data.frame.
#' @examples
#' \dontrun{
#' ts_Ee(tmp_tpm)
#' }
ts_Ee <- function(x)
{
  if(!all(is.na(x)))
  {
    x <- as.matrix(x)
    x[x<0] <- NA
    x <- cbind(x, r=rowSums(x, na.rm=FALSE))
    x <- rbind(x, c=colSums(x, na.rm=TRUE))
    x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <- x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] / (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x), which(x[nrow(x),]>0)] / x[nrow(x), ncol(x)])

    res <- apply(x[-nrow(x),-ncol(x)], c(1), FUN=max)
    res <- res/max(res, na.rm=TRUE) #Modification: To bring to normalized scale
  } else {
    res <- NA
    print("No data avalable.")
  }
  return(res)
}



#' Pem
#'
#' Function require a data frame with expression data, and give back a vector
#' with PEM scores. If expression for one tissue is not known, gene specificity
#' for this gene is NA. Transformation: x/(maxX).
#'
#' @param x data.frame.
#' @examples
#' \dontrun{
#' ts_Ee(tmp_tpm)
#' }
ts_Pem <- function(x)
{
  if(!all(is.na(x)))
  {
    x <- as.matrix(x)
    x[x<0] <- NA
    x <- cbind(x, r=rowSums(x, na.rm=FALSE)) #Add column with expression of gene per tissue
    x <- rbind(x, c=colSums(x, na.rm=TRUE))	#Add row with expression of all genes in a given tissue
    x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] <- x[which(x[,ncol(x)]!=0), which(x[nrow(x),]!=0)] / (x[which(x[,ncol(x)]>0), ncol(x)] %o% x[nrow(x), which(x[nrow(x),]>0)] / x[nrow(x), ncol(x)]) #calculate the score

    x[x<1] <- 1
    x <- log10(x)

    x<- abs(x)
    res <- apply(x[-nrow(x),-ncol(x)], c(1), FUN=max) #choose only the maximal score for each gene
    res <- res/max(res, na.rm=TRUE) #Modification: To bring to normalized scale from 0 to 1
  } else {
    res <- NA
    print("No data avalable.")
  }
  return(res)
}



#' Hg, entropy
#'
#' Function require a vector with expression of one gene in different tissues.
#' If expression for one tissue is not known, gene specificity for this gene is
#' NA. Transformation: 1-x/log2N.
#'
#' @param x A numeric vector.
#'
#' @examples
#' \dontrun{
#' ts_Hg(c(0,1,2,3,2,8))
#' }
ts_Hg <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(sum(x) !=0)
      {
        p <- x / sum(x)
        res <- -sum(p * log2(p), na.rm=TRUE)
        res <- 1 - (res / log2(length(p))) #Modification: To bring to normalized scale
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    }
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}



#' Z-score
#'
#' Function require a vector with expression of one gene in different tissues.
#' If expression for one tissue is not known, gene specificity for this gene is
#' NA. Transformation: x/n-1/sqrt(N).
#'
#' @param x data.frame.
#'
#' @examples
#' \dontrun{
#' ts_Z(c(0,1,2,3,2,8))
#' }
ts_Z <- function(x)
{
  if(all(!is.na(x)))
  {
    res <-  apply(scale(t(x), center=TRUE, scale=TRUE),2,max)/((length(x[1,])-1)/sqrt(length(x[1,])))
    res[is.na(res)] <- 0
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}



#' SPM score from TISGED
#'
#' Function require a vector with expression of one gene in different tissues.
#' If expression for one tissue is not known, gene specificity for this gene is
#' NA. Transformation: maxX.
#'
#' @param x A numeric vector.
#'
#' @examples
#' \dontrun{
#' ts_Spm(c(0,1,2,3,2,8))
#' }
ts_Spm <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(sum(x) !=0)
      {
        spm <- x^2/as.vector(x%*%x)
        res <- max(spm) #Modification:To bring to normalized scale. Choose max
      } else {
        res <- 0
      }
    } else {
      res <- NA
      #print("Expression values have to be positive!")
    }
  } else {
    res <- NA
    #print("No data for this gene avalable.")
  }
  return(res)
}


