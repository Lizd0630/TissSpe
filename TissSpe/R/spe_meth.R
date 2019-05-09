# Function require a vector with expression of one gene in different tissues.
# If expression for one tissue is not known, gene specificity for this gene is NA
# Minimum 2 tissues
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


# Function require a vector with expression of one gene in different tissues.
# If expression for one tissue is not known, gene specificity for this gene is NA
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


#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
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


#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
#Function requires setting of a treshold (PSI)
ts_Counts <- function(x, min)
{
  if(all(!is.na(x)))
  {
    res <- length(which(x > min))
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


#Function require a data frame with expression data, and give back a vector with EEi values for each gene
#If expression for one tissue is not known, gene specificity for this gene is NA
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


#Function require a data frame with expression data, and give back a vector with PEM scores
#If expression for one tissue is not known, gene specificity for this gene is NA
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


#Hg entropy
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
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


#Z-score
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
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


#SPM score from TISGED
#Function require a vector with expression of one gene in different tissues.
#If expression for one tissue is not known, gene specificity for this gene is NA
ts_Spm <- function(x)
{
  if(all(!is.na(x)))
  {
    if(min(x, na.rm=TRUE) >= 0)
    {
      if(sum(x) !=0)
      {
        spm <- x^2/(x%*%x)
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

