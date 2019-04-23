# Library (to use fitdistr)
library(MASS)

# Read all the files with signals corresponding to characters
file1 <- read.csv('crbmfjie')
file2 <- read.csv('fjekphqr')
file3 <- read.csv('joismlqr')
file4 <- read.csv('kuvypjaq')
file5 <- read.csv('mjpntxsi')

# Apply absolute value to get the distribution
file1$sample <- abs(file1$sample)
file2$sample <- abs(file2$sample)
file3$sample <- abs(file3$sample)
file4$sample <- abs(file4$sample)
file5$sample <- abs(file5$sample)

# Read the secret phrase
secret <- read.csv('secret.csv')
secret$sample <- abs(secret$sample)

#######################
# Find right N
#######################

# Plot one example to get fixed number of samples for one character
plot(abs(file5$sample))
file_len <- as.numeric(as.character(length(file1$sample)))
N <- file_len/8


#######################
# Find known letters
#######################

# Matrix to put in the letters
# Dimension of each letter's parameter (both shapes and rates) to the maximum appearance frequency (5)
shapes <- double(26)
rates <- double(26)

mapping <- data.frame(letters, shapes, rates, row.names=letters)

# Letters in files
file1_char <- c("c", "r", "b", "m", "f", "j", "i", "e")
file2_char <- c("f", "j", "e", "k", "p", "h", "q", "r")
file3_char <- c("j", "o", "i", "s", "m", "l", "q", "r")
file4_char <- c("k", "u", "v", "y", "p", "j", "a", "q")
file5_char <- c("m", "j", "p", "n", "t", "x", "s", "i")

# Build matrix of known characrers and vector of files
ds <- matrix(c(file1_char, file2_char, file3_char, file4_char, file5_char), nrow=5, ncol=8, byrow=TRUE)
files <- c(file1, file2, file3, file4, file5)

# To avoid errors by limiting the optimizer of 'fitdistr' to c(0,0)
eps <- sqrt(.Machine$double.eps)

# Iterate over the various files
for (i in 1:5) {
  # Iterate over letters in each file
  for (j in 1:8){
    # Use indexes to correctly slice the data
    ind1 <- 10000*(j-1) + 1
    ind2 <- 10000*j
    samp <- files[i]$sample[ind1:ind2]
    
    # Condition to avoid applying mean on first value insertion (cause vector has 0 in it)
    # Fit distribution: lower is used to avoid warning (sometimes optimizer searches in negative space)
    if (mapping[ds[i, j],2]==0) {
      mapping[ds[i, j],2] = fitdistr(samp, densfun ="gamma", lower=c(eps,eps))$estimate[[1]]
    } else {
      mapping[ds[i, j],2] = mean(fitdistr(samp, densfun ="gamma", lower=c(eps,eps))$estimate[[1]])
    }
    
    if (mapping[ds[i, j],3]==0) {
      mapping[ds[i, j],3] = fitdistr(samp, densfun ="gamma", lower=c(eps,eps))$estimate[[2]]
    } else {
      mapping[ds[i, j],3] = mean(fitdistr(samp, densfun ="gamma", lower=c(eps,eps))$estimate[[2]])
    }
  }
}


#######################
# Decryptate message
#######################
# Variable storing the decrypted message
message <- c()

# Iterate over the representations of the letters
for (i in 1:49){
  # Instantiate indexes and a control variable (checks wether a letter can be decripted)
  ind1 <- 10000*(i-1) + 1
  ind2 <- 10000*i
  control <- FALSE
  
  # Find the parameters for the sample in secret message
  tmp_shape <- fitdistr(secret$sample[ind1:ind2], densfun="gamma", lower=c(eps,eps))$estimate[[1]]
  tmp_rate <- fitdistr(secret$sample[ind1:ind2], densfun="gamma", lower=c(eps,eps))$estimate[[2]]
  
  # Iterate over the discovered letters in the map and check with empirical threshold if the parameters are compatible
  for (l in letters){
    shape_diff <- abs(mapping[l, 2] - tmp_shape)
    rate_diff <- abs(mapping[l, 3] - tmp_rate)
    if (shape_diff <= 0.2){
      if (rate_diff <= 0.03){
        message <- c(message, l)
        control <- TRUE
        break
      }
    }
  }
  
  # If no letter is found append an asterisk
  if (!control){
    message <- c(message, "*")
  }
}

print(message)


#####################################
# Inferring letters from the message
#####################################
# Letter in position 7 is a 'w'
mapping['w', 2] <- fitdistr(secret$sample[60001:70000], densfun="gamma", lower=c(eps,eps))$estimate[[1]]
mapping['w', 3] <- fitdistr(secret$sample[60001:70000], densfun="gamma", lower=c(eps,eps))$estimate[[2]]

# Letter in position 35 is a 'z'
mapping['z', 2] <- fitdistr(secret$sample[340001:350000], densfun="gamma", lower=c(eps,eps))$estimate[[1]]
mapping['z', 3] <- fitdistr(secret$sample[340001:350000], densfun="gamma", lower=c(eps,eps))$estimate[[2]]

# Letter in position 49 is a 'd'
mapping['d', 2] <- fitdistr(secret$sample[480001:490000], densfun="gamma", lower=c(eps,eps))$estimate[[1]]
mapping['d', 3] <- fitdistr(secret$sample[480001:490000], densfun="gamma", lower=c(eps,eps))$estimate[[2]]

