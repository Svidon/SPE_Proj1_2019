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

# To avoid errors by limiting the optimizer of 'fitdistr' to c(0,0)
# It still does not work perfectly with c(0,0) so a really small value is used (epsilon)
eps <- sqrt(.Machine$double.eps)


################################
# FUNCTIONS
################################

#################################################################################
# Function that gets known letters
# params: "ds" = matrix containing the known letters per file
#         "files" = vector of all the files needed
# returns: matrix "tmp_mapping" with the parameters corresponding to each letter
#################################################################################

get_mapping <- function(ds, files){
  
  # Matrix to put in the letters
  shapes <- double(26)
  rates <- double(26)
  tmp_mapping <- data.frame(letters, shapes, rates, row.names=letters)
  
  # Iterate over the various files
  for (i in 1:length(files)) {
    # Iterate over letters in each file
    for (j in 1:ncol(ds)){
      # Use indexes to correctly slice the data (get the N samples representing one letter)
      ind1 <- N*(j-1) + 1
      ind2 <- N*j
      samp <- files[i]$sample[ind1:ind2]
      
      # Fit distribution: lower is used to avoid warning (sometimes optimizer searches in negative space)
      # If condition used to avoid using mean when a letter is first encountered
      if (tmp_mapping[ds[i, j],2] == 0.0){
        tmp_mapping[ds[i, j],2] <- fitdistr(samp, densfun ="gamma", lower=c(eps,eps))$estimate[[1]]
      } else {
        tmp_mapping[ds[i, j],2] <- mean(c(fitdistr(samp, densfun ="gamma", lower=c(eps,eps))$estimate[[1]], tmp_mapping[ds[i, j],2]))
      }
      
      if (tmp_mapping[ds[i, j],3] == 0.0){
        tmp_mapping[ds[i, j],3] = fitdistr(samp, densfun ="gamma", lower=c(eps,eps))$estimate[[2]]
      } else {
        tmp_mapping[ds[i, j],3] = mean(c(fitdistr(samp, densfun ="gamma", lower=c(eps,eps))$estimate[[2]], tmp_mapping[ds[i, j],3]))
      }
    }
  }
  
  return(tmp_mapping)
}


############################################################
# Function that decodes the message
# params: "sec" = secret message
#         "map" = the mapping character to gamma parameters
# returns: vector "tmp_message" with the decoded message
############################################################

decrypt <- function(sec, map){
  # Variable storing the decrypted message
  tmp_message <- c()
  
  # Iterate over the representations of the letters in the secret message
  for (i in 1:(length(sec$sample)/N)){
    # Instantiate indexes and a control variable (checks wether a letter can be decripted)
    ind1 <- N*(i-1) + 1
    ind2 <- N*i
    control <- FALSE
    
    # Find the parameters for the sample in secret message
    tmp_shape <- fitdistr(sec$sample[ind1:ind2], densfun="gamma", lower=c(eps,eps))$estimate[[1]]
    tmp_rate <- fitdistr(sec$sample[ind1:ind2], densfun="gamma", lower=c(eps,eps))$estimate[[2]]
    
    # Iterate over the discovered letters in the map and check with empirical threshold if the parameters are compatible
    for (l in letters){
      shape_diff <- abs(map[l, 2] - tmp_shape)
      rate_diff <- abs(map[l, 3] - tmp_rate)
      if (shape_diff <= 0.3){
        if (rate_diff <= 0.04){
          tmp_message <- c(tmp_message, l)
          control <- TRUE
          break
        }
      }
    }
    
    # If no letter is found append an asterisk
    if (!control){
      tmp_message <- c(tmp_message, "*")
    }
  }
  
  return(tmp_message)
}



##########################################
# Apply functions to solve the assignment
##########################################


#######################
# Find right N
#######################

# Since it's a constant and there are 8 letters and 80k samples in a file it's 10k
file_len <- as.numeric(as.character(length(file1$sample)))
N <- file_len/8


#######################
# Find known letters
#######################

# Letters in files
file1_char <- c("c", "r", "b", "m", "f", "j", "i", "e")
file2_char <- c("f", "j", "e", "k", "p", "h", "q", "r")
file3_char <- c("j", "o", "i", "s", "m", "l", "q", "r")
file4_char <- c("k", "u", "v", "y", "p", "j", "a", "q")
file5_char <- c("m", "j", "p", "n", "t", "x", "s", "i")

# Build matrix of known characrers and vector of files
known_characters <- matrix(c(file1_char, file2_char, file3_char, file4_char, file5_char), nrow=5, ncol=8, byrow=TRUE)
files_list <- c(file1, file2, file3, file4, file5)

mapping <- get_mapping(known_characters, files_list)


#######################
# Decryptate message
#######################
message <- decrypt(secret, mapping)


###########################################
# Inferring other letters from the message
###########################################

# Letter in position 7 is a 'w'
mapping['w', 2] <- fitdistr(secret$sample[60001:70000], densfun="gamma", lower=c(eps,eps))$estimate[[1]]
mapping['w', 3] <- fitdistr(secret$sample[60001:70000], densfun="gamma", lower=c(eps,eps))$estimate[[2]]

# Letter in position 35 is a 'z'
mapping['z', 2] <- fitdistr(secret$sample[340001:350000], densfun="gamma", lower=c(eps,eps))$estimate[[1]]
mapping['z', 3] <- fitdistr(secret$sample[340001:350000], densfun="gamma", lower=c(eps,eps))$estimate[[2]]

# Letter in position 49 is a 'd'
mapping['d', 2] <- fitdistr(secret$sample[480001:490000], densfun="gamma", lower=c(eps,eps))$estimate[[1]]
mapping['d', 3] <- fitdistr(secret$sample[480001:490000], densfun="gamma", lower=c(eps,eps))$estimate[[2]]


#########
# Plots
#########

# Plot of file5, to visually see that the constant N is 10000
plot(abs(file5$sample), xlab="Index", ylab="Value", cex.lab=1.5, cex.axis=1.5)

# Plot the two example letters to prove that even if one parameter is close (shape),
# as the other (rate) is very different the resulting curves are very different
# Generate two samples with the refined parameters obtained
a_sample <- rgamma(N, mapping['a', 2], mapping['a', 3])
b_sample <- rgamma(N, mapping['b', 2], mapping['b', 3])
l_sample <- rgamma(N, mapping['l', 2], mapping['l', 3])
r_sample <- rgamma(N, mapping['r', 2], mapping['r', 3])

# Plot them in the same figure
plot(density(l_sample), col = "red",
     xlab="Value", ylab="Probability", main=NA,
     cex.lab=1.5, cex.axis=1.5, lwd=2)
lines(density(b_sample), col = "blue", lwd=2)
lines(density(a_sample), col = "green", lwd=2)
lines(density(r_sample), col = "purple", lwd=2)
legend(6, 2.0, legend=c("Letter l", "Letter b", "Letter a", "Letter r"), col=c("red", "blue", "green", "purple"), lty=1:2, cex=1.5)


############################################
# Obtain csv to get latex table with script
############################################

# Transform the matrix into a dataframe and change column names accordingly
mapping_df <- as.data.frame(mapping)
names(mapping_df) <- c("letter", "alpha", "lambda")

# Reverse rate and get the parameter scale, which is the lambda in the task
mapping_df$lambda <- 1/mapping_df$lambda

# Drop unknown letters' rows
mapping_df <- mapping_df[mapping_df$alpha>0, ]

# Write to csv the example table with only 2 characters
write.csv(mapping_df[mapping_df$letter == "a" || mapping_df$letter == "b" || mapping_df$letter == "l"
                     || mapping_df$letter == "r", ], file="example.csv")

# Write to csv the complete table
write.csv(mapping_df, file="mapping.csv")
