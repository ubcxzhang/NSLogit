getKL <- function(data) {
  res_max <- rep(0, dim(data$count)[1])
  res_mean <- rep(0, dim(data$count)[1])
  res_diff <- rep(0, dim(data$count)[1])
  for (i in 1:dim(data$count)[1]) {
    grp1 <- data$count[i, ][data$condt==unique(data$condt)[1]]
    grp2 <- data$count[i, ][data$condt==unique(data$condt)[2]]
    de1 <- density(grp1, from = 0, to = max(data$count[i, ]))
    de2 <- density(grp2, from = 0, to = max(data$count[i, ]))
    if (all(de1$x==de2$x) == FALSE) {return(0)}
    KL <- kld(de1$y, de2$y)
    res_max[i] <- KL$max
    res_mean[i] <- KL$mean
    res_diff[i] <- sum(abs(de1$y - de2$y))
  }
  list(max = res_max, mean = res_mean, diff = res_diff)
}

getlogfc <- function(data) {
  cells.1 <- data$condt==unique(data$condt)[1]
  cells.2 <- data$condt==unique(data$condt)[2]
  data.1 <- apply(
    X = data$count[, cells.1, drop = FALSE],
    MARGIN = 1,
    FUN = function(x) {
      return(log(x = mean(x = expm1(x = x)) + 1))
    }
  )
  data.2 <- apply(
    X = data$count[, cells.2, drop = FALSE],
    MARGIN = 1,
    FUN = function(x) {
      return(log(x = mean(x = expm1(x = x)) + 1))
    }
  )
  total.diff <- (data.1 - data.2)
}

calkl <- function(y1, y2) {
  kl <- 0
  for (i in 1:length(y1)) {
    if (y2[i] == 0|y1[i] == 0) {next()}
    kl <- kl + y2[i]*log(y2[i]/y1[i])
  }
  abs(kl)
}

kld <- function(y1, y2) {
  list(
    mean = mean(calkl(y1, y2), calkl(y2, y1)),
    max = max(calkl(y1, y2), calkl(y2, y1)))
}

args = commandArgs(trailingOnly=TRUE)
filename <- as.character(args[1])
work_dir <- as.character(args[2]) # the dir containing all files should not end with "/"
work_dir <- paste0(work_dir, "/")

# set R library
.libPaths(paste0(work_dir, "3.5"))

realdata_part_file <- paste0(work_dir, "GSE135893_count/")
sct <- readRDS(paste0(realdata_part_file, filename, ".rds"))

# use the scale in paper for filtering data
sct$count <- log(sct$count + 1)

fc <- getlogfc(sct)
kl <- getKL(sct)
res <- cbind(
  as.numeric(abs(fc)),
  as.numeric(kl$max),
  as.numeric(kl$mean),
  as.numeric(kl$diff))
colnames(res) <- c("fc", "kl_max", "kl_mean", "kl_diff")
rownames(res) <- rownames(sct$count)
res <- as.data.frame(res)
saveRDS(res, file = paste0(realdata_part_file, filename, "_fckl.rds"))
