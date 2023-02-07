source("functions.R")

# 1. Generate different source populations 
xx <- seq(from = 0, to = 1, by = 0.001)
mm <- expand.grid(rep(list(xx), 2))
mm <- mm[(mm[,1] + mm[,2] < 1) & (mm[,1] > 0) & (mm[,2] > 0),]
(nn <- nrow(mm))
mm[,3] <- round(runif(nn, min = -5, max = 5) , 4)
mm[,4] <- round(runif(nn, min = -5, max = 5) , 4)
mm[,5] <- round(runif(nn, min = -5, max = 5) , 4)
mm[,6] <- round(runif(nn, min = -5, max = 5) , 4)
colnames(mm) <- c("a", "b", "b0", "b1", "b2", "b3")


# 2. for matching factors X and H, obtain two distributions of (B | H) and cfb 
triples <-  matrix(NA, nrow = nn, ncol = 9)
cfb <- matrix(NA, nrow = nn, ncol = 8)
for(i in 1:nn){
  a <-  mm[i,1]
  b <-  mm[i,2]
  coeff <- as.numeric(mm[i,3:6])
  
  # x = 0, x = 1, x = 2
  px <- c(a, b, 1 - a - b) 
  # h = -1, h = 1
  ph <- c(a + b, 1 - a - b) 
  # x | h = -1
  px_h1 <- c(a/(a + b), b/(a + b), 0)
  # x | h = 1
  px_h3 <- c(0, 0, 1)
  
  xx <- c(0,1,2)
  
  f0 <- coeff[1] + coeff[2]*xx
  py0 <- expit(f0)
  
  f1 <- (coeff[1] + coeff[3]) + (coeff[2] + coeff[4])*xx
  py1 <- expit(f1)
  
  out_h1x <- out_h3x <- matrix(NA, nrow = 3, ncol = 3)
  out_h1h <- out_h3h <- matrix(NA, nrow = 3*3, ncol = 3)
  
  # 1) matching factor X
  for(r in 1:3){
    out_h1x[r,] <- px_h1[r]*c(py1[r]*(1 - py0[r]),
                              py1[r]*py0[r] + (1 - py1[r])*(1 - py0[r]),
                              (1 - py1[r])* py0[r])
    
    out_h3x[r,] <- px_h3[r]*c(py1[r]*(1 - py0[r]),
                              py1[r]*py0[r] + (1 - py1[r])*(1 - py0[r]),
                              (1 - py1[r])* py0[r])
  }
  
  out_h1x <- colSums(out_h1x)
  out_h3x <- colSums(out_h3x)
  triple1 <- c(out_h1x, out_h3x)
  names(triple1) <- c("p1", "p2", "p3", "q1", "q2", "q3")
  
  
  # 2) matching factor H
  for(k in 1:3){
    for(j in 1:3){
      out_h1h[(k-1)*3+j,] <- px_h1[k]*px_h1[j]*c(py1[k]*(1 - py0[j]),
                                                 py1[k]*py0[j] + (1 - py1[k])*(1 - py0[j]),
                                                 (1 - py1[k])* py0[j])
      
      out_h3h[(k-1)*3+j,] <- px_h3[k]*px_h3[j]*c(py1[k]*(1 - py0[j]),
                                                 py1[k]*py0[j] + (1 - py1[k])*(1 - py0[j]),
                                                 (1 - py1[k])* py0[j])
    }
  }
  
  out_h1h <- colSums(out_h1h)
  out_h3h <- colSums(out_h3h)
  triple2 <- c(out_h1h, out_h3h)
  names(triple2) <- c("p1", "p2", "p3", "q1", "q2", "q3")
  
  # 3) cfb calculations
  cfb[i, 1] <- cal_cfb(ph[2], triple1[1], triple1[3], triple1[4], triple1[6])[[2]]
  cfb[i, 2] <- cal_cfb(ph[2], triple2[1], triple2[3], triple2[4], triple2[6])[[2]]
  cfb[i, 3:8] <- c(a, b, coeff)
  triples[i, ] <- c(ph[2], triple1[1], triple1[3], triple1[4], triple1[6],
                    triple2[1], triple2[3], triple2[4], triple2[6])
}