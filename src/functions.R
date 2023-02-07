# 1. Calculated cfb^* for binary X and Y
# p(H = 0) = 1 - c (or P(X = 0) = 1 - c)
# p(H = 1) = c (or P(X = 1)=c)
# p(B = -1 | X = 0) = p1
# p(B = 0 | X = 0) = p2
# p(B = 1 | X = 0) = p3


# p(B = -1 | X = 1) = q1
# p(B = 0 | X = 1) = q2
# p(B = 1 | X = 1) = q3


# @param p1,p3,q1,q3
# @param c
# @output cfb and joint distribution of (B1, B2, H1, H2)
cal_cfb <- function(c, p1, p3, q1, q3){
  # P(X1 = X2)
  xe <- c^2 + (1-c)^2
  # P(X1 > X2) = P(X1 < X2)
  xl <- xg <- c*(1-c)
  
  # p(B = 0 | X = 0) = p2
  p2 <- 1 - p1 - p3
  # p(B = 0 | X = 1) = q2
  q2 <- 1 - q1 - q3
  
  # p(B = -1)
  b1 <- p1*(1-c) + q1*c
  # p(B = 0) 
  b2 <- p2*(1-c) + q2*c
  # p(B = 1)
  b3 <- p3*(1-c) + q3*c
  
  # p(B1 = B2)
  be <- b1^2 + b2^2 + b3^2
  # p(B1 > B2) = p(B1 < B2)
  bl <- bg <- b1*b2 + b1*b3 + b2*b3
  
  # P(B1 = B2, X1 = X2)
  bexe <- (p1*(1-c))^2 + (p2*(1-c))^2 + (p3*(1-c))^2 + (q1*c)^2 + (q2*c)^2 + (q3*c)^2
  # P(B1 = B2, X1 > X2) =  P(B1 = B2, X1 < X2)
  bexl <- bexg <- (p1*q1 + p2*q2 + p3*q3)*(1-c)*c
  be <- c(bexe, bexg, bexl)
  
  # P(B1 > B2, X1 = X2)
  bgxe <-  p3*p2*(1-c)^2 + p3*p1*(1-c)^2 + p2*p1*(1-c)^2 + q3*q2*c^2 + q3*q1*c^2 + q2*q1*c^2 
  # P(B1 > B2, X1 > X2)
  bgxg <- (q3*p2 + q3*p1 + q2*p1)*(1-c)*c
  # P(B1 > B2, X1 < X2)
  bgxl <- (p3*q2 + p3*q1 + p2*q1)*(1-c)*c
  bg <- c(bgxe, bgxg, bgxl)
  
  # P(B1 < B2, X1 = X2)
  blxe <-  p1*p2*(1-c)^2 + p1*p3*(1-c)^2 + p2*p3*(1-c)^2 + q1*q2*c^2 + q1*q3*c^2 + q2*q3*c^2 
  # P(B1 < B2, X1 > X2)
  blxg <- (q1*p2 + q1*p3 + q2*p3)*(1-c)*c
  # P(B1 < B2, X1 < X2)
  blxl <- (p1*q2 + p1*q3 + p2*q3)*(1-c)*c
  bl <- c(blxe, blxg, blxl)
  
  cfb <- (0.5*bgxe + bgxg)/sum(bg)
  out <- list()
  table <- data.frame(be, bg, bl)
  colnames(table) <- c("B1 = B2", "B1 > B2", "B1 < B2")
  row.names(table) <- c("H1 = H2", "H1 > H2", "H1 < H2")
  
  out[[1]] <- table
  out[[2]] <- cfb
  return(out)
}


# example

cal_cfb(1/4, 80/765, 184/765, 8/585, 64/585)
cal_cfb(1/4, (8/135 + 4/459 + 8/225 + 4/765), (16/135 + 32/459 + 8/225 + 16/765), 8/585, 64/585)

# 2. Calculate the probability of a outcome for X with three levels based on logistic regressions
# h(x) = x^2 - x - 1
# x = {0, 1, 2}
# h = {-1, 1}
# py0_x = Pr(y(0) = 1 | x),
# py1_x = Pr(y(1) = 1 | x)
# log(py0_x/(1 - py0_x)) = b0 + b1x
# log(py1_x/(1 - py1_x)) = b0 + b1x + b2A + b3xA

expit <- function(x) 1/(1 + exp(-x))

# @import expit(x)
# @param a: treatment assignment group
# @param coeff: (b0, b1, b2, b3)
# @param x: c(0,1,2)
# @output Pr(y(a) = 1 | x)
py_x <- function(a, coeff){
  x = c(0,1,2)
  if(a == 0){
    f <- coeff[1] + coeff[2]*x
    p <- expit(f)
  }else{
    f <- (coeff[1] + coeff[3]) + (coeff[2] + coeff[4])*x
    p <- expit(f)
  }
  return(p)
}

# example
coeff <- c(log(2),log(2),log(2),log(2))
py_x(1, coeff)


# @import expit(x)
# @import py_x(a, coeff)
# @param coeff: (b0, b1, b2, b3)
# @param match_by: {"x", "h"}
# @output c(p1,p2,p3,q1,q2,q3) for (B | H)
pb_h <- function(coeff, px_h1, px_h3, match_by){
  py0 <- py_x(a = 0, coeff)
  py1 <- py_x(a = 1, coeff)
  if(match_by == "x"){
    out_h1 <- out_h3 <- matrix(NA, nrow = 3, ncol = 3)
    for(i in 1:3){
      out_h1[i,] <- px_h1[i]*c(py1[i]*(1 - py0[i]),
                               py1[i]*py0[i] + (1 - py1[i])*(1 - py0[i]),
                               (1 - py1[i])* py0[i])
      
      out_h3[i,] <- px_h3[i]*c(py1[i]*(1 - py0[i]),
                               py1[i]*py0[i] + (1 - py1[i])*(1 - py0[i]),
                               (1 - py1[i])* py0[i])
    }
  }else{
    out_h1 <- out_h3 <- matrix(NA, nrow = 3*3, ncol = 3)
    for(i in 1:3){
      for(j in 1:3){
        out_h1[(i-1)*3+j,] <- px_h1[i]*px_h1[j]*c(py1[i]*(1 - py0[j]),
                                                  py1[i]*py0[j] + (1 - py1[i])*(1 - py0[j]),
                                                  (1 - py1[i])* py0[j])
        
        out_h3[(i-1)*3+j,] <- px_h3[i]*px_h3[j]*c(py1[i]*(1 - py0[j]),
                                                  py1[i]*py0[j] + (1 - py1[i])*(1 - py0[j]),
                                                  (1 - py1[i])* py0[j])
      }
    }
  }
  
  
  out_h1 <- colSums(out_h1)
  out_h3 <- colSums(out_h3)
  out <- c(out_h1, out_h3)
  names(out) <- c("p1", "p2", "p3", "q1", "q2", "q3")
  return(out)
}

# example
a <- 0.1000000
b <- 0.1000000 
ceoff <- c(-1.0461, -0.1270, 0.6424, 2.7752)
px_h1 <- c(a/(a + b), b/(a + b), 0)
px_h3 <- c(0, 0, 1)
pb_h(coeff, px_h1, px_h3, match_by = "h")






