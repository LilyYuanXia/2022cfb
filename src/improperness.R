library(hrbrthemes)
library(tidyverse)

# 1. Screening probability triples for that give cfb* < 0.5 binary X

# generate 101 possible values from [0, 1] for each probability 
# consider 101^4 possible combinations
p <- seq(0, 1, by = 0.01)
ps <- t(combn(p, m = 4))
nameMatrix <- rbind(c("p1", "p3", "q1", "q3"),
                    c("p1", "p3", "q3", "q1"),
                    c("p1", "q1", "p3", "q3") ,
                    c("p1", "q3", "p3", "q1"),
                    c("p1", "q3", "q1", "p3"),
                    c("p1", "q1", "q3", "p3") ,
                    c("p3", "p1", "q1", "q3"),
                    c("p3", "p1", "q3", "q1"),
                    c("p3", "q1", "p1", "q3") ,
                    c("p3", "q3", "p1", "q1") ,
                    c("p3", "q3", "q1", "p1") ,
                    c("p3", "q1", "q3", "p1") ,
                    c("q1", "p3", "p1", "q3") ,
                    c("q1", "p3", "q3", "p1") ,
                    c("q1", "p1", "p3", "q3") ,
                    c("q1", "q3", "p3", "p1") ,
                    c("q1", "q3", "p1", "p3") ,
                    c("q1", "p1", "q3", "p3") ,
                    c("q1", "p3", "p1", "q3"),
                    c("q3", "p3", "q1", "p1") ,
                    c("q3", "p3", "p1", "q1") ,
                    c("q3", "q1", "p3", "p1") ,
                    c("q3", "p1", "p3", "q1") ,
                    c("q3", "p1", "q1", "p3") ,
                    c("q3", "q1", "p1", "p3") )


outp1 <- outp3 <- outq1 <- outq3 <- NA

for(i in 1:nrow(nameMatrix)){
  colnames(ps) <- nameMatrix[i,]
  ps <- as.data.frame(ps)
  
  ps1 <- ps %>% filter(( p1 + p3 <= 1) &  (q1 + q3 <= 1))
  if(dim(ps1)[1] > 0){
    ps2 <- ps1 %>% filter(q1 - p1 - q1*p3 > q3 - p3 - q3*p1) # (q−1 - p−1 - q−1p+1 > q+1 - p+1 - q+1p−1)
    if(dim(ps2)[1] > 0){
      ps3 <- ps2 %>% filter(q1 - p1 < q3 - p3) #  (q+1 − q−1 > p+1 − p−1)
      if(dim(ps3)[1] > 0){
        message("found one")
        for(j in 1:nrow(ps3)){
          currentrow <- ps3[j,]
          outp1 <- c(outp1,currentrow$p1)
          outp3 <- c(outp3,currentrow$p3)
          outq1 <- c(outq1,currentrow$q1)
          outq3 <- c(outq3,currentrow$q3)
        }
      }else{
        message("no rows in level 3")
      }
    }else{
      message("no rows in level 2")
    }
  }else{
    message("no rows in level 1")
  }
}


out <- cbind(outp1, outp3, outq1, outq3)
out <- as.data.frame(na.omit(out))



# 2. Calculate cfb^* for probability triples

#out <- out[round(out$outq3 - out$ outp3,2) != round(out$outq1 - out$outp1,2), ]

#P(H1 > H2, B1 > B2)
nom1 <- 1/4 * (out$outq3*(1 - out$outp3 - out$outp1) + out$outq3*out$outp1 + (1 - out$outq3 - out$outq1)*out$outp1)
#P(H1 = H2, B1 > B2) 
nom2 <- 1/4 * (out$outq3*(1 - out$outq3 - out$outq1) +
                 out$outq3*out$outq1 +
                 (1 - out$outq3 - out$outq1)*out$outq1 +
                 out$outp3*(1 - out$outp3 - out$outp1) +
                 out$outp3*out$outp1 +
                 (1 - out$outp3 - out$outp1)*out$outp1)


#P(B1 > B2)
denom <- 1/4* ((out$outp3 + out$outq3)*(out$outp1 + out$outq1) + 
                 ((1 - out$outp3 - out$outp1) + (1 - out$outq3 - out$outq1))*(out$outp1 + out$outq1) +  
                 (out$outp3 + out$outq3)*((1 - out$outp3 - out$outp1) + (1 - out$outq3 - out$outq1)))

cfb <-  (nom1 + 0.5*nom2)/denom



# 3. Histogram of cfb (Figure 1)
summary(cfb)
cfb.total <- as.data.frame(cfb)

p <- cfb.total %>%  ggplot(aes(x = cfb)) +
  geom_histogram( binwidth=0.0003, fill="#69b3a2", color="#69b3a2", alpha=3) +
  theme(
    plot.title = element_text(size=15)
  ) + 
  geom_vline(xintercept = 0.4886, linetype="dashed", 
             color = "red", size=1) +
  geom_vline(xintercept = 0.4916, linetype="dashed", 
             color = "blue", size=1) +
  annotate("text", x=0.4188, y=10, label= "0.4188", col = "red", size = 3)+
  annotate("text", x=0.5, y=100, label= "0.5", col = "red", size = 3)+
  annotate("text", x=0.48, y=4000, label= "median = 0.4916", col = "blue", size = 3) +
  annotate("text", x=0.48, y=3500, label= "mean = 0.4886", col = "red", size = 3)

p + labs(x = "cfb",
         y = "Counts") +
  theme_ipsum() 

# 4. Improper scenarios with continuous H* (involve randomness)
triple1 <- c(0.08,0,0.92,0,0.15,0.85)
triple2 <- c(0.54,0.37,0.09,0.68,0.01,0.31)
sim_cfb(triple1)
sim_cfb(triple2)
