## motivating (hypotetical) example for Rb

s2_typ <- function(vi){
  I <- length(vi)
  wi <- 1/vi
  c("I2" = (I-1)*sum(wi)/(sum(wi)^2 - sum(wi^2)), "RI" = I/sum(wi))   
}

vi_a <- c(5.0, 5.2, 4.9, 5.3, 4.8)
#s2_typ(vi_a)
vi_b <- c(4, 17, 15, 2)
I <- length(vi_b)
last <- (I*sum(1/vi_b) - s2_typ(vi_a)[1]*sum(1/vi_b)^2 + s2_typ(vi_a)[1]*sum(1/vi_b^2))/
  (2*s2_typ(vi_a)[1]*sum(1/vi_b) - I)

vi_b <- unname(c(vi_b, round(1/last, 1)))
#vi_b

s2_typ(vi_b)
s2_typ(vi_a)

hyp_vi <- data.frame(
  analysis = c("A", "B"),
  v = c(paste(vi_a, collapse = ", "), paste(vi_b, collapse = ", ")),
  CV = c(sd(vi_a)/mean(vi_a), sd(vi_b)/mean(vi_b)),
  s2I2 = c(s2_typ(vi_a)[1], s2_typ(vi_b)[1]),
  s2RI = c(s2_typ(vi_a)[2], s2_typ(vi_b)[2])
)

save(hyp_vi, file = "data/hyp_vi.RData")
