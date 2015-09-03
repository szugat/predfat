depth3Par <- function(res, type = c("full", "1", "2")) {
  type <- match.arg(type)
  signs <- sign(res)
  K <- 3 # Number of parameters
  if (type == "full") {
    combs <- combn(signs, K+1)
    return(mean(apply(combs, 2, function(x) all(x == (-1)**(0:K)) || all(x == (-1)**(1:(K+1))))))
  }
  
  if (type == "1") {
    starts <- seq(1, length(signs) - K, by = K + 1)
    ends <- seq(K+1, length(signs), by = K + 1)
    return(mean(sapply(seq_along(starts), function(x) all(signs[starts[x]:ends[x]] == (-1)**(0:K)) || all(signs[starts[x]:ends[x]] == (-1)**(1:(K+1))))))
  }
  
  if (type == "2") {
    starts <- 1:(length(signs)-K)
    ends <- starts + K
    return(mean(sapply(seq_along(starts), function(x) all(signs[starts[x]:ends[x]] == (-1)**(0:K)) || all(signs[starts[x]:ends[x]] == (-1)**(1:(K+1))))))
  }
}