#' Solver for quadratic equations a*x^2 + b*x + c = 0
#' 
#' @param a real number
#' @param b real number
#' @param c real number
#' 
#' @return solution of the quadratic equation a*x^2 + b*x + c = 0
#' @author Ted Harding, \email{Ted.Harding@@nessie.mcc.ac.uk}
qs <- function(a, b, c){
  a<-as.complex(a); b<-as.complex(b); c<-as.complex(c)
  i2<-(a!=0); i1<-((a==0)&(b!=0));
  solns<-as.complex(rep(NA,length(a)))
  solns<-cbind(solns,solns); colnames(solns)<-c("soln 1","soln 2")
  a2<-a[i2]; b2<-b[i2]; c2<-c[i2]
  solns[i2,1]<-(-b2 + sqrt(b2^2 - 4*a2*c2))/(2*a2)
  solns[i2,2]<-(-b2 - sqrt(b2^2 - 4*a2*c2))/(2*a2)
  b1<-b[i1]; c1<-c[i1]
  solns[i1,1]<-(-c1)/b1
  solns
}