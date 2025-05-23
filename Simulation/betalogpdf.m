%compute the log pdf of beta density

function logpdf = betalogpdf(x,s,m)
logpdf = ((s*m)-1)*log(x)+((s*(1-m))-1)*log(1-x)-betaln(s*m,s-s*m);