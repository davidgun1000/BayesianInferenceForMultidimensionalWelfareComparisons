% compute the dirichelt pdf 

function [log_density] = logdir3pdf_general(x,phi)
%p1 = x(1,1);
%p2 = x(2,1);
%p3 = x(3,1);

size_x = size(x,1);
second_term = gammaln(phi(1));
third_term = (phi(1)-1)*log(x(1));
for i=2:size_x
    second_term = second_term +gammaln(phi(i));
    third_term = third_term+(phi(i)-1)*log(x(i));
end
log_density = gammaln(sum(phi)) - second_term + third_term;

%log_density = gammaln(phi1+phi2+phi3)-(gammaln(phi1)+gammaln(phi2)+gammaln(phi3))+(phi1-1)*log(p1)+(phi2-1)*log(p2)+(phi3-1)*log(p3);
density = exp(log_density);

end