%This function is used to sample parameters of the education and happiness using the algorithm given in the online supplement section B5.

function [delta_educ,scale_delta_educ,accept_delta_educ,educ_tilde]=sample_categorical_parameter(data,delta_educ,theta_gauss,prior,iter,scale_delta_educ,V1_delta_educ,...
                                                                                        accept_delta_educ,target_accept,income_tilde,health_tilde,happiness_tilde,dim_cop,dim_educ)


% construct the correlation matrix of the Gaussian copula.                                                                                    
chol_gauss=diag([ones(dim_cop,1)]);
number_gauss = 1;
for i=1:dim_cop-1
    for j=i:dim_cop-1
        chol_gauss(i,j+1) = theta_gauss(number_gauss,1); 
        number_gauss=number_gauss+1;
    end

end
sigma = chol_gauss'*chol_gauss;
diag_mat = diag([1./sqrt(diag(sigma))]);
covmat = diag_mat*sigma*diag_mat;                                                                                      


num_obs = length(data(:,1));
y_cts = [income_tilde,health_tilde,happiness_tilde];
cf = covmat(dim_educ, [1:(dim_educ-1) (dim_educ+1):end]) / covmat([1:(dim_educ-1) (dim_educ+1):end],...
        [1:(dim_educ-1) (dim_educ+1):end]);

cond_mean = cf*(y_cts)';    
cond_mean = cond_mean';
cond_var = covmat(dim_educ,dim_educ) - cf * covmat([1:(dim_educ-1) (dim_educ+1):end], dim_educ);
cond_sd = sqrt(cond_var);

num_category = length(unique(data(:,dim_educ)));
for k=1:num_category
    id(:,k) = (data(:,dim_educ)==k);    
end

for k=1:num_category-1
    if k==1
        limit_p(1,k) = ((-4 + exp(delta_educ(1,k))));
    else
        limit_p(1,k) = limit_p(1,k-1) + exp(delta_educ(1,k)); 
    end
    stand_limit_p(:,k) = (limit_p(1,k) - cond_mean)./cond_sd;
end

prior_delta = log(mvnpdf(delta_educ,zeros(1,num_category-1),prior.delta_var_educ));
stand_limit_p(:,num_category)=norminv(0.9999)*ones(num_obs,1);
comp_sum=0;
for k=1:num_category
    if k==1
       comp_sum = comp_sum + sum(id(:,k).*log(normcdf(stand_limit_p(:,k))));
    else
       comp_sum = comp_sum + sum(id(:,k).*log(normcdf(stand_limit_p(:,k)) - normcdf(stand_limit_p(:,k-1))));
    end   
end

post = prior_delta+comp_sum;
%random walk proposals
R1 = mvnrnd(delta_educ,scale_delta_educ.*V1_delta_educ);
delta_educ_star = R1;
for k=1:num_category-1
    if k==1
        limit_p_star(1,k) = ((-4 + exp(delta_educ_star(1,k))));
    else
        limit_p_star(1,k) = limit_p_star(1,k-1) + exp(delta_educ_star(1,k)); 
    end
    stand_limit_p_star(:,k) = (limit_p_star(1,k) - cond_mean)./cond_sd;
end

prior_delta_star = log(mvnpdf(delta_educ_star,zeros(1,num_category-1),prior.delta_var_educ));
stand_limit_p_star(:,num_category)=norminv(0.9999)*ones(num_obs,1);
comp_sum_star=0;
for k=1:num_category
    if k==1
       comp_sum_star = comp_sum_star + sum(id(:,k).*log(normcdf(stand_limit_p_star(:,k))));
    else
       comp_sum_star = comp_sum_star + sum(id(:,k).*log(normcdf(stand_limit_p_star(:,k)) - normcdf(stand_limit_p_star(:,k-1))));
    end   
end
post_star = prior_delta_star+comp_sum_star;
r1 = exp(post_star - post);
C1 = min(1,r1);
A1 = rand;
%accept/reject the new proposals
if A1<=C1
   accept_delta_educ=accept_delta_educ+1;
   delta_educ=delta_educ_star;
end

if iter>1000     
   scale_delta_educ = update_sigma(scale_delta_educ,C1,target_accept,iter,num_category-1);
end


for k=1:num_category-1
    if k==1
       limit_educ(1,k)= -4 + exp(delta_educ(1,k));
    else
       limit_educ(1,k)= limit_educ(1,k-1) + exp(delta_educ(1,k)); 
    end    
end
limit_educ(1,num_category) = norminv(0.9999);

%sampling latent variable
educ_tilde = zeros(num_obs,1);
for k=1:num_category
    if k==1
       educ_tilde(id(:,k),1) = trandn(-inf*ones(sum(id(:,k)),1),(limit_educ(1,k)-cond_mean(id(:,k),1))./cond_sd)*cond_sd+cond_mean(id(:,k),1);
    
    else
       educ_tilde(id(:,k),1) = trandn((limit_educ(1,k-1)-cond_mean(id(:,k),1))./cond_sd,(limit_educ(1,k)-cond_mean(id(:,k),1))./cond_sd)*cond_sd+cond_mean(id(:,k),1); 
    end
end




    
end
