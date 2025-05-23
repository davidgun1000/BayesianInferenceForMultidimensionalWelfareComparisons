%gamma3 mixture trivariate blocking

function [w_gamma,m_gamma,v_gamma,scale_gamma_com,accept_w_gamma,accept_gamma_com,income_tilde] = sampling_gamma_parameter(data,...
    w_gamma,m_gamma,v_gamma,theta_gauss,prior,iter,scale_gamma_com,V1_gamma_com,...
    accept_w_gamma,accept_gamma_com,num_component,target_accept,health_tilde,educ_tilde,happiness_tilde,dim_cop)


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
  


A1=rand;
n=length(data(:,1));
cdf1 = zeros(n,1);
dens1 = zeros(n,1);
for k=1:num_component
    dens1=dens1+w_gamma(k)*gampdf(data(:,1),v_gamma(k),m_gamma(k)/v_gamma(k));
    cdf1=cdf1+w_gamma(k)*gamcdf(data(:,1),v_gamma(k),m_gamma(k)/v_gamma(k));
end
lik_weight=sum(log(dens1));
income_tilde=norminv(cdf1);
prior_weight = logdir3pdf_general(w_gamma',[prior.phi*ones(1,num_component)]);
post_weight = prior_weight + lik_weight +...
        sum(log(mvnpdf([income_tilde,health_tilde,educ_tilde,happiness_tilde],zeros(1,dim_cop),covmat))) - ... 
               sum(log(normpdf(income_tilde,0,1))+log(normpdf(health_tilde,0,1))+log(normpdf(educ_tilde,0,1))+log(normpdf(happiness_tilde,0,1)));
w_gamma_star = dirich_rnd(prior.c1.*w_gamma')';
if sum((w_gamma_star>0.999)>0)
else
cdf1_star = 0;
dens1_star = 0;
for k=1:num_component
    dens1_star=dens1_star+w_gamma_star(k)*gampdf(data(:,1),v_gamma(k),m_gamma(k)/v_gamma(k));
    cdf1_star=cdf1_star+w_gamma_star(k)*gamcdf(data(:,1),v_gamma(k),m_gamma(k)/v_gamma(k));
end
lik_weight_star = sum(log(dens1_star));
income_tilde_star = norminv(cdf1_star);
prior_weight_star = logdir3pdf_general(w_gamma_star',[prior.phi*ones(1,num_component)]);
post_weight_star = prior_weight_star + lik_weight_star + sum(log(mvnpdf([income_tilde_star,health_tilde,educ_tilde,happiness_tilde],zeros(1,dim_cop),covmat))) - ... 
               sum(log(normpdf(income_tilde_star,0,1))+log(normpdf(health_tilde,0,1))+log(normpdf(educ_tilde,0,1))+log(normpdf(happiness_tilde,0,1)));
logwstar_w = logdir3pdf_general(w_gamma',[prior.c1.*w_gamma_star']);
%logwstar_w_alt = log(dir3pdf(w_gamma',prior.c1*w_gamma_star(1),prior.c1*w_gamma_star(2)));
logw_wstar = logdir3pdf_general(w_gamma_star',[prior.c1.*w_gamma']);
r1 = exp(post_weight_star - post_weight + logwstar_w - logw_wstar);
C1 = min(1,r1);
if A1<=C1
   w_gamma=w_gamma_star; 
   accept_w_gamma=accept_w_gamma+1; 
end
end
for k=1:num_component
    A1=rand();
    dens1=0;
    cdf1=0;
    for j=1:num_component
        dens1=dens1+w_gamma(j)*gampdf(data(:,1),v_gamma(j),m_gamma(j)/v_gamma(j));
        cdf1=cdf1+w_gamma(j)*gamcdf(data(:,1),v_gamma(j),m_gamma(j)/v_gamma(j));
    end
    lik_com=sum(log(dens1));
    income_tilde=norminv(cdf1);
    %prior_com = log_IG_PDF_used(m_gamma(k),prior.alpha_gamma,prior.beta_gamma)+log(exppdf(v_gamma(k),1/prior.theta_gamma));
    prior_com = log(gampdf(m_gamma(k),prior.alpha_gamma,prior.beta_gamma(k)))+log(exppdf(v_gamma(k),1/prior.theta_gamma));
    post_com = prior_com+lik_com+...     
               sum(log(mvnpdf([income_tilde,health_tilde,educ_tilde,happiness_tilde],zeros(1,dim_cop),covmat))) - ... 
               sum(log(normpdf(income_tilde,0,1))+log(normpdf(health_tilde,0,1))+log(normpdf(educ_tilde,0,1))+log(normpdf(happiness_tilde,0,1))); 
    jac_com = log(1/m_gamma(k)) + log(1/v_gamma(k));
    theta_com = [log(m_gamma(k)),log(v_gamma(k))];
    R1_com = mvnrnd(theta_com,scale_gamma_com(k).*V1_gamma_com(:,:,k));
    m_gamma_star = m_gamma;
    m_gamma_star(k) = exp(R1_com(1,1));
    v_gamma_star = v_gamma;
    v_gamma_star(k) = exp(R1_com(1,2));
    dens1_star=0;
    cdf1_star=0;
    for j=1:num_component
        dens1_star=dens1_star+w_gamma(j)*gampdf(data(:,1),v_gamma_star(j),m_gamma_star(j)/v_gamma_star(j));
        cdf1_star=cdf1_star+w_gamma(j)*gamcdf(data(:,1),v_gamma_star(j),m_gamma_star(j)/v_gamma_star(j));
    end
    lik_com_star = sum(log(dens1_star));
    income_tilde_star = norminv(cdf1_star);
    %prior_com_star = log_IG_PDF_used(m_gamma_star(k),prior.alpha_gamma,prior.beta_gamma)+log(exppdf(v_gamma_star(k),1/prior.theta_gamma));    
    prior_com_star = log(gampdf(m_gamma_star(k),prior.alpha_gamma,prior.beta_gamma(k)))+log(exppdf(v_gamma_star(k),1/prior.theta_gamma));
    post_com_star = prior_com_star+lik_com_star+...     
               sum(log(mvnpdf([income_tilde_star,health_tilde,educ_tilde,happiness_tilde],zeros(1,dim_cop),covmat))) - ... 
               sum(log(normpdf(income_tilde_star,0,1))+log(normpdf(health_tilde,0,1))+log(normpdf(educ_tilde,0,1))+log(normpdf(happiness_tilde,0,1)));  
    jac_com_star = log(1/m_gamma_star(k)) + log(1/v_gamma_star(k));
    r1 = exp(post_com_star - post_com +jac_com -jac_com_star);
    C1_com(k) = min(1,r1);
    if A1 <=C1_com(k)
       m_gamma(k) = m_gamma_star(k);
       v_gamma(k) = v_gamma_star(k);
       accept_gamma_com(k) = accept_gamma_com(k)+1;
             
    end

end

if iter>1000
     for k=1:num_component
         scale_gamma_com(k) = update_sigma(scale_gamma_com(k),C1_com(k),target_accept,iter,2);
%          if scale_gamma_com(k)>1
%             scale_gamma_com(k)=1; 
%          end
     end
end

%col = [m_gamma',v_gamma',w_gamma'];
%col_sort=sortrows(col,2);
%m_gamma = col_sort(:,1)';
%v_gamma = col_sort(:,2)';
%w_gamma = col_sort(:,3)';

dens1=0;
cdf1=0;
for j=1:num_component
    dens1=dens1+w_gamma(j)*gampdf(data(:,1),v_gamma(j),m_gamma(j)/v_gamma(j));
    cdf1=cdf1+w_gamma(j)*gamcdf(data(:,1),v_gamma(j),m_gamma(j)/v_gamma(j));
end
income_tilde=norminv(cdf1);
end
% theta_gaussian = [1 theta12 theta13; 0 1 theta23; 0 0 1];
% sig = (theta_gaussian'*theta_gaussian);
% diag_mat = [(1/sqrt(sig(1,1))) 0 0 ; 0 (1/sqrt(sig(2,2))) 0 ; 0 0 (1/sqrt(sig(3,3)))];
% cov_mat = diag_mat*sig*diag_mat;
% 
% 
%      copula_data_income = w1_income*gamcdf(x,v1_income,m1_income/v1_income) + ...
%     w2_income*gamcdf(x,v2_income,m2_income/v2_income) + ...
%     w3_income*gamcdf(x,v3_income,m3_income/v3_income);
%         income_star = norminv(copula_data_income,0,1);
%      %  sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%         %       sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)))
%         
%     A_w1_income = rand();
%      c1 = 10000;
%      c2 = 10000;
%      weight_income_star = dirich_rnd([c1*w1_income;c1*w2_income;c1*(1-w1_income-w2_income)]);
%       w1_income_star = weight_income_star(1,1);
%       w2_income_star = weight_income_star(2,1);
%       w3_income_star = weight_income_star(3,1);
%      prior_w1_income = dir3pdf([w1_income; w2_income; w3_income],phi1,phi2,phi3);
%      lik_w1_income = sum(log(w1_income*gampdf(x,v1_income,m1_income/v1_income) + ...
%         w2_income*gampdf(x,v2_income,m2_income/v2_income) + ...
%         w3_income*gampdf(x,v3_income,m3_income/v3_income)));
%      copula_data_income = w1_income*gamcdf(x,v1_income,m1_income/v1_income) + ...
%      w2_income*gamcdf(x,v2_income,m2_income/v2_income) + ...
%      w3_income*gamcdf(x,v3_income,m3_income/v3_income);
%      income_star = norminv(copula_data_income,0,1);
%      post_w1_income = prior_w1_income + lik_w1_income + ...
%          sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      prior_w1_income_star = dir3pdf([w1_income_star; w2_income_star; w3_income_star],phi1,phi2,phi3);
%      lik_w1_income_star = sum(log(w1_income_star*gampdf(x,v1_income,m1_income/v1_income) + ...
%         w2_income_star*gampdf(x,v2_income,m2_income/v2_income) + ...
%         w3_income_star*gampdf(x,v3_income,m3_income/v3_income)));
%      copula_data_income_w1_star = w1_income_star*gamcdf(x,v1_income,m1_income/v1_income) + ...
%      w2_income_star*gamcdf(x,v2_income,m2_income/v2_income) + ...
%      w3_income_star*gamcdf(x,v3_income,m3_income/v3_income);
%      income_star_star_w1 = norminv(copula_data_income_w1_star,0,1);
%      post_w1_income_star = prior_w1_income_star + lik_w1_income_star + ...
%                sum(log(mvnpdf([income_star_star_w1 health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star_star_w1,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
% %      logw1star_w1 = log(betaden(w1_income,c1,w1_income_star));
% %      logw1_w1star = log(betaden(w1_income_star,c1,w1_income));
% %      logw2star_w2 = log(betaden(w2_income,c2,w2_income_star));
% %      logw2_w2star = log(betaden(w2_income_star,c2,w2_income));
%      logwstar_w = log(dir3pdf([w1_income;w2_income;w3_income],c1*w1_income_star,c1*w2_income_star,c1*w3_income_star));
%      logw_wstar = log(dir3pdf([w1_income_star;w2_income_star;w3_income_star],c1*w1_income,c1*w2_income,c1*w3_income));
%      rat_w1_income = exp(post_w1_income_star - post_w1_income + logwstar_w - logw_wstar);
%      MH_w1_income = min(1,rat_w1_income);
%      if A_w1_income <= MH_w1_income
%          w1_income = w1_income_star;
%          w2_income = w2_income_star;
%          w3_income = 1-w1_income-w2_income;
%          income_star = income_star_star_w1;
%          accept0 = 1;
%      else
%          w1_income = w1_income;
%          w2_income = w2_income;
%          w3_income = 1-w1_income-w2_income;
%          income_star = income_star;
%          accept0 = 0;
%      end
%      
% 
%      
%      A_m1_income = rand();
%      r_m1_income = 500;
%      m1_income_star =gamrnd(r_m1_income,m1_income/r_m1_income);
%      
%      prior_m1_income = log(IG_PDF_used(m1_income,alpha1,beta1));
%      lik_m1_income = sum(log(w1_income*gampdf(x,v1_income,m1_income/v1_income) + ...
%         w2_income*gampdf(x,v2_income,m2_income/v2_income) + ...
%         w3_income*gampdf(x,v3_income,m3_income/v3_income)));
%      post_m1_income = prior_m1_income + lik_m1_income + ...
%          sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      prior_m1_income_star = log(IG_PDF_used(m1_income_star,alpha1,beta1));
%      lik_m1_income_star = sum(log(w1_income*gampdf(x,v1_income,m1_income_star/v1_income) + ...
%         w2_income*gampdf(x,v2_income,m2_income/v2_income) + ...
%         w3_income*gampdf(x,v3_income,m3_income/v3_income)));
%     
%      copula_data_income_m1_star = w1_income*gamcdf(x,v1_income,m1_income_star/v1_income) + ...
%       w2_income*gamcdf(x,v2_income,m2_income/v2_income) + ...
%       w3_income*gamcdf(x,v3_income,m3_income/v3_income);
%      income_star_star_m1 = norminv(copula_data_income_m1_star,0,1);
%      post_m1_income_star = prior_m1_income_star + lik_m1_income_star + ...
%                sum(log(mvnpdf([income_star_star_m1 health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star_star_m1,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      logm1star_m1 = log(gampdf(m1_income,r_m1_income,m1_income_star/r_m1_income));
%      logm1_m1star = log(gampdf(m1_income_star,r_m1_income,m1_income/r_m1_income));
%      rat_m1_income = exp(post_m1_income_star - post_m1_income +logm1star_m1 - logm1_m1star);
%      MH_m1_income = min(1,rat_m1_income);
%      if A_m1_income<=MH_m1_income
%          m1_income = m1_income_star;
%          income_star = income_star_star_m1;
%          accept1=1;
%      else
%          m1_income = m1_income;
%          income_star = income_star;
%          accept1=0;
%      end
%      
%      A_m2_income = rand();
%      r_m2_income = 1000;
%      m2_income_star =gamrnd(r_m2_income,m2_income/r_m2_income);
%      
%      prior_m2_income = log(IG_PDF_used(m2_income,alpha2,beta2));
%      lik_m2_income = sum(log(w1_income*gampdf(x,v1_income,m1_income/v1_income) + ...
%         w2_income*gampdf(x,v2_income,m2_income/v2_income) + ...
%         w3_income*gampdf(x,v3_income,m3_income/v3_income)));
%      post_m2_income = prior_m2_income + lik_m2_income + ...
%                sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      prior_m2_income_star = log(IG_PDF_used(m2_income_star,alpha2,beta2));
%      lik_m2_income_star = sum(log(w1_income*gampdf(x,v1_income,m1_income/v1_income) + ...
%         w2_income*gampdf(x,v2_income,m2_income_star/v2_income) + ...
%         w3_income*gampdf(x,v3_income,m3_income/v3_income)));
%      copula_data_income_m2_star = w1_income*gamcdf(x,v1_income,m1_income/v1_income) + ...
%       w2_income*gamcdf(x,v2_income,m2_income_star/v2_income) + ...
%       w3_income*gamcdf(x,v3_income,m3_income/v3_income);
%      income_star_star_m2 = norminv(copula_data_income_m2_star,0,1);
%      post_m2_income_star = prior_m2_income_star + lik_m2_income_star + ...
%                sum(log(mvnpdf([income_star_star_m2 health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star_star_m2,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      logm2star_m2 = log(gampdf(m2_income,r_m2_income,m2_income_star/r_m2_income));
%      logm2_m2star = log(gampdf(m2_income_star,r_m2_income,m2_income/r_m2_income));
%      rat_m2_income = exp(post_m2_income_star - post_m2_income +logm2star_m2 - logm2_m2star);
%      MH_m2_income = min(1,rat_m2_income);
%      if A_m2_income<=MH_m2_income
%          m2_income = m2_income_star;
%          income_star = income_star_star_m2;
%          accept2 = 1;
%      else
%          m2_income = m2_income;
%          income_star = income_star;
%          accept2=0;
%      end
%      
%      A_m3_income = rand();
%      r_m3_income = 1000;
%      m3_income_star =gamrnd(r_m3_income,m3_income/r_m3_income);
%      
%      prior_m3_income = log(IG_PDF_used(m3_income,alpha3,beta3));
%      lik_m3_income = sum(log(w1_income*gampdf(x,v1_income,m1_income/v1_income) + ...
%         w2_income*gampdf(x,v2_income,m2_income/v2_income) + ...
%         w3_income*gampdf(x,v3_income,m3_income/v3_income)));
%      post_m3_income = prior_m3_income + lik_m3_income + ...
%                sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      prior_m3_income_star = log(IG_PDF_used(m3_income_star,alpha3,beta3));
%      lik_m3_income_star = sum(log(w1_income*gampdf(x,v1_income,m1_income/v1_income) + ...
%         w2_income*gampdf(x,v2_income,m2_income/v2_income) + ...
%         w3_income*gampdf(x,v3_income,m3_income_star/v3_income)));
%      copula_data_income_m3_star = w1_income*gamcdf(x,v1_income,m1_income/v1_income) + ...
%       w2_income*gamcdf(x,v2_income,m2_income/v2_income) + ...
%       w3_income*gamcdf(x,v3_income,m3_income_star/v3_income);
%      income_star_star_m3 = norminv(copula_data_income_m3_star,0,1);
%      post_m3_income_star = prior_m3_income_star + lik_m3_income_star + ...
%                sum(log(mvnpdf([income_star_star_m3 health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star_star_m3,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      logm3star_m3 = log(gampdf(m3_income,r_m3_income,m3_income_star/r_m3_income));
%      logm3_m3star = log(gampdf(m3_income_star,r_m3_income,m3_income/r_m3_income));
%      rat_m3_income = exp(post_m3_income_star - post_m3_income +logm3star_m3 - logm3_m3star);
%      MH_m3_income = min(1,rat_m3_income);
%      if A_m3_income<=MH_m3_income
%          m3_income = m3_income_star;
%          income_star = income_star_star_m3;
%      else
%          m3_income = m3_income;
%          income_star = income_star;
%      end
%      
%    
%      
%      
%      A_v1_income = rand();
%      r_v1_income = 500;
%      v1_income_star = gamrnd(r_v1_income,v1_income/r_v1_income);     
%      prior_v1_income = log(exppdf(v1_income,1/theta1));
%      lik_v1_income = sum(log(w1_income*gampdf(x,v1_income,m1_income/v1_income) + ...
%         w2_income*gampdf(x,v2_income,m2_income/v2_income) + ...
%         w3_income*gampdf(x,v3_income,m3_income/v3_income)));
%      post_v1_income = prior_v1_income + lik_v1_income + ...
%                sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      prior_v1_income_star = log(exppdf(v1_income_star,1/theta1));
%      lik_v1_income_star = sum(log(w1_income*gampdf(x,v1_income_star,m1_income/v1_income_star) + ...
%         w2_income*gampdf(x,v2_income,m2_income/v2_income) + ...
%         w3_income*gampdf(x,v3_income,m3_income/v3_income)));
%      copula_data_income_v1_star = w1_income*gamcdf(x,v1_income_star,m1_income/v1_income_star) + ...
%       w2_income*gamcdf(x,v2_income,m2_income/v2_income) + ...
%       w3_income*gamcdf(x,v3_income,m3_income/v3_income);
%      income_star_star_v1 = norminv(copula_data_income_v1_star,0,1);
%      post_v1_income_star = prior_v1_income_star + lik_v1_income_star + ...
%                sum(log(mvnpdf([income_star_star_v1 health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star_star_v1,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      logv1star_v1 = log(gampdf(v1_income,r_v1_income,v1_income_star/r_v1_income));
%      logv1_v1star = log(gampdf(v1_income_star,r_v1_income,v1_income/r_v1_income));
%      rat_v1_income = exp(post_v1_income_star - post_v1_income + logv1star_v1 - logv1_v1star);
%      MH_v1_income = min(1,rat_v1_income);
%      if A_v1_income <= MH_v1_income
%          v1_income = v1_income_star;
%          income_star = income_star_star_v1;
%      else
%          v1_income = v1_income;
%          income_star = income_star;
%      end
%      
%      A_v2_income = rand();
%      r_v2_income = 500;
%      v2_income_star = gamrnd(r_v2_income,v2_income/r_v2_income);
%      
%      prior_v2_income = log(exppdf(v2_income,1/theta2));
%      lik_v2_income = sum(log(w1_income*gampdf(x,v1_income,m1_income/v1_income) + ...
%         w2_income*gampdf(x,v2_income,m2_income/v2_income) + ...
%         w3_income*gampdf(x,v3_income,m3_income/v3_income)));
%      post_v2_income = prior_v2_income + lik_v2_income + ...
%                sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      prior_v2_income_star = log(exppdf(v2_income_star,1/theta2));
%      lik_v2_income_star = sum(log(w1_income*gampdf(x,v1_income,m1_income/v1_income) + ...
%         w2_income*gampdf(x,v2_income_star,m2_income/v2_income_star) + ...
%         w3_income*gampdf(x,v3_income,m3_income/v3_income)));
%     copula_data_income_v2_star = w1_income*gamcdf(x,v1_income,m1_income/v1_income) + ...
%       w2_income*gamcdf(x,v2_income_star,m2_income/v2_income_star) + ...
%       w3_income*gamcdf(x,v3_income,m3_income/v3_income);
%     income_star_star_v2 = norminv(copula_data_income_v2_star,0,1);
%      post_v2_income_star = prior_v2_income_star + lik_v2_income_star+ ...
%                sum(log(mvnpdf([income_star_star_v2 health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star_star_v2,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      logv2star_v2 = log(gampdf(v2_income,r_v2_income,v2_income_star/r_v2_income));
%      logv2_v2star = log(gampdf(v2_income_star,r_v2_income,v2_income/r_v2_income));
%      rat_v2_income = exp(post_v2_income_star - post_v2_income + logv2star_v2 - logv2_v2star);
%      MH_v2_income = min(1,rat_v2_income);
%      if A_v2_income <= MH_v2_income
%          v2_income = v2_income_star;
%          income_star = income_star_star_v2;
%      else
%          v2_income = v2_income;
%          income_star = income_star;
%      end
%      
%      A_v3_income = rand();
%      r_v3_income = 500;
%      v3_income_star = gamrnd(r_v3_income,v3_income/r_v3_income);
%      
%      prior_v3_income = log(exppdf(v3_income,1/theta3));
%      lik_v3_income = sum(log(w1_income*gampdf(x,v1_income,m1_income/v1_income) + ...
%         w2_income*gampdf(x,v2_income,m2_income/v2_income) + ...
%         w3_income*gampdf(x,v3_income,m3_income/v3_income)));
%      post_v3_income = prior_v3_income + lik_v3_income + ...
%                sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      prior_v3_income_star = log(exppdf(v3_income_star,1/theta3));
%      lik_v3_income_star = sum(log(w1_income*gampdf(x,v1_income,m1_income/v1_income) + ...
%         w2_income*gampdf(x,v2_income,m2_income/v2_income) + ...
%         w3_income*gampdf(x,v3_income_star,m3_income/v3_income_star)));
%      copula_data_income_v3_star = w1_income*gamcdf(x,v1_income,m1_income/v1_income) + ...
%       w2_income*gamcdf(x,v2_income,m2_income/v2_income) + ...
%       w3_income*gamcdf(x,v3_income_star,m3_income/v3_income_star);
%      income_star_star_v3 = norminv(copula_data_income_v3_star,0,1);
%      post_v3_income_star = prior_v3_income_star + lik_v3_income_star + ...
%           sum(log(mvnpdf([income_star_star_v3 health_star educ_star],[0 0 0],cov_mat ))) - ... 
%                sum(log(normpdf(income_star_star_v3,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%      logv3star_v3 = log(gampdf(v3_income,r_v3_income,v3_income_star/r_v3_income));
%      logv3_v3star = log(gampdf(v3_income_star,r_v3_income,v3_income/r_v3_income));
%      rat_v3_income = exp(post_v3_income_star - post_v3_income + logv3star_v3 - logv3_v3star);
%      MH_v3_income = min(1,rat_v3_income);
%      if A_v3_income <= MH_v3_income
%          v3_income = v3_income_star;
%          income_star = income_star_star_v3;
%      else
%          v3_income = v3_income;
%          income_star = income_star;
%      end    
%         
%              copula_data_income = w1_income*gamcdf(x,v1_income,m1_income/v1_income) + ...
%     w2_income*gamcdf(x,v2_income,m2_income/v2_income) + ...
%     w3_income*gamcdf(x,v3_income,m3_income/v3_income);
%         income_star = norminv(copula_data_income,0,1);


















%     
%     if A11<MH1_copula
%         m1_income = m1_income_new;
%         v1_income = v1_income_new;
%         income_star = income_star_star1;
%         accept1=1;
%     else
%         m1_income = m1_income;
%         v1_income = v1_income;
%         income_star = income_star;
%         accept1=0;
%     end
%     
% %-------------------------------------------------------------------------------------------------------------
% m2_income_new=1/(gamrnd(alpha2+N2*v2_income,1/(beta2+S2*v2_income)));
% 
% r2=600;
% 
% A2=rand(1);
% 
% V2_income=gamrnd(r2,v2_income/r2);
% 
% logvv2=v2_income*N2*log(v2_income)-N2*gammaln(v2_income)-v2_income*(theta2+(S2/m2_income_new)+N2*log(m2_income_new)-logP2);
% logvc2=V2_income*N2*log(V2_income)-N2*gammaln(V2_income)-V2_income*(theta2+(S2/m2_income_new)+N2*log(m2_income_new)-logP2);
% logV_v2=log(gampdf(v2_income,r2,V2_income/r2));
% logv_V2=log(gampdf(V2_income,r2,v2_income/r2));
% MH2=exp(logvc2+logV_v2-logvv2-logv_V2);
% C2=min(1,MH2);
%     if A2<=C2
%        v2_income_new=V2_income;
%     else
%        v2_income_new=v2_income;
%     end
%     
%     %--------------------------------------------------------------------------------------------------------------------------------
%     %accept/reject 2nd component
%     
%     copula_data_income_star2 = w1_income*gamcdf(x,v1_income,m1_income/v1_income) + ...
%     w2_income*gamcdf(x,v2_income_new,m2_income_new/v2_income_new) + ...
%     w3_income*gamcdf(x,v3_income,m3_income/v3_income);
% %     id2_1 = copula_data_income_star2 <= 0;
% %     copula_data_income_star2(id2_1,1) = 0.0001;
% %     id2_2 = copula_data_income_star2 >= 0.9999;
% %     copula_data_income_star2(id2_2,1) = 0.9999;    
%     income_star_star2 = norminv(copula_data_income_star2,0,1);
% 
%     A11 = rand();
%         %I = diag([1;1;1]);
%     log_lik2 = sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%              sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%     log_lik_star2 = sum(log(mvnpdf([income_star_star2 health_star educ_star],[0 0 0],cov_mat )))-...
%             sum(log(normpdf(income_star_star2,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%         
%         ratio_MH2 = exp(log_lik_star2 - log_lik2);
%     
%         MH2_copula = min(1,ratio_MH2);
%     
%     if A11<MH2_copula
%         m2_income = m2_income_new;
%         v2_income = v2_income_new;
%         income_star = income_star_star2;
%         accept2=1;
%     else
%         m2_income = m2_income;
%         v2_income = v2_income;
%         income_star = income_star;
%         accept2=0;
%     end
%     
% %     %------------------------------------------------------------------------------------------------------------------------------
% %     
%         m3_income_new=1/(gamrnd(alpha3+N3*v3_income,1/(beta3+S3*v3_income)));
%         r3=600;
%         A3=rand(1);
%         V3_income=gamrnd(r3,v3_income/r3); 
%         logvv3=v3_income*N3*log(v3_income)-N3*gammaln(v3_income)-v3_income*(theta3+(S3/m3_income_new)+N3*log(m3_income_new)-logP3);
%         logvc3=V3_income*N3*log(V3_income)-N3*gammaln(V3_income)-V3_income*(theta3+(S3/m3_income_new)+N3*log(m3_income_new)-logP3);
%         logV_v3=log(gampdf(v3_income,r3,V3_income/r3));
%         logv_V3=log(gampdf(V3_income,r3,v3_income/r3));
%         MH3=exp(logvc3+logV_v3-logvv3-logv_V3);
%         C3=min(1,MH3);
%             if A3<=C3
%                 v3_income_new=V3_income;
%             else
%                 v3_income_new=v3_income;
%             end
%             
%             %--------------------------------------------------------------------------------------------------------------------------------
%     %accept/reject 3rd component
%     
%     copula_data_income_star3 = w1_income*gamcdf(x,v1_income,m1_income/v1_income) + ...
%     w2_income*gamcdf(x,v2_income,m2_income/v2_income) + ...
%     w3_income*gamcdf(x,v3_income_new,m3_income_new/v3_income_new);
%     
% 
% %     id3_1 = copula_data_income_star3 <= 0;
% %     copula_data_income_star3(id3_1,1) = 0.0001;
% %     id3_2 = copula_data_income_star3 >=0.9999;
% %     copula_data_income_star3(id3_2,1) = 0.999999;
%     income_star_star3 = norminv(copula_data_income_star3,0,1);
%     
%     A11 = rand();
%         %I = diag([1;1;1]);
%     log_lik3 = sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%              sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%     log_lik_star3 = sum(log(mvnpdf([income_star_star3 health_star educ_star],[0 0 0],cov_mat )))-...
%             sum(log(normpdf(income_star_star3,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%         
%         ratio_MH3 = exp(log_lik_star3 - log_lik3);
%     
%         MH3_copula = min(1,ratio_MH3);
%     
%     if A11<MH3_copula
%         m3_income = m3_income_new;
%         v3_income = v3_income_new;
%         income_star = income_star_star3;
%         accept3=1;
%     else
%         m3_income = m3_income;
%         v3_income = v3_income;
%         income_star = income_star;
%         accept3=0;
%     end
    %------------------------------------------------------------------------------------------------------------------------------
    
%  




















































%------------------------------------------------------------------------------------------------------------------------------------
% %weight_new = dirich_rnd([phi1+N1 ; phi2+N2 ; phi3+N3]);
% w1_income_new = weight_new(1,1);
% w2_income_new = weight_new(2,1);
% w3_income_new = weight_new(3,1);
% 
% %--------------------------------------------------------------------------------------------------------------------------------
%     %accept/reject weight component
%     
%     copula_data_income_star1 = w1_income_new*gamcdf(x,v1_income,m1_income/v1_income) + ...
%     w2_income_new*gamcdf(x,v2_income,m2_income/v2_income) + ...
%     w3_income_new*gamcdf(x,v3_income,m3_income/v3_income);
%         income_star_star1 = norminv(copula_data_income_star1,0,1);
% %     id1_1 = copula_data_income_star1 <= 0;
% %     copula_data_income_star1(id1_1,1) = 0.0001;
% %     id1_2 = copula_data_income_star1 >= 0.9999;
% %     copula_data_income_star1(id1_2,1) = 0.9999;    
%         A11 = rand();
%         %I = diag([1;1;1]);
%          log_lik1 = sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%              sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%          log_lik_star1 = sum(log(mvnpdf([income_star_star1 health_star educ_star],[0 0 0],cov_mat )))-...
%             sum(log(normpdf(income_star_star1,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%         
%         ratio_MH1 = exp(log_lik_star1 - log_lik1);
%     
%         MH1_copula = min(1,ratio_MH1);
%     
%     if A11<MH1_copula
%         w1_income = w1_income_new;
%         w2_income = w2_income_new;
%         w3_income = w3_income_new;
%         income_star = income_star_star1;
%         accept1=1;
%     else
%         w1_income = w1_income;
%         w2_income = w2_income;
%         w3_income = w3_income;
%         income_star = income_star;
%         accept1=0;
%     end
% %-----------------------------------------------------------------------------------------
% %second latent variable generation
%  Z1=w1_income*gampdf(x,v1_income,m1_income/v1_income);
%  Z2=w2_income*gampdf(x,v2_income,m2_income/v2_income);
%  Z3=w3_income*gampdf(x,v3_income,m3_income/v3_income);
%  PZ1=Z1./(Z1+Z2+Z3);
%  PZ2=Z2./(Z1+Z2+Z3);
%  PZ3=Z3./(Z1+Z2+Z3);
%  Z=rand(n,1);
%  z1 = (Z <=PZ1);
%  z2 = (Z > PZ1) & (Z <= PZ1+PZ2);
%  logx=log(x);
% N1=sum(z1);
% N2=sum(z2);
% N3=n-N1-N2;
% S1=sum(x.*(z1));
% S2=sum(x.*(z2));
% S3=sum(x.*((1-z1-z2)));
% logP1=sum(logx.*(z1));
% logP2=sum(logx.*(z2));
% logP3=sum(logx.*((1-z1-z2)));
% %     
% 
%     %--------------------------------------------------------------------------------------------------------------
% 
% m1_income_new=1/(gamrnd(alpha1+N1*v1_income,1/(beta1+S1*v1_income)));
% m2_income_new=1/(gamrnd(alpha2+N2*v2_income,1/(beta2+S2*v2_income)));
% m3_income_new=1/(gamrnd(alpha3+N3*v3_income,1/(beta3+S3*v3_income)));
% 
%     %--------------------------------------------------------------------------------------------------------------------------------
%     %accept/reject 2nd component
%     
%     copula_data_income_star2 = w1_income*gamcdf(x,v1_income,m1_income_new/v1_income) + ...
%     w2_income*gamcdf(x,v2_income,m2_income_new/v2_income) + ...
%     w3_income*gamcdf(x,v3_income,m3_income_new/v3_income);
% %     id2_1 = copula_data_income_star2 <= 0;
% %     copula_data_income_star2(id2_1,1) = 0.0001;
% %     id2_2 = copula_data_income_star2 >= 0.9999;
% %     copula_data_income_star2(id2_2,1) = 0.9999;    
%     income_star_star2 = norminv(copula_data_income_star2,0,1);
% 
%     A11 = rand();
%         %I = diag([1;1;1]);
%     log_lik2 = sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%              sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%     log_lik_star2 = sum(log(mvnpdf([income_star_star2 health_star educ_star],[0 0 0],cov_mat )))-...
%             sum(log(normpdf(income_star_star2,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%         
%         ratio_MH2 = exp(log_lik_star2 - log_lik2);
%     
%         MH2_copula = min(1,ratio_MH2);
%     
%     if A11<MH2_copula
%         m1_income = m1_income_new;
%         m2_income = m2_income_new;
%         m3_income = m3_income_new;
%         income_star = income_star_star2;
%         accept2=1;
%     else
%         m1_income = m1_income;
%         m2_income = m2_income;
%         m3_income = m3_income;
%         income_star = income_star;
%         accept2=0;
%     end
%     
%     %------------------------------------------------------------------------------------------------------------------------------
% %third latent variable generation
%  Z1=w1_income*gampdf(x,v1_income,m1_income/v1_income);
%  Z2=w2_income*gampdf(x,v2_income,m2_income/v2_income);
%  Z3=w3_income*gampdf(x,v3_income,m3_income/v3_income);
%  PZ1=Z1./(Z1+Z2+Z3);
%  PZ2=Z2./(Z1+Z2+Z3);
%  PZ3=Z3./(Z1+Z2+Z3);
%  Z=rand(n,1);
%  z1 = (Z <=PZ1);
%  z2 = (Z > PZ1) & (Z <= PZ1+PZ2);
%  logx=log(x);
% N1=sum(z1);
% N2=sum(z2);
% N3=n-N1-N2;
% S1=sum(x.*(z1));
% S2=sum(x.*(z2));
% S3=sum(x.*((1-z1-z2)));
% logP1=sum(logx.*(z1));
% logP2=sum(logx.*(z2));
% logP3=sum(logx.*((1-z1-z2)));
%     
% 
%     %--------------------------------------------------------------------------------------------------------------
% r1=3000;
% 
% A1=rand(1);
% 
% V1_income=gamrnd(r1,v1_income/r1);
% 
% logvv1=v1_income*N1*log(v1_income)-N1*gammaln(v1_income)-v1_income*(theta1+(S1/m1_income)+N1*log(m1_income)-logP1);
% logvc1=V1_income*N1*log(V1_income)-N1*gammaln(V1_income)-V1_income*(theta1+(S1/m1_income)+N1*log(m1_income)-logP1);
% logV_v1=log(gampdf(v1_income,r1,V1_income/r1));
% logv_V1=log(gampdf(V1_income,r1,v1_income/r1));
% MH1=exp(logvc1+logV_v1-logvv1-logv_V1);
% C1=min(1,MH1);
%     if A1<=C1
%        v1_income_new=V1_income;
%     else
%        v1_income_new=v1_income;
%     end
% 
% r2=3000;
% 
% A2=rand(1);
% 
% V2_income=gamrnd(r2,v2_income/r2);
% 
% logvv2=v2_income*N2*log(v2_income)-N2*gammaln(v2_income)-v2_income*(theta2+(S2/m2_income)+N2*log(m2_income)-logP2);
% logvc2=V2_income*N2*log(V2_income)-N2*gammaln(V2_income)-V2_income*(theta2+(S2/m2_income)+N2*log(m2_income)-logP2);
% logV_v2=log(gampdf(v2_income,r2,V2_income/r2));
% logv_V2=log(gampdf(V2_income,r2,v2_income/r2));
% MH2=exp(logvc2+logV_v2-logvv2-logv_V2);
% C2=min(1,MH2);
%     if A2<=C2
%        v2_income_new=V2_income;
%     else
%        v2_income_new=v2_income;
%     end
% 
%         r3=3000;
%         A3=rand(1);
%         V3_income=gamrnd(r3,v3_income/r3); 
%         logvv3=v3_income*N3*log(v3_income)-N3*gammaln(v3_income)-v3_income*(theta3+(S3/m3_income)+N3*log(m3_income)-logP3);
%         logvc3=V3_income*N3*log(V3_income)-N3*gammaln(V3_income)-V3_income*(theta3+(S3/m3_income)+N3*log(m3_income)-logP3);
%         logV_v3=log(gampdf(v3_income,r3,V3_income/r3));
%         logv_V3=log(gampdf(V3_income,r3,v3_income/r3));
%         MH3=exp(logvc3+logV_v3-logvv3-logv_V3);
%         C3=min(1,MH3);
%             if A3<=C3
%                 v3_income_new=V3_income;
%             else
%                 v3_income_new=v3_income;
%             end
% 
%             %--------------------------------------------------------------------------------------------------------------------------------
%     %accept/reject v component
%     
%     copula_data_income_star3 = w1_income*gamcdf(x,v1_income_new,m1_income/v1_income_new) + ...
%     w2_income*gamcdf(x,v2_income_new,m2_income/v2_income_new) + ...
%     w3_income*gamcdf(x,v3_income_new,m3_income/v3_income_new);
%     
% 
% %     id3_1 = copula_data_income_star3 <= 0;
% %     copula_data_income_star3(id3_1,1) = 0.0001;
% %     id3_2 = copula_data_income_star3 >=0.9999;
% %     copula_data_income_star3(id3_2,1) = 0.999999;
%     income_star_star3 = norminv(copula_data_income_star3,0,1);
%     
%     A11 = rand();
%         %I = diag([1;1;1]);
%     log_lik3 = sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%              sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%     log_lik_star3 = sum(log(mvnpdf([income_star_star3 health_star educ_star],[0 0 0],cov_mat )))-...
%             sum(log(normpdf(income_star_star3,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%         
%         ratio_MH3 = exp(log_lik_star3 - log_lik3);
%     
%         MH3_copula = min(1,ratio_MH3);
%     
%     if A11<MH3_copula
%         v1_income = v1_income_new;
%         v2_income = v2_income_new;
%         v3_income = v3_income_new;
%         income_star = income_star_star3;
%         accept3=1;
%     else
%         v1_income = v1_income;
%         v2_income = v2_income;
%         v3_income = v3_income;
%         income_star = income_star;
%         accept3=0;
%     end
    
    %------------------------------------------------------------------------------------------------------------------------------
    















 






























% m1_income_new=1/(gamrnd(alpha1+N1*v1_income,1/(beta1+S1*v1_income)));
% 
% r1=600;
% 
% A1=rand(1);
% 
% V1_income=gamrnd(r1,v1_income/r1);
% 
% logvv1=v1_income*N1*log(v1_income)-N1*gammaln(v1_income)-v1_income*(theta1+(S1/m1_income_new)+N1*log(m1_income_new)-logP1);
% logvc1=V1_income*N1*log(V1_income)-N1*gammaln(V1_income)-V1_income*(theta1+(S1/m1_income_new)+N1*log(m1_income_new)-logP1);
% logV_v1=log(gampdf(v1_income,r1,V1_income/r1));
% logv_V1=log(gampdf(V1_income,r1,v1_income/r1));
% MH1=exp(logvc1+logV_v1-logvv1-logv_V1);
% C1=min(1,MH1);
%     if A1<=C1
%        v1_income_new=V1_income;
%     else
%        v1_income_new=v1_income;
%     end
% %--------------------------------------------------------------------------------------------------------------------------------
%     %accept/reject 1st component
%     
%     copula_data_income_star1 = w1_income_new*gamcdf(x,v1_income_new,m1_income_new/v1_income_new) + ...
%     w2_income*gamcdf(x,v2_income,m2_income/v2_income) + ...
%     w3_income*gamcdf(x,v3_income,m3_income/v3_income);
%         income_star_star1 = norminv(copula_data_income_star1,0,1);
%     id1_1 = copula_data_income_star1 <= 0;
%     copula_data_income_star1(id1_1,1) = 0.0001;
%     id1_2 = copula_data_income_star1 >= 0.9999;
%     copula_data_income_star1(id1_2,1) = 0.9999;    
%         A11 = rand();
%         %I = diag([1;1;1]);
%          log_lik1 = sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%              sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%          log_lik_star1 = sum(log(mvnpdf([income_star_star1 health_star educ_star],[0 0 0],cov_mat )))-...
%             sum(log(normpdf(income_star_star1,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%         
%         ratio_MH1 = exp(log_lik_star1 - log_lik1);
%     
%         MH1_copula = min(1,ratio_MH1);
%     
%     if A11<MH1_copula
%         w1_income = w1_income_new;
%         m1_income = m1_income_new;
%         v1_income = v1_income_new;
%         income_star = income_star_star1;
%         accept1=1;
%     else
%         w1_income = w1_income;
%         m1_income = m1_income;
%         v1_income = v1_income;
%         income_star = income_star;
%         accept1=0;
%     end
%     
%     %--------------------------------------------------------------------------------------------------------------
%     
% m2_income_new=1/(gamrnd(alpha2+N2*v2_income,1/(beta2+S2*v2_income)));
% 
% r2=600;
% 
% A2=rand(1);
% 
% V2_income=gamrnd(r2,v2_income/r2);
% 
% logvv2=v2_income*N2*log(v2_income)-N2*gammaln(v2_income)-v2_income*(theta2+(S2/m2_income_new)+N2*log(m2_income_new)-logP2);
% logvc2=V2_income*N2*log(V2_income)-N2*gammaln(V2_income)-V2_income*(theta2+(S2/m2_income_new)+N2*log(m2_income_new)-logP2);
% logV_v2=log(gampdf(v2_income,r2,V2_income/r2));
% logv_V2=log(gampdf(V2_income,r2,v2_income/r2));
% MH2=exp(logvc2+logV_v2-logvv2-logv_V2);
% C2=min(1,MH2);
%     if A2<=C2
%        v2_income_new=V2_income;
%     else
%        v2_income_new=v2_income;
%     end
%     
%     %--------------------------------------------------------------------------------------------------------------------------------
%     %accept/reject 2nd component
%     
%     copula_data_income_star2 = w1_income*gamcdf(x,v1_income,m1_income/v1_income) + ...
%     w2_income_new*gamcdf(x,v2_income_new,m2_income_new/v2_income_new) + ...
%     w3_income*gamcdf(x,v3_income,m3_income/v3_income);
%     id2_1 = copula_data_income_star2 <= 0;
%     copula_data_income_star2(id2_1,1) = 0.0001;
%     id2_2 = copula_data_income_star2 >= 0.9999;
%     copula_data_income_star2(id2_2,1) = 0.9999;    
%     income_star_star2 = norminv(copula_data_income_star2,0,1);
% 
%     A11 = rand();
%         %I = diag([1;1;1]);
%     log_lik2 = sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%              sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%     log_lik_star2 = sum(log(mvnpdf([income_star_star2 health_star educ_star],[0 0 0],cov_mat )))-...
%             sum(log(normpdf(income_star_star2,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%         
%         ratio_MH2 = exp(log_lik_star2 - log_lik2);
%     
%         MH2_copula = min(1,ratio_MH2);
%     
%     if A11<MH2_copula
%         w2_income = w2_income_new;
%         m2_income = m2_income_new;
%         v2_income = v2_income_new;
%         income_star = income_star_star2;
%         accept2=1;
%     else
%         w2_income = w2_income;
%         m2_income = m2_income;
%         v2_income = v2_income;
%         income_star = income_star;
%         accept2=0;
%     end
%     
%     %------------------------------------------------------------------------------------------------------------------------------
%     
%         m3_income_new=1/(gamrnd(alpha3+N3*v3_income,1/(beta3+S3*v3_income)));
%         r3=1000;
%         A3=rand(1);
%         V3_income=gamrnd(r3,v3_income/r3); 
%         logvv3=v3_income*N3*log(v3_income)-N3*gammaln(v3_income)-v3_income*(theta3+(S3/m3_income_new)+N3*log(m3_income_new)-logP3);
%         logvc3=V3_income*N3*log(V3_income)-N3*gammaln(V3_income)-V3_income*(theta3+(S3/m3_income_new)+N3*log(m3_income_new)-logP3);
%         logV_v3=log(gampdf(v3_income,r3,V3_income/r3));
%         logv_V3=log(gampdf(V3_income,r3,v3_income/r3));
%         MH3=exp(logvc3+logV_v3-logvv3-logv_V3);
%         C3=min(1,MH3);
%             if A3<=C3
%                 v3_income_new=V3_income;
%             else
%                 v3_income_new=v3_income;
%             end
%             
%             %--------------------------------------------------------------------------------------------------------------------------------
%     %accept/reject 3rd component
%     
%     copula_data_income_star3 = w1_income*gamcdf(x,v1_income,m1_income/v1_income) + ...
%     w2_income*gamcdf(x,v2_income,m2_income/v2_income) + ...
%     w3_income_new*gamcdf(x,v3_income_new,m3_income_new/v3_income_new);
%     
% 
%     id3_1 = copula_data_income_star3 <= 0;
%     copula_data_income_star3(id3_1,1) = 0.0001;
%     id3_2 = copula_data_income_star3 >=0.9999;
%     copula_data_income_star3(id3_2,1) = 0.999999;
%     income_star_star3 = norminv(copula_data_income_star3,0,1);
%     
%     A11 = rand();
%         %I = diag([1;1;1]);
%     log_lik3 = sum(log(mvnpdf([income_star health_star educ_star],[0 0 0],cov_mat ))) - ... 
%              sum(log(normpdf(income_star,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%     log_lik_star3 = sum(log(mvnpdf([income_star_star3 health_star educ_star],[0 0 0],cov_mat )))-...
%             sum(log(normpdf(income_star_star3,0,1))+log(normpdf(health_star,0,1))+log(normpdf(educ_star,0,1)));
%         
%         ratio_MH3 = exp(log_lik_star3 - log_lik3);
%     
%         MH3_copula = min(1,ratio_MH3);
%     
%     if A11<MH3_copula
%         w3_income = w3_income_new;
%         m3_income = m3_income_new;
%         v3_income = v3_income_new;
%         income_star = income_star_star3;
%         accept3=1;
%     else
%         w3_income = w3_income;
%         m3_income = m3_income;
%         v3_income = v3_income;
%         income_star = income_star;
%         accept3=0;
%     end
%     
%     %------------------------------------------------------------------------------------------------------------------------------
%     
%  