function [theta_gauss,scale_theta_copula,accept_theta_copula]=sampling_copula_parameter(data,theta_gauss,prior,iter,scale_theta_copula,...
                                                            accept_theta_copula,target_accept,income_tilde,health_tilde,educ_tilde,happiness_tilde,dim_cop)

                                                        

                                                        
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

y_cts = [income_tilde,health_tilde,educ_tilde,happiness_tilde];
loglik = sum(log(mvnpdf(y_cts,zeros(1,dim_cop),covmat)));

for k=1:length(theta_gauss)
    A1 = rand();
    R1 = mvnrnd(theta_gauss(k,1),scale_theta_copula(k));
    theta_gauss_star = theta_gauss;
    theta_gauss_star(k,1) = R1(1,1);
    chol_gauss_star=diag([ones(dim_cop,1)]);
    number_gauss = 1;
    
    for i=1:dim_cop-1
        for j=i:dim_cop-1
            chol_gauss_star(i,j+1) = theta_gauss_star(number_gauss,1); 
            number_gauss=number_gauss+1;
        end
    end
    sigma_star = chol_gauss_star'*chol_gauss_star;
    diag_mat_star = diag([1./sqrt(diag(sigma_star))]);
    covmat_star = diag_mat_star*sigma_star*diag_mat_star;  
    loglik_star = sum(log(mvnpdf(y_cts,zeros(1,dim_cop),covmat_star)));
    r1 = exp(loglik_star - loglik);
    C1(1,k) = min(1,r1);
    if A1<=C1(1,k)
       accept_theta_copula(k) = accept_theta_copula(k)+1; 
       theta_gauss = theta_gauss_star;
        
    end
    
    
    
end

if iter>1000
   for k=1:length(theta_gauss)  
   scale_theta_copula(1,k) = update_sigma(scale_theta_copula(1,k),C1(1,k),target_accept,iter,1);
   end
end



end

% function [theta12,theta13,theta23,accept_theta12,accept_theta13,accept_theta23] = gaussian_cop_tri(pseudo_income_data,pseudo_health_data,pseudo_educ_data,...
%              w1_income,w2_income,w3_income,m1_income,m2_income,m3_income,v1_income,v2_income,v3_income,...
%              w1_health,w2_health,w3_health,m1_health,m2_health,m3_health,s1_health,s2_health,s3_health,...           
%              p1_educ,p2_educ,p3_educ,p4_educ,p5_educ,p6_educ,p7_educ,p8_educ,p9_educ,p10_educ,p11_educ,p12_educ,p13_educ,...
%              income_star,health_star,educ_star,...
%              theta12,theta13,theta23)
% theta_gaussian = [1 theta12 theta13; 0 1 theta23; 0 0 1];
% sig = (theta_gaussian'*theta_gaussian);
% diag_mat = [(1/sqrt(sig(1,1))) 0 0 ; 0 (1/sqrt(sig(2,2))) 0 ; 0 0 (1/sqrt(sig(3,3)))];
% cov_mat = diag_mat*sig*diag_mat;
%  twice_trans_data = [income_star health_star educ_star];
%  u1 = rand();
%  u2 = rand();
%  u3 = rand();
% % 
% %MH for theta12
% log_lik = sum(log(mvnpdf(twice_trans_data,[0 0 0],cov_mat)));
% theta12_star = normrnd(theta12,0.005);
% 
% rstar1 = [1 theta12_star theta13; 0 1 theta23; 0 0 1];
% sig_star1 = (rstar1'*rstar1);
% diag_mat_star1 = [(1/sqrt(sig_star1(1,1))) 0 0 ; 0 (1/sqrt(sig_star1(2,2))) 0 ; 0 0 (1/sqrt(sig_star1(3,3)))];
% cov_mat_star1 = diag_mat_star1*sig_star1*diag_mat_star1;
% log_lik1_star = sum(log(mvnpdf(twice_trans_data,[0 0 0],cov_mat_star1)));
% sum(log(mvnpdf(y_cts,zeros(1,dim_cop),covmat)));
%     
%     MH1 = min(1,ratio1);
%     
%     if u1<MH1
%         theta12 = theta12_star;
%         accept_theta12 = 1;
%     else
%         theta12 = theta12;
%         accept_theta12 = 0;
%     end
% 
% 
% %MH for theta13
% theta13_star = normrnd(theta13,0.005);
% rstar2 = [1 theta12 theta13_star;0 1 theta23; 0 0 1];
% sig_star2 = (rstar2'*rstar2);
% diag_mat_star2 = [(1/sqrt(sig_star2(1,1))) 0 0 ; 0 (1/sqrt(sig_star2(2,2))) 0 ; 0 0 (1/sqrt(sig_star2(3,3)))];
% cov_mat_star2 = diag_mat_star2*sig_star2*diag_mat_star2;
% log_lik2_star = sum(log(mvnpdf(twice_trans_data,[0 0 0],cov_mat_star2)));
% ratio2 = exp(log_lik2_star - log_lik);
% MH2 = min(1,ratio2);
% if u2<MH2
%         theta13 = theta13_star;
%         accept_theta13 = 1;
%     else
%         theta13 = theta13;
%         accept_theta13 = 0;
% end
% 
% %MH for theta23
% theta23_star = normrnd(theta23,0.005);
% rstar3 = [1 theta12 theta13;0 1 theta23_star; 0 0 1];
% sig_star3 = (rstar3'*rstar3);
% diag_mat_star3 = [(1/sqrt(sig_star3(1,1))) 0 0 ; 0 (1/sqrt(sig_star3(2,2))) 0 ; 0 0 (1/sqrt(sig_star3(3,3)))];
% cov_mat_star3 = diag_mat_star3*sig_star3*diag_mat_star3;
% log_lik3_star = sum(log(mvnpdf(twice_trans_data,[0 0 0],cov_mat_star3)));
% ratio3 = exp(log_lik3_star - log_lik);
% MH3 = min(1,ratio3);
% if u3<MH3
%         theta23 = theta23_star;
%         accept_theta23 = 1;
%     else
%         theta23 = theta23;
%         accept_theta23 = 0;







% 
% r1 = [1 r12 r13 ; 0 1 r23 ; 0 0 1];
%     sig1 = (r1'*r1);
%     diag_mat1 = [(1/sqrt(sig1(1,1))) 0 0 ; 0 (1/sqrt(sig1(2,2))) 0 ; 0 0 (1/sqrt(sig1(3,3)))];
%     covmat1 = diag_mat1*sig1*diag_mat1;     
%     log_lik1 = sum(log(mvnpdf(x,[0 0 0],covmat1)))-sum(log(normpdf(xstar.^2,0,1))+log(normpdf(ystar.^2,0,1))+log(normpdf(zstar.^2,0,1)));
%     r12_star = normrnd(r12,0.01);
% %     theta1 = [r12 r13 r23];
% %     k1 = 0.5;
% %     R1 = mvnrnd(theta1,k1.*(cov_dens'*cov_dens));
% %     r12_star = R1(1,1);
% %     r13_star = R1(1,2);
% %     r23_star = R1(1,3);
%     
%     %r_star = [1 r12_star r13_star; 0 1 r23_star ; 0 0 1];
%      rstar1 = [1 r12_star r13;0 1 r23; 0 0 1];
%      sig_star1 = (rstar1'*rstar1);
%      diag_mat_star1 = [(1/sqrt(sig_star1(1,1))) 0 0 ; 0 (1/sqrt(sig_star1(2,2))) 0 ; 0 0 (1/sqrt(sig_star1(3,3)))];
%      cov_mat_star1 = diag_mat_star1*sig_star1*diag_mat_star1;
%      log_lik1_star = sum(log(mvnpdf(x,[0 0 0],covmat_star1)))-sum(log(normpdf(xstar.^2,0,1))+log(normpdf(ystar.^2,0,1))+log(normpdf(zstar.^2,0,1)));
%     
%     ratio1 = exp(log_lik1_star - log_lik1);
%     
%     MH1 = min(1,ratio1);
%     
%     if u1<MH1
%         r12 = r12_star;
%         accept1 = accept1+1;
%     else
%         r12 = r12;
%     end














%      A1_gaussian = rand(1);
%      prior_gaussian1 = log(unifpdf(theta_gaussian,-1,1));
%      lik_gaussian1 =  sum(log(normal_copula_pdf(copula_data_income,copula_data_health,theta_gaussian)));
%      p1_gaussian1 = prior_gaussian1+lik_gaussian1;
%      theta_gaussian_star = normrnd(theta_gaussian,0.01);
%      prior_gaussian_star1 = log(unifpdf(theta_gaussian_star,-1,1));
%      lik_gaussian_star1 =  sum(log(normal_copula_pdf(copula_data_income,copula_data_health,theta_gaussian_star)));
%      p1_gaussian_star1 = prior_gaussian_star1+lik_gaussian_star1;
%      r_gaussian1 = exp(p1_gaussian_star1 - p1_gaussian1);
%      C1_gaussian1 = min(1,r_gaussian1);
%      if A1_gaussian<=C1_gaussian1
%         theta_gaussian = theta_gaussian_star;
%         
%      else
%          theta_gaussian = theta_gaussian;
%      end

% r = [1 r12 r13 ; 0 1 r23 ; 0 0 1];
%     sig = (r'*r);
%     diag_mat = [(1/sqrt(sig(1,1))) 0 0 ; 0 (1/sqrt(sig(2,2))) 0 ; 0 0 (1/sqrt(sig(3,3)))];
%     covmat = diag_mat*sig*diag_mat;

% if educ_category_2006(j,1) == 1
%         mu = 0 + [covmat_used(3,1) covmat_used(3,2)]*(inv([1 covmat_used(1,2); covmat_used(2,1) 1]))*[xstar(j,1);ystar(j,1)];
%         cov = sqrt(1 - [covmat_used(3,1) covmat_used(3,2)]*(inv([1 covmat_used(1,2);covmat_used(2,1) 1]))*[covmat_used(1,3);covmat_used(2,3)]);
%         zstar(j,1) = truncated_normal(mu,cov,-inf,latent1_educ,1);
%         
% end
