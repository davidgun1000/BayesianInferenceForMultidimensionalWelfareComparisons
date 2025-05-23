%copula_main_function
%function out = copula_main_function_tri(k_P)
% seednum = str2double(k_P);
% rng(seednum);
function [prior,Post,V1_gamma_com]=InitialSamplerGammaMix(data,num_component_gamma)
prior.phi=1;
prior.theta_gamma=0.01;
prior.alpha_gamma=10;
prior.beta_gamma=[20,50,100,200];
prior.c1=1000;
burn = 5000;
nit = 10000;
s=burn+nit;

%get initial values

sort_income = sort(data(:,1));
for k=1:num_component_gamma
    if k==1
       X = sort_income(1:round((k/num_component_gamma)*length(sort_income(:,1))));         
    else
       X = sort_income(round(((k-1)/num_component_gamma)*length(sort_income(:,1))):round(((k)/num_component_gamma)*length(sort_income(:,1))));  
    end
    phat(k,:) = mle(X,'distribution','gamma');
end
m_gamma = phat(:,1)'.*phat(:,2)';
%v_gamma = phat(:,2)';

w_gamma = ones(1,num_component_gamma)./num_component_gamma;
%m_gamma = [50,250,500];
v_gamma = 2*ones(1,num_component_gamma);

for j=1:num_component_gamma
    V1_gamma_com(:,:,j)=0.01*eye(2);
    scale_gamma_com(1,j)=0.1;
end

accept_w_gamma=0;
accept_gamma_com=zeros(num_component_gamma,1);
target_accept=0.2;

for i=1:s
    %i
    %m_gamma'
    [w_gamma,m_gamma,v_gamma,scale_gamma_com,accept_w_gamma,accept_gamma_com]=sampling_gamma_parameter_only_Mix(data,w_gamma,m_gamma,v_gamma,prior,i,scale_gamma_com,V1_gamma_com,...
                                                                                                       accept_w_gamma,accept_gamma_com,num_component_gamma,target_accept);
    
    for k=1:num_component_gamma
        thetasave_gamma{k,1}(i,:) = [log(m_gamma(k)),log(v_gamma(k))];
    end
    
    if i>1000
        for k=1:num_component_gamma
            V1_gamma_com(:,:,k) = cov(thetasave_gamma{k,1}(1:i,:));
            V1_gamma_com(:,:,k) = jitChol(V1_gamma_com(:,:,k));
        end
    end
    temp.w_gamma(i,:) = w_gamma;
    temp.m_gamma(i,:) = m_gamma;
    temp.v_gamma(i,:) = v_gamma;
    
    
                                                                    
                                                                    
end

Post.w_gamma = mean(temp.w_gamma(burn+1:end,:));
Post.m_gamma = mean(temp.m_gamma(burn+1:end,:));
Post.v_gamma = mean(temp.v_gamma(burn+1:end,:));

end

% for k=1:num_category
%     prop_init(1,k)=mean(data(:,2)==k);
% end
% 
% cum_marginal = cumsum(prop_init(1:num_category-1));
% limit = norminv(cum_marginal);
% for j=1:num_category-1
%     if j==1
%         delta(1,j)=log(limit(j)+4);
%     else
%         delta(1,j)=log(limit(j)-limit(j-1));
%     end
% end



% load('overall_2001_income_health.mat');
% 
% income_data = real_atifdip;
% weight_data = ahhwtrps;
% health_data = aghmh_scaled;
% educ_data = years_aeduc;
% n=size(income_data,1);
% 
% b=1000;	%number of burn-in replications
% nit=10000;	%number of retained replications / sample size
% tot=b+nit;	%total
% k=3;
% W1_income=zeros(nit,1);
% W2_income=zeros(nit,1);
% W3_income=zeros(nit,1);
% M1_income=zeros(nit,1);
% M2_income=zeros(nit,1);
% M3_income=zeros(nit,1);
% V1_income=zeros(nit,1);
% V2_income=zeros(nit,1);
% V3_income=zeros(nit,1);
% 
% W1_health=zeros(nit,1);
% W2_health=zeros(nit,1);
% W3_health=zeros(nit,1);
% M1_health=zeros(nit,1);
% M2_health=zeros(nit,1);
% M3_health=zeros(nit,1);
% S1_health=zeros(nit,1);
% S2_health=zeros(nit,1);
% S3_health=zeros(nit,1);
% 
% p1_educ_store = zeros(nit,1);
% p2_educ_store = zeros(nit,1);
% p3_educ_store = zeros(nit,1);
% p4_educ_store = zeros(nit,1);
% p5_educ_store = zeros(nit,1);
% p6_educ_store = zeros(nit,1);
% p7_educ_store = zeros(nit,1);
% p8_educ_store = zeros(nit,1);
% p9_educ_store = zeros(nit,1);
% p10_educ_store = zeros(nit,1);
% p11_educ_store = zeros(nit,1);
% p12_educ_store = zeros(nit,1);
% p13_educ_store = zeros(nit,1);
% 
% 
% 
% theta12_store = zeros(nit,1);
% theta13_store = zeros(nit,1);
% theta23_store = zeros(nit,1);
% rho12_store = zeros(nit,1);
% rho13_store = zeros(nit,1);
% rho23_store = zeros(nit,1);
% com = [income_data health_data educ_data weight_data];
% com(:,5) = com(:,4)./sum(com(:,4));
% com(:,6) = cumsum(com(:,5));
% com(n,6) = 1;
% 
% %----------------------------------------------------------------------------------------------------------------------
% %initialisation
%    w1_income = 0.13;
%    w2_income = 0.42;
%    w3_income = 1-w1_income-w2_income;
%    m1_income =  128.1318;
%    m2_income =  286.3545;
%    m3_income = 389.7359;
%    v1_income =  33.6899;
%    v2_income =  0.86;
%    v3_income = 5.39;
%    
%    w1_health = 0.0762;
%    w2_health = 0.5742;
%    w3_health = 1-w1_health-w2_health;
%    m1_health =  0.7525;
%    m2_health =  0.6721;
%    m3_health = 0.8527;
%    s1_health =  1.2818;
%    s2_health =  8.4229;
%    s3_health = 27.4908;
%    
%    theta12 = 0.0970;
%    theta13 = 0.3392;
%    theta23 = 0.0834;
%    %--------------------------------------------------------------------------
%    %educ
% %    p1_educ = 0.0029;
% %    p2_educ = 0.0101;
% %    p3_educ = 0.0107;
% %    p4_educ  = 0.0201;
% %    p5_educ = 0.0373;
% %    p6_educ = 0.1300;
% %    p7_educ = 0.0887;
% %    p8_educ = 0.1705;
% %    p9_educ = 0.2110;
% %    p10_educ = 0.0895;
% %    p11_educ = 0.1339;
% %    p12_educ = 0.0553;
% %    p13_educ = 0.0401;
% %    
% %    limit_p1 = norminv(p1_educ,0,1);
% %   limit_p2 = norminv(p1_educ+p2_educ,0,1);
% %   limit_p3 = norminv(p1_educ+p2_educ+p3_educ,0,1);
% %   limit_p4 = norminv(p1_educ+p2_educ+p3_educ+p4_educ,0,1);
% %   limit_p5 = norminv(p1_educ+p2_educ+p3_educ+p4_educ+p5_educ,0,1);
% %   limit_p6 = norminv(p1_educ+p2_educ+p3_educ+p4_educ+p5_educ+p6_educ,0,1);
% %   limit_p7 = norminv(p1_educ+p2_educ+p3_educ+p4_educ+p5_educ+p6_educ+p7_educ,0,1);
% %   limit_p8 = norminv(p1_educ+p2_educ+p3_educ+p4_educ+p5_educ+p6_educ+p7_educ+p8_educ,0,1);
% %   limit_p9 = norminv(p1_educ+p2_educ+p3_educ+p4_educ+p5_educ+p6_educ+p7_educ+p8_educ+p9_educ,0,1);
% %   limit_p10 = norminv(p1_educ+p2_educ+p3_educ+p4_educ+p5_educ+p6_educ+p7_educ+p8_educ+p9_educ+p10_educ,0,1);
% %   limit_p11 = norminv(p1_educ+p2_educ+p3_educ+p4_educ+p5_educ+p6_educ+p7_educ+p8_educ+p9_educ+p10_educ+p11_educ,0,1);
% %   limit_p12 = norminv(p1_educ+p2_educ+p3_educ+p4_educ+p5_educ+p6_educ+p7_educ+p8_educ+p9_educ+p10_educ+p11_educ+p12_educ,0,1);
% %   limit_p13 = norminv(0.9999,0,1);
%   
%   copula_data_income_init = w1_income*gamcdf(income_data,v1_income,m1_income/v1_income) + ...
%     w2_income*gamcdf(income_data,v2_income,m2_income/v2_income) + ...
%     w3_income*gamcdf(income_data,v3_income,m3_income/v3_income);
%    income_star = norminv(copula_data_income_init,0,1);
%    
%    copula_data_health_init = w1_health*betacdf(health_data,s1_health*m1_health,s1_health*(1-m1_health)) + ...
%     w2_health*betacdf(health_data,s2_health*m2_health,s2_health*(1-m2_health)) + ...
%     w3_health*betacdf(health_data,s3_health*m3_health,s3_health*(1-m3_health));
%    health_star = norminv(copula_data_health_init,0,1);
%    
%    theta_gaussian = [1 theta12 theta13; 0 1 theta23; 0 0 1];
%    sig = (theta_gaussian'*theta_gaussian);
%    diag_mat = [(1/sqrt(sig(1,1))) 0 0 ; 0 (1/sqrt(sig(2,2))) 0 ; 0 0 (1/sqrt(sig(3,3)))];
%    cov_mat = diag_mat*sig*diag_mat;
%    cov =  sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
% mu =  [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star';health_star'];
% mu = mu';
% var_prior_delta = 10^2*eye(12);
% mu_prior_delta = [0.2148 -0.6251 -1.4164 -1.2922 -1.1707 -0.5190 -1.2814 -0.7975 -0.6055 -1.3091 -0.5669 -0.8178]; 
% 
% 
% 
% [estimator,cov_educ] = copula_ml('ordered',educ_data);
% 
% delta1 = estimator(1,1);
% delta2 = estimator(1,2);
% delta3 = estimator(1,3);
% delta4 = estimator(1,4);
% delta5 = estimator(1,5);
% delta6 = estimator(1,6);
% delta7 = estimator(1,7);
% delta8 = estimator(1,8);
% delta9 = estimator(1,9);
% delta10 = estimator(1,10);
% delta11 = estimator(1,11);
% delta12 = estimator(1,12);
% sig = chol(cov_educ);
% prior1 = log(mvnpdf([delta1 delta2 delta3 delta4 delta5 delta6 delta7 delta8 delta9 delta10 delta11 delta12],mu_prior_delta,var_prior_delta));
%   limit_p1 = -4 + exp(delta1);
%   limit_p2 = -4 + exp(delta1)+ exp(delta2);
%   limit_p3 = -4 + exp(delta1)+ exp(delta2) + exp(delta3);
%   limit_p4 = -4 + exp(delta1)+ exp(delta2) + exp(delta3) + exp(delta4);
%   limit_p5 = -4 + exp(delta1)+ exp(delta2) + exp(delta3) + exp(delta4) + exp(delta5);
%   limit_p6 = -4 + exp(delta1)+ exp(delta2) + exp(delta3) + exp(delta4) + exp(delta5) + exp(delta6);
%   limit_p7 = -4 + exp(delta1)+ exp(delta2) + exp(delta3) + exp(delta4) + exp(delta5) + exp(delta6) + exp(delta7);
%   limit_p8 = -4 + exp(delta1)+ exp(delta2) + exp(delta3) + exp(delta4) + exp(delta5) + exp(delta6) + exp(delta7) + exp(delta8);
%   limit_p9 = -4 + exp(delta1)+ exp(delta2) + exp(delta3) + exp(delta4) + exp(delta5) + exp(delta6) + exp(delta7) + exp(delta8) + exp(delta9);
%   limit_p10 = -4 + exp(delta1)+ exp(delta2) + exp(delta3) + exp(delta4) + exp(delta5) + exp(delta6) + exp(delta7) + exp(delta8) + exp(delta9) + exp(delta10);
%   limit_p11 = -4 + exp(delta1)+ exp(delta2) + exp(delta3) + exp(delta4) + exp(delta5) + exp(delta6) + exp(delta7) + exp(delta8) + exp(delta9) + exp(delta10) + exp(delta11);
%   limit_p12 = -4 + exp(delta1)+ exp(delta2) + exp(delta3) + exp(delta4) + exp(delta5) + exp(delta6) + exp(delta7) + exp(delta8) + exp(delta9) + exp(delta10) + exp(delta11) + exp(delta12);
%   limit_p13 = norminv(0.9999,0,1);
%    for s = 1:n
%     
%      if educ_data(s,1) == 3
%         %mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
%         %cov(s,1) = sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%         educ_star(s,1) = truncated_normal(mu(s,1),cov,-20,limit_p1,1); 
%     end
%     
%     if educ_data(s,1) == 6
%         %mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
%         %cov(s,1) = sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%         educ_star(s,1) = truncated_normal(mu(s,1),cov,limit_p1,limit_p2,1); 
%     end
%     
%     if educ_data(s,1) == 7
%         %mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
%         %cov(s,1) = sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%         educ_star(s,1) = truncated_normal(mu(s,1),cov,limit_p2,limit_p3,1); 
%     end
%     
%     if educ_data(s,1) == 8
%         %mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
%         %cov(s,1) = sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%         educ_star(s,1) = truncated_normal(mu(s,1),cov,limit_p3,limit_p4,1); 
%     end
%     
%     if educ_data(s,1) == 9
%         %mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
%         %cov(s,1) = sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%         educ_star(s,1) = truncated_normal(mu(s,1),cov,limit_p4,limit_p5,1); 
%     end
%     
%     if educ_data(s,1) == 10
%         %mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
%         %cov(s,1) = sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%         educ_star(s,1) = truncated_normal(mu(s,1),cov,limit_p5,limit_p6,1);
%     end
%     
%     if educ_data(s,1) == 11
%         %mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
%         %cov(s,1) = sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%         educ_star(s,1) = truncated_normal(mu(s,1),cov,limit_p6,limit_p7,1); 
%     end
%     
%     if educ_data(s,1) == 12
%         %mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
%         %cov(s,1) = sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%         educ_star(s,1) = truncated_normal(mu(s,1),cov,limit_p7,limit_p8,1); 
%     end
%     
%     if educ_data(s,1) == 13
%         %mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
%         %cov(s,1) = sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%         educ_star(s,1) = truncated_normal(mu(s,1),cov,limit_p8,limit_p9,1);
%     end
%     
%     if educ_data(s,1) == 14
%         %mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
%         %cov(s,1) = sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%         educ_star(s,1) = truncated_normal(mu(s,1),cov,limit_p9,limit_p10,1); 
%     end
%     
%     if educ_data(s,1) == 16
%         %mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
%        % cov(s,1) = sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%         educ_star(s,1) = truncated_normal(mu(s,1),cov,limit_p10,limit_p11,1);
%     end
%     
%     if educ_data(s,1) == 17
% %         mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
% %         cov(s,1) = sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%          educ_star(s,1) = truncated_normal(mu(s,1),cov,limit_p11,limit_p12,1); 
%     end
%     
%     if educ_data(s,1) == 20
% %         mu(s,1) = 0 + [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2); cov_mat(2,1) 1]))*[income_star(s,1);health_star(s,1)];
% %         cov(s,1)= sqrt(1 - [cov_mat(3,1) cov_mat(3,2)]*(inv([1 cov_mat(1,2);cov_mat(2,1) 1]))*[cov_mat(1,3);cov_mat(2,3)]);
%          educ_star(s,1) = truncated_normal(mu(s,1),cov,limit_p12,limit_p13,1);
%     end
%     
%    end
%    %log_lik_educ = sum(log_lik_educ_indv);
%    %p1= log_lik_educ+prior1;
%    %------------------------------------------------------------------------------------------------------------------------
%    count0_income = 0;
%    count1_income = 0;
%    count2_income = 0;
%    count3_income = 0;
%    count1_health = 0;
%    count2_health = 0;
%    count3_health = 0;
%    count_educ = 0;
%    count_theta12 = 0;
%    count_theta13 = 0;
%    count_theta23 = 0;
%    
% for i=1:tot;
% %    w1_income = 0.16;
% %    w2_income = 0.4415;
% %    w3_income = 1-w1_income-w2_income;
% %    m1_income =  166.1318;
% %    m2_income =  408.3545;
% %    m3_income = 458.7359;
% %    v1_income =  12.6899;
% %    v2_income =  0.8397;
% %    v3_income = 7.3662;
% 
%    for j = 1:n
%      u1 = rand();
%      index = find(u1 <= com(:,6),1,'first');
%      weighted_data_income(j,1) = com(index(1,1),1);    
%      weighted_data_health(j,1) = com(index(1,1),2);
%      weighted_data_educ(j,1) = com(index(1,1),3);
%    end
%    
%      pseudo_income_data = weighted_data_income;
%      pseudo_health_data = weighted_data_health;
%      pseudo_educ_data = weighted_data_educ;
%    
%       copula_data_income = w1_income*gamcdf(pseudo_income_data,v1_income,m1_income/v1_income) + ...
%     w2_income*gamcdf(pseudo_income_data,v2_income,m2_income/v2_income) + ...
%     w3_income*gamcdf(pseudo_income_data,v3_income,m3_income/v3_income);
%         income_star = norminv(copula_data_income,0,1);
%         
%       copula_data_health = w1_health*betacdf(pseudo_health_data,s1_health*m1_health,s1_health*(1-m1_health)) + ...
%     w2_health*betacdf(pseudo_health_data,s2_health*m2_health,s2_health*(1-m2_health)) + ...
%     w3_health*betacdf(pseudo_health_data,s3_health*m3_health,s3_health*(1-m3_health));
%        health_star = norminv(copula_data_health,0,1);
%      
%        [educ_star]  = get_latent_educ(pseudo_educ_data,delta1,delta2,delta3,delta4,delta5,delta6,delta7,delta8,delta9,delta10,delta11,delta12,...
%                                             theta12,theta13,theta23,income_star,health_star);
%       
% 
%     
%          
%           [w1_income,w2_income,w3_income,m1_income,m2_income,m3_income,v1_income,v2_income,v3_income,accept0_income,accept1_income,accept2_income,income_star] = gamma3_mixture_draws_tri_block_new(pseudo_income_data,...
%              w1_income,w2_income,w3_income,m1_income,m2_income,m3_income,v1_income,v2_income,v3_income,...
%              health_star,educ_star,...
%              theta12,theta13,theta23);
%          count0_income = count0_income+accept0_income;
%          count1_income = count1_income+accept1_income;
%          count2_income = count2_income+accept2_income;
%          %count3_income = count3_income+accept3_income;
%          
%          [w1_health,w2_health,w3_health,m1_health,m2_health,m3_health,s1_health,s2_health,s3_health,accept1_health,accept2_health,accept3_health,health_star] = beta3_mixture_draws_tri_block_new(pseudo_health_data,...
%              w1_health,w2_health,w3_health,m1_health,m2_health,m3_health,s1_health,s2_health,s3_health,...
%              income_star,educ_star,...
%              theta12,theta13,theta23);
%          count1_health = count1_health+accept1_health;
%          count2_health = count2_health+accept2_health;
%          count3_health = count3_health+accept3_health;
%          
%          [p1_educ,p2_educ,p3_educ,p4_educ,p5_educ,p6_educ,p7_educ,p8_educ,p9_educ,p10_educ,p11_educ,p12_educ,p13_educ,...
%           delta1,delta2,delta3,delta4,delta5,delta6,delta7,delta8,delta9,delta10,delta11,delta12,...
%           accept_educ,educ_star] = educ_draws_discrete_alternative3(pseudo_educ_data,...
%             delta1,delta2,delta3,delta4,delta5,delta6,delta7,delta8,delta9,delta10,delta11,delta12,...
%             income_star,health_star,...
%             theta12,theta13,theta23,...
%             cov_educ); 
%         
%          count_educ = count_educ+accept_educ;   
%         
%   [theta12,theta13,theta23,accept_theta12,accept_theta13,accept_theta23] = gaussian_cop_tri(pseudo_income_data,pseudo_health_data,pseudo_educ_data,...
%              w1_income,w2_income,w3_income,m1_income,m2_income,m3_income,v1_income,v2_income,v3_income,...
%              w1_health,w2_health,w3_health,m1_health,m2_health,m3_health,s1_health,s2_health,s3_health,...
%              p1_educ,p2_educ,p3_educ,p4_educ,p5_educ,p6_educ,p7_educ,p8_educ,p9_educ,p10_educ,p11_educ,p12_educ,p13_educ,...
%              income_star,health_star,educ_star,...
%              theta12,theta13,theta23);
%         
%          count_theta12 = count_theta12+accept_theta12;
%          count_theta13 = count_theta13+accept_theta13;
%          count_theta23 = count_theta23+accept_theta23;
% 
% 
%          
%      
%          est_theta_gaussian = [1 theta12 theta13; 0 1 theta23; 0 0 1];
%          est_sig = (est_theta_gaussian'*est_theta_gaussian);
%          est_diag_mat = [(1/sqrt(est_sig(1,1))) 0 0 ; 0 (1/sqrt(est_sig(2,2))) 0 ; 0 0 (1/sqrt(est_sig(3,3)))];
%          est_cov_mat = est_diag_mat*est_sig*est_diag_mat;
%      if m1_income>m3_income
%         qq_income=m3_income;
%         m3_income=m1_income;
%         m1_income=qq_income;
%         qq_income=v3_income;
%         v3_income=v1_income;
%         v1_income=qq_income;
%         qq_income=w3_income;
%         w3_income=w1_income;
%         w1_income=qq_income;
%      end
%     if m1_income>m2_income
%         qq_income=m2_income;
%         m2_income=m1_income;
%         m1_income=qq_income;
%         qq_income=v2_income;
%         v2_income=v1_income;
%         v1_income=qq_income;
%         qq_income=w2_income;
%         w2_income=w1_income;
%         w1_income=qq_income;
%     end
%     if s2_health>s3_health
%         qq_health=m3_health;
%         m3_health=m2_health;
%         m2_health=qq_health;
%         qq_health=s3_health;
%         s3_health=s2_health;
%         s2_health=qq_health;
%         qq_health=w3_health;
%         w3_health=w2_health;
%         w2_health=qq_health;
%     end
% 
%     if s1_health>s3_health
%         qq_health=m3_health;
%         m3_health=m1_health;
%         m1_health=qq_health;
%         qq_health=s3_health;
%         s3_health=s1_health;
%         s1_health=qq_health;
%         qq_health=w3_health;
%         w3_health=w1_health;
%         w1_health=qq_health;
%     end
%     if s1_health>s2_health
%         qq_health=m2_health;
%         m2_health=m1_health;
%         m1_health=qq_health;
%         qq_health=s2_health;
%         s2_health=s1_health;
%         s1_health=qq_health;
%         qq_health=w2_health;
%         w2_health=w1_health;
%         w1_health=qq_health;
%     end   
%     
%     
%        if i>b
%       W1_income(i-b,1)=w1_income;
%       W2_income(i-b,1)=w2_income;
%       W3_income(i-b,1)=w3_income;
%       M1_income(i-b,1)=m1_income;
%       M2_income(i-b,1)=m2_income;
%       M3_income(i-b,1)=m3_income;
%       V1_income(i-b,1)=v1_income;
%       V2_income(i-b,1)=v2_income;
%       V3_income(i-b,1)=v3_income;
% %       
%       W1_health(i-b,1)=w1_health;
%       W2_health(i-b,1)=w2_health;
%       W3_health(i-b,1)=w3_health;
%       M1_health(i-b,1)=m1_health;
%       M2_health(i-b,1)=m2_health;
%       M3_health(i-b,1)=m3_health;
%       S1_health(i-b,1)=s1_health;
%       S2_health(i-b,1)=s2_health;
%       S3_health(i-b,1)=s3_health;
%       
%       theta12_store(i-b,1) = theta12;
%       theta13_store(i-b,1) = theta13;
%       theta23_store(i-b,1) = theta23;
%       rho12_store(i-b,1) = est_cov_mat(1,2);
%       rho13_store(i-b,1) = est_cov_mat(1,3);
%       rho23_store(i-b,1) = est_cov_mat(2,3);
%       
%       p1_educ_store(i-b,1) = p1_educ;
%       p2_educ_store(i-b,1) = p2_educ;
%       p3_educ_store(i-b,1) = p3_educ;
%       p4_educ_store(i-b,1) = p4_educ;
%       p5_educ_store(i-b,1) = p5_educ;
%       p6_educ_store(i-b,1) = p6_educ;
%       p7_educ_store(i-b,1) = p7_educ;
%       p8_educ_store(i-b,1) = p8_educ;
%       p9_educ_store(i-b,1) = p9_educ;
%       p10_educ_store(i-b,1) = p10_educ;
%       p11_educ_store(i-b,1) = p11_educ;
%       p12_educ_store(i-b,1) = p12_educ;
%       p13_educ_store(i-b,1) = p13_educ;
%     end        
%   %  nm = ['output_copula_2001_city',k_P,'.mat']
%   %  save(nm,'W1_income','W2_income','W3_income',...
% %     'M1_income','M2_income','M3_income',...
% %     'V1_income','V2_income','V3_income',...
% %     'W1_health','W2_health','W3_health',...
% %     'M1_health','M2_health','M3_health',...
% %     'S1_health','S2_health','S3_health',...
% %     'p1_educ_store','p2_educ_store','p3_educ_store','p4_educ_store',...
% %     'p5_educ_store','p6_educ_store','p7_educ_store','p8_educ_store',...
% %     'p9_educ_store','p10_educ_store','p11_educ_store','p12_educ_store',...
% %     'p13_educ_store',...
% %     'theta12_store','theta13_store','theta23_store',...
% %     'rho12_store','rho13_store','rho23_store');
% 
% 
% W1_income_overall_2001 = W1_income(:,1);
% W2_income_overall_2001 = W2_income(:,1);
% W3_income_overall_2001 = W3_income(:,1);
% M1_income_overall_2001 = M1_income(:,1);
% M2_income_overall_2001 = M2_income(:,1);
% M3_income_overall_2001 = M3_income(:,1);
% V1_income_overall_2001 = V1_income(:,1);
% V2_income_overall_2001 = V2_income(:,1);
% V3_income_overall_2001 = V3_income(:,1);
%     
% W1_health_overall_2001 = W1_health(:,1);
% W2_health_overall_2001 = W2_health(:,1);
% W3_health_overall_2001 = W3_health(:,1);
% M1_health_overall_2001 = M1_health(:,1);
% M2_health_overall_2001 = M2_health(:,1);
% M3_health_overall_2001 = M3_health(:,1);
% S1_health_overall_2001 = S1_health(:,1);
% S2_health_overall_2001 = S2_health(:,1);
% S3_health_overall_2001 = S3_health(:,1);
% 
% p1_educ_overall_2001 = p1_educ_store(:,1);
% p2_educ_overall_2001 = p2_educ_store(:,1);
% p3_educ_overall_2001 = p3_educ_store(:,1);
% p4_educ_overall_2001 = p4_educ_store(:,1);
% p5_educ_overall_2001 = p5_educ_store(:,1);
% p6_educ_overall_2001 = p6_educ_store(:,1);
% p7_educ_overall_2001 = p7_educ_store(:,1);
% p8_educ_overall_2001 = p8_educ_store(:,1);
% p9_educ_overall_2001 = p9_educ_store(:,1);
% p10_educ_overall_2001 = p10_educ_store(:,1);
% p11_educ_overall_2001 = p11_educ_store(:,1);
% p12_educ_overall_2001 = p12_educ_store(:,1);
% p13_educ_overall_2001 = p13_educ_store(:,1);
% % 
% theta12_overall_2001 = theta12_store(:,1);
% theta13_overall_2001 = theta13_store(:,1);
% theta23_overall_2001 = theta23_store(:,1);
% rho12_overall_2001 = rho12_store(:,1);
% rho13_overall_2001 = rho13_store(:,1);
% rho23_overall_2001 = rho23_store(:,1);
% 
% 
% % nm = ['output_copula_overall',k_P,'.mat']
% %    save(nm,'W1_income_overall_2001','W2_income_overall_2001','W3_income_overall_2001',...
% %     'M1_income_overall_2001','M2_income_overall_2001','M3_income_overall_2001',...
% %     'V1_income_overall_2001','V2_income_overall_2001','V3_income_overall_2001',...
% %     'W1_health_overall_2001','W2_health_overall_2001','W3_health_overall_2001',...
% %     'M1_health_overall_2001','M2_health_overall_2001','M3_health_overall_2001',...
% %     'S1_health_overall_2001','S2_health_overall_2001','S3_health_overall_2001',...
% %     'p1_educ_overall_2001','p2_educ_overall_2001','p3_educ_overall_2001','p4_educ_overall_2001',...
% %     'p5_educ_overall_2001','p6_educ_overall_2001','p7_educ_overall_2001','p8_educ_overall_2001',...
% %     'p9_educ_overall_2001','p10_educ_overall_2001','p11_educ_overall_2001','p12_educ_overall_2001',...
% %     'p13_educ_overall_2001',...
% %     'theta12_overall_2001','theta13_overall_2001','theta23_overall_2001',...
% %     'rho12_overall_2001','rho13_overall_2001','rho23_overall_2001');
% 
% end
% %    save('copula_overall2001_without_weight.mat','W1_income_overall_2010','W2_income_overall_2010','W3_income_overall_2010',...
% %     'M1_income_overall_2010','M2_income_overall_2010','M3_income_overall_2010',...
% %     'V1_income_overall_2010','V2_income_overall_2010','V3_income_overall_2010',...
% %     'W1_health_overall_2010','W2_health_overall_2010','W3_health_overall_2010',...
% %     'M1_health_overall_2010','M2_health_overall_2010','M3_health_overall_2010',...
% %     'S1_health_overall_2010','S2_health_overall_2010','S3_health_overall_2010',...
% %     'p1_educ_overall_2010','p2_educ_overall_2010','p3_educ_overall_2010','p4_educ_overall_2010',...
% %     'p5_educ_overall_2010','p6_educ_overall_2010','p7_educ_overall_2010','p8_educ_overall_2010',...
% %     'p9_educ_overall_2010','p10_educ_overall_2010','p11_educ_overall_2010','p12_educ_overall_2010',...
% %     'p13_educ_overall_2010',...
% %     'theta12_overall_2010','theta13_overall_2010','theta23_overall_2010',...
% %     'rho12_overall_2010','rho13_overall_2010','rho23_overall_2010');
% 
% % subplot(3,3,1);plot(W1_income_overall_2001(1:5000));title('W1 income');
% % subplot(3,3,2);plot(W2_income_overall_2001(1:5000));title('W2 income');
% % subplot(3,3,3);plot(W3_income_overall_2001(1:5000));title('W3 income');
% % subplot(3,3,4);plot(M1_income_overall_2001(1:5000));title('M1 income');
% % subplot(3,3,5);plot(M2_income_overall_2001(1:5000));title('M2 income');
% % subplot(3,3,6);plot(M3_income_overall_2001(1:5000));title('M3 income');
% % subplot(3,3,7);plot(V1_income_overall_2001(1:5000));title('V1 income');
% % subplot(3,3,8);plot(V2_income_overall_2001(1:5000));title('V2 income');
% % subplot(3,3,9);plot(V3_income_overall_2001(1:5000));title('V3 income');
% 
% % subplot(3,3,1);plot(W1_health_overall_2001(1:5000));title('W1 health');
% % subplot(3,3,2);plot(W2_health_overall_2001(1:5000));title('W2 health');
% % subplot(3,3,3);plot(W3_health_overall_2001(1:5000));title('W3 health');
% % subplot(3,3,4);plot(M1_health_overall_2001(1:5000));title('M1 health');
% % subplot(3,3,5);plot(M2_health_overall_2001(1:5000));title('M2 health');
% % subplot(3,3,6);plot(M3_health_overall_2001(1:5000));title('M3 health');
% % subplot(3,3,7);plot(S1_health_overall_2001(1:5000));title('S1 health');
% % subplot(3,3,8);plot(S2_health_overall_2001(1:5000));title('S2 health');
% % subplot(3,3,9);plot(S3_health_overall_2001(1:5000));title('S3 health');
% 
% % subplot(3,3,1);plot(p1_educ_overall_2001(1:5000));title('p1 educ');
% % subplot(3,3,2);plot(p2_educ_overall_2001(1:5000));title('p2 educ');
% % subplot(3,3,3);plot(p3_educ_overall_2001(1:5000));title('p3 educ');
% % subplot(3,3,4);plot(p4_educ_overall_2001(1:5000));title('p4 educ');
% % subplot(3,3,5);plot(p5_educ_overall_2001(1:5000));title('p5 educ');
% % subplot(3,3,6);plot(p6_educ_overall_2001(1:5000));title('p6 educ');
% % subplot(3,3,7);plot(p7_educ_overall_2001(1:5000));title('p7 educ');
% % subplot(3,3,8);plot(p8_educ_overall_2001(1:5000));title('p8 educ');
% % subplot(3,3,9);plot(p9_educ_overall_2001(1:5000));title('p9 educ');
% 
% % subplot(1,3,1);plot(rho12_overall_2010(1:5000));title('rho12');
% % subplot(1,3,2);plot(rho13_overall_2010(1:5000));title('rho13');
% % subplot(1,3,3);plot(rho23_overall_2010(1:5000));title('rho23');