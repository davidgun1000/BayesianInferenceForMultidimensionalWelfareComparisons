%generate simulated data
load('param2001.mat');
num_group_educ=5;
num_group_happiness=5;
num_component=2;
dim_cop=4;
prop_marginal_educ_true = mean(prop_educ);
cum_marginal_educ = cumsum(prop_marginal_educ_true);
cum_marginal_educ(1,end)=1;
prop_marginal_happiness_true = mean(prop_happiness);
cum_marginal_happiness = cumsum(prop_marginal_happiness_true);
cum_marginal_happiness(1,end)=1;
% theta_gauss=0.8*ones(dim_cop,dim_cop);
% for i=1:dim_cop
%     theta_gauss(i,i)=1;
% end
theta_gauss_true = mean(theta_gauss)';
chol_gauss=diag([ones(dim_cop,1)]);
number_gauss = 1;
for i=1:dim_cop-1
    for j=i:dim_cop-1
        chol_gauss(i,j+1) = theta_gauss_true(number_gauss,1); 
        number_gauss=number_gauss+1;
    end

end
sigma = chol_gauss'*chol_gauss;
diag_mat = diag([1./sqrt(diag(sigma))]);
covmat = diag_mat*sigma*diag_mat;    
%theta_gaussian = [1 theta12 theta13; 0 1 theta23; 0 0 1];
%sig = (theta_gaussian'*theta_gaussian);
%diag_mat = [(1/sqrt(sig(1,1))) 0 0 ; 0 (1/sqrt(sig(2,2))) 0 ; 0 0 (1/sqrt(sig(3,3)))];
%cov_mat = diag_mat*sig*diag_mat;



n=20000;
u_mix=rand(n,1);
y_cop = mvnrnd(zeros(1,dim_cop),covmat,n);
u_cop = normcdf(y_cop);

w_gamma_true = [0.3, 0.7];
m_gamma_true = [50, 500];
v_gamma_true = mean(v_gamma(:,1:2));

w_beta_true = [0.3, 0.7];
m_beta_true = [0.5, 0.8];
s_beta_true = mean(s_beta(:,1:2));

num_rep=100000;
u1=rand(num_rep,1);
obs_income=[];
for j=1:num_component
    z_income=(u1>sum(w_gamma_true(1:(j-1)))) & (u1<=sum(w_gamma_true(1:j)));
    n_income=sum(z_income);
    obs_income_temp=gamrnd(v_gamma_true(j),m_gamma_true(j)/v_gamma_true(j),n_income,1);
    obs_income=[obs_income;obs_income_temp];
end
obs_income=sort(obs_income,1);
cdf_income=sum(repmat(w_gamma_true(1:num_component), num_rep, 1) .* gamcdf(repmat(obs_income, 1, num_component), repmat(v_gamma_true(1:num_component), num_rep, 1),repmat(m_gamma_true(1:num_component), num_rep, 1)./repmat(v_gamma_true(1:num_component), num_rep, 1)), 2);
cdf_income(end,1)=1;
table_income=[obs_income,cdf_income];
for i=1:n
    indx_income=find(u_cop(i,1)<=table_income(:,2),1,'first');
    if indx_income == 1
       data(i,1) = table_income(indx_income(1,1),1); 
    else
       data(i,1) = unifrnd(table_income(indx_income(1,1)-1,1),table_income(indx_income(1,1),1)); 
    end
end

u1=rand(num_rep,1);
obs_health=[];
for j=1:num_component
    z_health=(u1>sum(w_beta_true(1:(j-1)))) & (u1<=sum(w_beta_true(1:j)));
    n_health = sum(z_health);
    obs_health_temp=betarnd(s_beta_true(j)*m_beta_true(j),s_beta_true(j)*(1-m_beta_true(j)),n_health,1);
    obs_health=[obs_health;obs_health_temp];
    
end
obs_health = sort(obs_health,1);
cdf_health = sum(repmat(w_beta_true(1:num_component), num_rep, 1) .* ... 
        betacdf(repmat(obs_health, 1, num_component), repmat(m_beta_true(1:num_component), num_rep, 1).*repmat(s_beta_true(1:num_component), num_rep, 1),...
        repmat(s_beta_true(1:num_component), num_rep, 1).*(1-repmat(m_beta_true(1:num_component), num_rep, 1))), 2);
cdf_health(end,1)=1;
table_health=[obs_health,cdf_health];
for i=1:n
    indx_health=find(u_cop(i,2)<=table_health(:,2),1,'first');
    if indx_health == 1
       data(i,2) = table_health(indx_health(1,1),1); 
    else
       data(i,2) = unifrnd(table_health(indx_health(1,1)-1,1),table_health(indx_health(1,1),1)); 
    end
end
id_data2 = data(:,2)<=0.001;
data(id_data2,2) = 0.001;

id_data2 = data(:,2)>=0.999;
data(id_data2,2) = 0.999;

for i=1:num_group_educ
    if i==1
       id = (u_cop(:,3)<=cum_marginal_educ(i));
       data(id,3) = i;
    else
       id = (cum_marginal_educ(i-1)<=u_cop(:,3)) & (u_cop(:,3)<=cum_marginal_educ(i));
       data(id,3) = i;
    end
end

for i=1:num_group_happiness
    if i==1
       id = (u_cop(:,4)<=cum_marginal_happiness(i));
       data(id,4) = i;
    else
       id = (cum_marginal_happiness(i-1)<=u_cop(:,4)) & (u_cop(:,4)<=cum_marginal_happiness(i));
       data(id,4) = i;
    end
end



% theta_gauss_true=theta_gauss;
% prop_marginal_educ_true=prop_marginal_educ;
% prop_marginal_happiness_true=prop_marginal_happiness;
% 
% w_gamma_true=w_gamma;
% m_gamma_true=m_gamma;
% v_gamma_true=v_gamma;
% 
% w_beta_true=w_beta;
% m_beta_true=m_beta;
% s_beta_true=s_beta;
y_cop_true=y_cop;
u_cop_true=u_cop;
%theta_gauss = theta_gauss_true;

save('simulated_data.mat','theta_gauss_true','data','prop_marginal_educ_true','prop_marginal_happiness_true',...
     'w_gamma_true','m_gamma_true','v_gamma_true','w_beta_true','m_beta_true','s_beta_true','y_cop_true','u_cop_true','theta_gauss_true');

