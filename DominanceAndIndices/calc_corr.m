%calculate pairwise correlation between income, mental health, education,
%and happiness. 
rng(123456)
load('param2019.mat');
dim_cop=4;
length_draws = length(m_beta(:,1));
for m=1:length_draws

chol_gauss=diag([ones(dim_cop,1)]);
number_gauss = 1;
for i=1:dim_cop-1
    for j=i:dim_cop-1
        chol_gauss(i,j+1) = theta_gauss(m,number_gauss); 
        number_gauss=number_gauss+1;
    end

end
sigma = chol_gauss'*chol_gauss;
diag_mat = diag([1./sqrt(diag(sigma))]);
covmat = diag_mat*sigma*diag_mat;          

corr_income_health(m,1) = covmat(2,1);
corr_income_educ(m,1) = covmat(3,1);
corr_income_happiness(m,1) = covmat(4,1);
corr_health_educ(m,1) = covmat(2,3);
corr_health_happiness(m,1) = covmat(2,4);
corr_educ_happiness(m,1) = covmat(3,4);

Kendall_income_health(m,1) = copulastat('Gaussian',corr_income_health(m,1));
Kendall_income_educ(m,1) = copulastat('Gaussian',corr_income_educ(m,1));
Kendall_income_happiness(m,1) = copulastat('Gaussian',corr_income_happiness(m,1));
Kendall_health_educ(m,1) = copulastat('Gaussian',corr_health_educ(m,1));
Kendall_health_happiness(m,1) = copulastat('Gaussian',corr_health_happiness(m,1));
Kendall_educ_happiness(m,1) = copulastat('Gaussian',corr_educ_happiness(m,1));


end

save('corr_2019.mat','corr_income_health','corr_income_educ','corr_income_happiness',...
                     'corr_health_educ','corr_health_happiness','corr_educ_happiness',...
                     'Kendall_income_health','Kendall_income_educ','Kendall_income_happiness',...
                     'Kendall_health_educ','Kendall_health_happiness','Kendall_educ_happiness');
                 
mean(Kendall_income_health)
std(Kendall_income_health)
mean(Kendall_income_educ)
std(Kendall_income_educ)
mean(Kendall_income_happiness)
std(Kendall_income_happiness)
mean(Kendall_health_educ)
std(Kendall_health_educ)
mean(Kendall_health_happiness)
std(Kendall_health_happiness)
mean(Kendall_educ_happiness)
std(Kendall_educ_happiness)
