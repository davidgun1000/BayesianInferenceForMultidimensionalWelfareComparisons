%parpool(28)

%model selection
%user needs to supply the data
N=600000;
num_FPBB=200;
year_x='2001';
load(['wellbeing_data_train_test_',num2str(year_x),'.mat']);



data_x = data_train;
n_x = length(data_x);
samplingweight_x = weight_train./sum(weight_train).*N;

num_component_gamma=3;
num_component_beta=2;

[prior_income,init_param_income,V1_gamma_com_income]=InitialSamplerGammaMix(data_x(:,1),num_component_gamma);
[prior_health,init_param_health,V1_beta_com_health]=InitialSamplerBetaMix(data_x(:,2),num_component_beta);
parfor j=1:num_FPBB
    obs_FPBB_x = generate_pseudo_representative(data_x,samplingweight_x,n_x);
    
    obs_income = obs_FPBB_x(:,1);
    obs_health = obs_FPBB_x(:,2);
    
    [w_gamma_temp_x(:,:,j),m_gamma_temp_x(:,:,j),v_gamma_temp_x(:,:,j)]=MainSamplerGammaMix(obs_income,prior_income,init_param_income,...
                                                             V1_gamma_com_income,num_component_gamma);

    [w_beta_temp_x(:,:,j),m_beta_temp_x(:,:,j),s_beta_temp_x(:,:,j)]=MainSamplerBetaMix(obs_health,prior_health,init_param_health,...
                                                             V1_beta_com_health,num_component_beta);                                                     
end

w_gamma_x=[];
m_gamma_x=[];
v_gamma_x=[];

w_beta_x=[];
m_beta_x=[];
s_beta_x=[];
for j=1:num_FPBB
    w_gamma_x=[w_gamma_x;w_gamma_temp_x(:,:,j)];
    m_gamma_x=[m_gamma_x;m_gamma_temp_x(:,:,j)];
    v_gamma_x=[v_gamma_x;v_gamma_temp_x(:,:,j)];   
    
    w_beta_x=[w_beta_x;w_beta_temp_x(:,:,j)];
    m_beta_x=[m_beta_x;m_beta_temp_x(:,:,j)];
    s_beta_x=[s_beta_x;s_beta_temp_x(:,:,j)]; 
end

normalised_weight_test = weight_test./sum(weight_test);

com_income = [data_test(:,1),normalised_weight_test];
com_sort_income = sortrows(com_income,1);
com_sort_income(:,2) = cumsum(com_sort_income(:,2));
num_draws = length(w_gamma_x);
for k=1:length(data_test)
    cdf_gamma(:,k) = zeros(num_draws,1);
    pdf_gamma(:,k) = zeros(num_draws,1);
for j=1:num_component_gamma
    cdf_gamma(:,k) = cdf_gamma(:,k) + w_gamma_x(:,j).*gamcdf(com_sort_income(k,1),v_gamma_x(:,j),m_gamma_x(:,j)./v_gamma_x(:,j));
    pdf_gamma(:,k) = pdf_gamma(:,k) + w_gamma_x(:,j).*gampdf(com_sort_income(k,1),v_gamma_x(:,j),m_gamma_x(:,j)./v_gamma_x(:,j));
end
end
mean_cdf_gamma = (mean(cdf_gamma))';
mean_pdf_gamma = (mean(pdf_gamma))';
MAE_income = mean(abs(mean_cdf_gamma - com_sort_income(:,2)));
MSE_income = sqrt(mean((mean_cdf_gamma - com_sort_income(:,2)).^2));
Log_score_income = sum(log(mean_pdf_gamma));

com_health = [data_test(:,2),normalised_weight_test];
com_sort_health = sortrows(com_health,1);
com_sort_health(:,2) = cumsum(com_sort_health(:,2));
for k=1:length(data_test)
    cdf_beta(:,k) = zeros(num_draws,1);
    pdf_beta(:,k) = zeros(num_draws,1);
for j=1:num_component_beta
    cdf_beta(:,k) = cdf_beta(:,k) + w_beta_x(:,j).*betacdf(com_sort_health(k,1),s_beta_x(:,j).*m_beta_x(:,j),s_beta_x(:,j).*(1-m_beta_x(:,j)));
    pdf_beta(:,k) = pdf_beta(:,k) + w_beta_x(:,j).*betaden(com_sort_health(k,1),s_beta_x(:,j),m_beta_x(:,j));
end
end
mean_cdf_beta = (mean(cdf_beta))';
mean_pdf_beta = (mean(pdf_beta))';

MAE_health = mean(abs(mean_cdf_beta - com_sort_health(:,2)));
MSE_health = sqrt(mean((mean_cdf_beta - com_sort_health(:,2)).^2));
Log_score_health = sum(log(mean_pdf_beta));
save(['train_test_income_health','_',num2str(year_x),'_',num2str(num_component_gamma),'.mat'],'MAE_income','MSE_income','Log_score_income','w_gamma_x','m_gamma_x','v_gamma_x',...
                                                                                              'MAE_health','MSE_health','Log_score_health','w_beta_x','m_beta_x','s_beta_x');