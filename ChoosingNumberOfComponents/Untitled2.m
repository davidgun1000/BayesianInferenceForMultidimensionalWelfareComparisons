parpool(28)
N=600000;
num_FPBB=200;
year_x='2019';
load(['wellbeing_data_',num2str(year_x),'.mat']);
data_train = data_overall_2019(1:round(0.5*length(data_overall_2019(:,1))),2);
data_test = data_overall_2019(round(0.5*length(data_overall_2019(:,1)))+1:length(data_overall_2019(:,1)),2);

weight_train = weight_overall_2019(1:round(0.5*length(data_overall_2019(:,1))),1);
weight_test = weight_overall_2019(round(0.5*length(data_overall_2019(:,1)))+1:length(data_overall_2019(:,1)),1);

data_x = data_train;
n_x = length(data_x);
samplingweight_x = weight_train./sum(weight_train).*N;

num_component_beta=1;

[prior_x,init_param_x,V1_beta_com_x]=InitialSamplerBetaMix(data_x,num_component_beta);

parfor j=1:num_FPBB
    obs_FPBB_x = generate_pseudo_representative(data_x,samplingweight_x,n_x);
    [w_beta_temp_x(:,:,j),m_beta_temp_x(:,:,j),s_beta_temp_x(:,:,j)]=MainSamplerBetaMix(obs_FPBB_x,prior_x,init_param_x,...
                                                             V1_beta_com_x,num_component_beta);

end

w_beta_x=[];
m_beta_x=[];
s_beta_x=[];

for j=1:num_FPBB
    w_beta_x=[w_beta_x;w_beta_temp_x(:,:,j)];
    m_beta_x=[m_beta_x;m_beta_temp_x(:,:,j)];
    s_beta_x=[s_beta_x;s_beta_temp_x(:,:,j)];    
end

normalised_weight_test = weight_test./sum(weight_test);
com = [data_test,normalised_weight_test];
com_sort = sortrows(com,1);
com_sort(:,2) = cumsum(com_sort(:,2));

num_draws = length(w_beta_x);
for k=1:length(data_test)
    cdf_beta(:,k) = zeros(num_draws,1);
    pdf_beta(:,k) = zeros(num_draws,1);
for j=1:num_component_beta
    cdf_beta(:,k) = cdf_beta(:,k) + w_beta_x(:,j).*betacdf(com_sort(k,1),s_beta_x(:,j).*m_beta_x(:,j),s_beta_x(:,j).*(1-m_beta_x(:,j)));
    pdf_beta(:,k) = pdf_beta(:,k) + w_beta_x(:,j).*betaden(com_sort(k,1),s_beta_x(:,j),m_beta_x(:,j));
end
end

mean_cdf_beta = (mean(cdf_beta))';
mean_pdf_beta = (mean(pdf_beta))';

MAE = mean(abs(mean_cdf_beta - com_sort(:,2)));
MSE = sqrt(mean((mean_cdf_beta - com_sort(:,2)).^2));
Log_score = sum(log(mean_pdf_beta));
save(['train_test_health','_',num2str(year_x),'_',num2str(num_component_beta),'.mat'],'MAE','MSE','Log_score','w_beta_x','m_beta_x','s_beta_x');