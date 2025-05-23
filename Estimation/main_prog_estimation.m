%mainprogfourdimensionCopulaEstimation
num_FPBB=200;% number of pseudo representative samples 
year_x='2015'; 
load(['wellbeing_data_',num2str(year_x),'.mat']); %load the data

data_x = data_overall_2015; 
samplingweight_x = weight_overall_2015;
n_x = length(data_x);


num_component_gamma_x = 2; %number of components in the gamma mixture model
num_component_beta_x = 3; %number of components in the beta mixture model


obs_FPBB_x = generate_pseudo_representative(data_x,samplingweight_x,n_x); % generate pseudo representative sample
[prior_x,init_param_x,V1_gamma_com_x,V1_beta_com_x,V1_delta_educ_x,V1_delta_happiness_x,scale_theta_copula_x]=InitialSamplerCopula(obs_FPBB_x,num_component_gamma_x,num_component_beta_x); %to generate initial values for the main MCMC
parpool(48)
parfor j=1:num_FPBB
    obs_FPBB_x = generate_pseudo_representative(data_x,samplingweight_x,n_x); %generate pseudo representative sample
    %sample the copula model parameters using MCMC
    [w_gamma_temp_x(:,:,j),m_gamma_temp_x(:,:,j),v_gamma_temp_x(:,:,j),...
        w_beta_temp_x(:,:,j),m_beta_temp_x(:,:,j),s_beta_temp_x(:,:,j),delta_educ_temp_x(:,:,j),prop_educ_temp_x(:,:,j),delta_happiness_temp_x(:,:,j),prop_happiness_temp_x(:,:,j),theta_gauss_temp_x(:,:,j)]=MainSamplerCopula(obs_FPBB_x,prior_x,init_param_x,...
                                                             V1_gamma_com_x,V1_beta_com_x,V1_delta_educ_x,V1_delta_happiness_x,scale_theta_copula_x,num_component_gamma_x,num_component_beta_x);
                                                         

end

w_gamma_x=[];
m_gamma_x=[];
v_gamma_x=[];
w_beta_x=[];
m_beta_x=[];
s_beta_x=[];
theta_gauss_x=[];
prop_educ_x=[];
prop_happiness_x=[];

for j=1:num_FPBB
    w_gamma_x=[w_gamma_x;w_gamma_temp_x(:,:,j)];
    m_gamma_x=[m_gamma_x;m_gamma_temp_x(:,:,j)];
    v_gamma_x=[v_gamma_x;v_gamma_temp_x(:,:,j)];
    w_beta_x=[w_beta_x;w_beta_temp_x(:,:,j)];
    m_beta_x=[m_beta_x;m_beta_temp_x(:,:,j)];
    s_beta_x=[s_beta_x;s_beta_temp_x(:,:,j)];
    theta_gauss_x=[theta_gauss_x;theta_gauss_temp_x(:,:,j)];
    prop_educ_x=[prop_educ_x;prop_educ_temp_x(:,:,j)];
    prop_happiness_x=[prop_happiness_x;prop_happiness_temp_x(:,:,j)];

    
    
    
end
%save the results
save(['param',num2str(year_x),'.mat'],'w_gamma_x','m_gamma_x','v_gamma_x','w_beta_x','m_beta_x','s_beta_x','theta_gauss_x','prop_educ_x','prop_happiness_x');



                                                                                
                                        
                                        