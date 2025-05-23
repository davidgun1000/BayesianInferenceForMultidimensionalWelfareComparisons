%this function calculates multidimensional indices and dominance
%probabilities
rng(123456);

num_FPBB=200;
year_x='2001'; %load the data
year_y='2010';
%load(['wellbeing_data_',num2str(year_x),'.mat']);
%load(['wellbeing_data_',num2str(year_y),'.mat']);

% data_x = data_overall_2001;
% samplingweight_x = weight_overall_2001;
% n_x = length(data_x);
% 
% data_y = data_overall_2010;
% samplingweight_y = weight_overall_2010;
% n_y = length(data_y);

num_component_gamma_x = 3;
num_component_beta_x = 3;

num_component_gamma_y = 2;
num_component_beta_y = 3;


%load model parameters
load(['param',num2str(year_x),'.mat']);
w_gamma_x = w_gamma;
m_gamma_x = m_gamma;
v_gamma_x = v_gamma;
w_beta_x = w_beta;
m_beta_x = m_beta;
s_beta_x = s_beta;
prop_educ_x = prop_educ;
prop_happiness_x = prop_happiness;
theta_gauss_x = theta_gauss;


load(['param',num2str(year_y),'.mat']);
w_gamma_y = w_gamma;
m_gamma_y = m_gamma;
v_gamma_y = v_gamma;
w_beta_y = w_beta;
m_beta_y = m_beta;
s_beta_y = s_beta;
prop_educ_y = prop_educ;
prop_happiness_y = prop_happiness;
theta_gauss_y = theta_gauss;


%grids
grid_income = (4:0.0335:7.3132)';
grid_health = 0.01:0.01:0.99;
grid_income = (exp(grid_income));
grid_health = grid_health';
grid_educ = [1;2;3;4;5];
grid_happiness = [1;2;3;4;5];

%poverty lines
poverty_line.income=40;
poverty_line.health=50;
poverty_line.educ=2;
poverty_line.happiness=2;

num_param = length(w_gamma_x);

F_income_health_x = [];
F_income_educ_x = [];
F_income_happiness_x = []; 
H_income_health_x = [];
H_income_educ_x = [];
H_income_happiness_x = [];

F_income_health_y = [];
F_income_educ_y = [];
F_income_happiness_y = []; 
H_income_health_y = [];
H_income_educ_y = [];
H_income_happiness_y = [];
%parpool(48)
parfor j=1:num_param
    %calculate required quantities for dominance analysis and compute
    %indices.
    [output_x]=obtain_dominance_matrix_v3(w_gamma_x(j,:),m_gamma_x(j,:),v_gamma_x(j,:),w_beta_x(j,:),m_beta_x(j,:),s_beta_x(j,:),theta_gauss_x(j,:),prop_educ_x(j,:),prop_happiness_x(j,:),poverty_line,grid_income,grid_health,grid_educ,grid_happiness,...
                num_component_gamma_x,num_component_beta_x);
    mean_income_x(j,1) = output_x.mean_income;
    mean_health_x(j,1) = output_x.mean_health;
    headcount_income_x(j,1) = output_x.headcount_income;
    headcount_health_x(j,1) = output_x.headcount_health;
    headcount_educ_x(j,1) = output_x.headcount_educ;
    headcount_happiness_x(j,1) = output_x.headcount_happiness;
    multi_headcount_k1_x(j,1) = output_x.multi_headcount_k1;
    multi_headcount_k2_x(j,1) = output_x.multi_headcount_k2;
    multi_headcount_k3_x(j,1) = output_x.multi_headcount_k3;
    multi_headcount_k4_x(j,1) = output_x.multi_headcount_k4;
    multi_headcount_x(j,:)= [multi_headcount_k1_x(j,1),multi_headcount_k2_x(j,1),multi_headcount_k3_x(j,1),multi_headcount_k4_x(j,1)];
    multi_M0_k1_x(j,1) = output_x.multi_M0_k1;
    multi_M0_k2_x(j,1) = output_x.multi_M0_k2;
    multi_M0_k3_x(j,1) = output_x.multi_M0_k3;
    multi_M0_k4_x(j,1) = output_x.multi_M0_k4;
    multi_M0_x(j,:) = [multi_M0_k1_x(j,1),multi_M0_k2_x(j,1),multi_M0_k3_x(j,1),multi_M0_k4_x(j,1)];
    FGT1_income_x(j,1) = output_x.FGT1_income;
    FGT2_income_x(j,1) = output_x.FGT2_income;
    FGT1_health_x(j,1) = output_x.FGT1_health;
    FGT2_health_x(j,1) = output_x.FGT2_health;
    CF01_education_x(j,1) = output_x.CF01_education;
    CF09_education_x(j,1) = output_x.CF09_education;
    CF01_happiness_x(j,1) = output_x.CF01_happiness;
    CF09_happiness_x(j,1) = output_x.CF09_happiness;    
    Gini_income_x(j,1) = output_x.Gini_income;
    Gini_health_x(j,1) = output_x.Gini_health;
    F_income_x(j,:) = output_x.F_income;
    H_income_x(j,:) = output_x.H_income;
    F_health_x(j,:) = output_x.F_health;
    H_health_x(j,:) = output_x.H_health;
    F_educ_x(j,:) = output_x.F_educ;
    H_educ_x(j,:) = output_x.H_educ;
    F_happiness_x(j,:) = output_x.F_happiness;
    H_happiness_x(j,:) = output_x.H_happiness;
    F_income_health_x = [F_income_health_x;output_x.F_income_health];
    F_income_educ_x = [F_income_educ_x;output_x.F_income_educ];
    F_income_happiness_x = [F_income_happiness_x;output_x.F_income_happiness];
    H_income_health_x = [H_income_health_x;output_x.H_income_health];
    H_income_educ_x = [H_income_educ_x;output_x.H_income_educ];
    H_income_happiness_x = [H_income_happiness_x;output_x.H_income_happiness];
    
    [output_y]=obtain_dominance_matrix_v3(w_gamma_y(j,:),m_gamma_y(j,:),v_gamma_y(j,:),w_beta_y(j,:),m_beta_y(j,:),s_beta_y(j,:),theta_gauss_y(j,:),prop_educ_y(j,:),prop_happiness_y(j,:),poverty_line,grid_income,grid_health,grid_educ,grid_happiness,...
                num_component_gamma_y,num_component_beta_y);
    mean_income_y(j,1) = output_y.mean_income;
    mean_health_y(j,1) = output_y.mean_health;
    headcount_income_y(j,1) = output_y.headcount_income;
    headcount_health_y(j,1) = output_y.headcount_health;
    headcount_educ_y(j,1) = output_y.headcount_educ;
    headcount_happiness_y(j,1) = output_y.headcount_happiness;
    multi_headcount_k1_y(j,1) = output_y.multi_headcount_k1;
    multi_headcount_k2_y(j,1) = output_y.multi_headcount_k2;
    multi_headcount_k3_y(j,1) = output_y.multi_headcount_k3;
    multi_headcount_k4_y(j,1) = output_y.multi_headcount_k4;
    multi_headcount_y(j,:)= [multi_headcount_k1_y(j,1),multi_headcount_k2_y(j,1),multi_headcount_k3_y(j,1),multi_headcount_k4_y(j,1)];
    multi_M0_k1_y(j,1) = output_y.multi_M0_k1;
    multi_M0_k2_y(j,1) = output_y.multi_M0_k2;
    multi_M0_k3_y(j,1) = output_y.multi_M0_k3;
    multi_M0_k4_y(j,1) = output_y.multi_M0_k4;
    multi_M0_y(j,:) = [multi_M0_k1_y(j,1),multi_M0_k2_y(j,1),multi_M0_k3_y(j,1),multi_M0_k4_y(j,1)];
    FGT1_income_y(j,1) = output_y.FGT1_income;
    FGT2_income_y(j,1) = output_y.FGT2_income;
    FGT1_health_y(j,1) = output_y.FGT1_health;
    FGT2_health_y(j,1) = output_y.FGT2_health;
    CF01_education_y(j,1) = output_y.CF01_education;
    CF09_education_y(j,1) = output_y.CF09_education;
    CF01_happiness_y(j,1) = output_y.CF01_happiness;
    CF09_happiness_y(j,1) = output_y.CF09_happiness;
    Gini_income_y(j,1) = output_y.Gini_income;
    Gini_health_y(j,1) = output_y.Gini_health;
    F_income_y(j,:) = output_y.F_income;
    H_income_y(j,:) = output_y.H_income;
    F_health_y(j,:) = output_y.F_health;
    H_health_y(j,:) = output_y.H_health;
    F_educ_y(j,:) = output_y.F_educ;
    H_educ_y(j,:) = output_y.H_educ;
    F_happiness_y(j,:) = output_y.F_happiness;
    H_happiness_y(j,:) = output_y.H_happiness;
    F_income_health_y = [F_income_health_y;output_y.F_income_health];
    F_income_educ_y = [F_income_educ_y;output_y.F_income_educ];
    F_income_happiness_y = [F_income_happiness_y;output_y.F_income_happiness];
    H_income_health_y = [H_income_health_y;output_y.H_income_health];
    H_income_educ_y = [H_income_educ_y;output_y.H_income_educ];
    H_income_happiness_y = [H_income_happiness_y;output_y.H_income_happiness];
    
end

m=size(mean_income_x,1);

for j=1:1000
    
    R = (randperm(m))';
    for i=1:m
        F_income_y_perm(i,:) = F_income_y(R(i,1),:);
        H_income_y_perm(i,:) = H_income_y(R(i,1),:);
        F_health_y_perm(i,:) = F_health_y(R(i,1),:);
        H_health_y_perm(i,:) = H_health_y(R(i,1),:);
        F_educ_y_perm(i,:) = F_educ_y(R(i,1),:);
        H_educ_y_perm(i,:) = H_educ_y(R(i,1),:);
        F_happiness_y_perm(i,:) = F_happiness_y(R(i,1),:);
        H_happiness_y_perm(i,:) = H_happiness_y(R(i,1),:);
        F_income_health_y_perm(i,:,:) = F_income_health_y(R(i,1),:,:);
        F_income_educ_y_perm(i,:,:) = F_income_educ_y(R(i,1),:,:);
        F_income_happiness_y_perm(i,:,:) = F_income_happiness_y(R(i,1),:,:);
        H_income_health_y_perm(i,:,:) = H_income_health_y(R(i,1),:,:);
        H_income_educ_y_perm(i,:,:) = H_income_educ_y(R(i,1),:,:);
        H_income_happiness_y_perm(i,:,:) = H_income_happiness_y(R(i,1),:,:);
        multi_headcount_y_perm(i,:) = multi_headcount_y(R(i,1),:);
        multi_M0_y_perm(i,:) = multi_M0_y(R(i,1),:);
    end
    


    %First order dominance for income
    xFSDy_income = F_income_x<=F_income_y_perm;
    yFSDx_income = F_income_y_perm<=F_income_x;
    xFSDy_income = double(xFSDy_income);
    yFSDx_income = double(yFSDx_income);
    prop_xFSDy_income(j,:) = mean(xFSDy_income);
    prop_yFSDx_income(j,:) = mean(yFSDx_income);
    overall_xFSDy_income(j,1) = mean(prod(xFSDy_income,2));
    overall_yFSDx_income(j,1) = mean(prod(yFSDx_income,2));
    overall_xFSDy_income_lowest(j,1) = mean(prod(xFSDy_income(:,1:poverty_line.income),2));
    overall_yFSDx_income_lowest(j,1) = mean(prod(yFSDx_income(:,1:poverty_line.income),2));

    %second order dominance for income
    xSSDy_income = H_income_x<=H_income_y_perm;
    ySSDx_income = H_income_y_perm<=H_income_x;
    xSSDy_income = double(xSSDy_income);
    ySSDx_income = double(ySSDx_income);
    prop_xSSDy_income(j,:) = mean(xSSDy_income);
    prop_ySSDx_income(j,:) = mean(ySSDx_income);
    overall_xSSDy_income(j,1) = mean(prod(xSSDy_income,2));
    overall_ySSDx_income(j,1) = mean(prod(ySSDx_income,2));
    overall_xSSDy_income_lowest(j,1) = mean(prod(xSSDy_income(:,1:poverty_line.income),2));
    overall_ySSDx_income_lowest(j,1) = mean(prod(ySSDx_income(:,1:poverty_line.income),2));
    
    %first order dominance for health
    xFSDy_health = F_health_x<=F_health_y_perm;
    yFSDx_health = F_health_y_perm<=F_health_x;
    xFSDy_health = double(xFSDy_health);
    yFSDx_health = double(yFSDx_health);
    prop_xFSDy_health(j,:) = mean(xFSDy_health);
    prop_yFSDx_health(j,:) = mean(yFSDx_health);
    overall_xFSDy_health(j,1) = mean(prod(xFSDy_health,2));
    overall_yFSDx_health(j,1) = mean(prod(yFSDx_health,2));
    overall_xFSDy_health_lowest(j,1) = mean(prod(xFSDy_health(:,1:poverty_line.health),2));
    overall_yFSDx_health_lowest(j,1) = mean(prod(yFSDx_health(:,1:poverty_line.health),2));

    %second order dominance for health
    xSSDy_health = H_health_x<=H_health_y_perm;
    ySSDx_health = H_health_y_perm<=H_health_x;
    xSSDy_health = double(xSSDy_health);
    ySSDx_health = double(ySSDx_health);
    prop_xSSDy_health(j,:) = mean(xSSDy_health);
    prop_ySSDx_health(j,:) = mean(ySSDx_health);
    overall_xSSDy_health(j,1) = mean(prod(xSSDy_health,2));
    overall_ySSDx_health(j,1) = mean(prod(ySSDx_health,2));
    overall_xSSDy_health_lowest(j,1) = mean(prod(xSSDy_health(:,1:poverty_line.health),2));
    overall_ySSDx_health_lowest(j,1) = mean(prod(ySSDx_health(:,1:poverty_line.health),2));
   
    
    %first order dominance for education
    xFSDy_educ = F_educ_x<=F_educ_y_perm;
    yFSDx_educ = F_educ_y_perm<=F_educ_x;
    xFSDy_educ = double(xFSDy_educ);
    yFSDx_educ = double(yFSDx_educ);
    prop_xFSDy_educ(j,:) = mean(xFSDy_educ);
    prop_yFSDx_educ(j,:) = mean(yFSDx_educ);
    overall_xFSDy_educ(j,1) = mean(prod(xFSDy_educ(:,1:end-1),2));
    overall_yFSDx_educ(j,1) = mean(prod(yFSDx_educ(:,1:end-1),2));
    overall_xFSDy_educ_lowest(j,1) = mean(prod(xFSDy_educ(:,1:poverty_line.educ),2));
    overall_yFSDx_educ_lowest(j,1) = mean(prod(yFSDx_educ(:,1:poverty_line.educ),2));
    
    %second order dominance for education
    xSSDy_educ = H_educ_x<=H_educ_y_perm;
    ySSDx_educ = H_educ_y_perm<=H_educ_x;
    xSSDy_educ = double(xSSDy_educ);
    ySSDx_educ = double(ySSDx_educ);
    prop_xSSDy_educ(j,:) = mean(xSSDy_educ);
    prop_ySSDx_educ(j,:) = mean(ySSDx_educ);
    overall_xSSDy_educ(j,1) = mean(prod(xSSDy_educ(:,1:end-1),2));
    overall_ySSDx_educ(j,1) = mean(prod(ySSDx_educ(:,1:end-1),2));
    overall_xSSDy_educ_lowest(j,1) = mean(prod(xSSDy_educ(:,1:poverty_line.educ),2));
    overall_ySSDx_educ_lowest(j,1) = mean(prod(ySSDx_educ(:,1:poverty_line.educ),2));

    %first order dominance for happiness
    xFSDy_happiness = F_happiness_x<=F_happiness_y_perm;
    yFSDx_happiness = F_happiness_y_perm<=F_happiness_x;
    xFSDy_happiness = double(xFSDy_happiness);
    yFSDx_happiness = double(yFSDx_happiness);
    prop_xFSDy_happiness(j,:) = mean(xFSDy_happiness);
    prop_yFSDx_happiness(j,:) = mean(yFSDx_happiness);
    overall_xFSDy_happiness(j,1) = mean(prod(xFSDy_happiness(:,1:end-1),2));
    overall_yFSDx_happiness(j,1) = mean(prod(yFSDx_happiness(:,1:end-1),2));
    overall_xFSDy_happiness_lowest(j,1) = mean(prod(xFSDy_happiness(:,1:poverty_line.happiness),2));
    overall_yFSDx_happiness_lowest(j,1) = mean(prod(yFSDx_happiness(:,1:poverty_line.happiness),2));

    %second order dominance for happiness
    xSSDy_happiness = H_happiness_x<=H_happiness_y_perm;
    ySSDx_happiness = H_happiness_y_perm<=H_happiness_x;
    xSSDy_happiness = double(xSSDy_happiness);
    ySSDx_happiness = double(ySSDx_happiness);
    prop_xSSDy_happiness(j,:) = mean(xSSDy_happiness);
    prop_ySSDx_happiness(j,:) = mean(ySSDx_happiness);
    overall_xSSDy_happiness(j,1) = mean(prod(xSSDy_happiness(:,1:end-1),2));
    overall_ySSDx_happiness(j,1) = mean(prod(ySSDx_happiness(:,1:end-1),2));
    overall_xSSDy_happiness_lowest(j,1) = mean(prod(xSSDy_happiness(:,1:poverty_line.happiness),2));
    overall_ySSDx_happiness_lowest(j,1) = mean(prod(ySSDx_happiness(:,1:poverty_line.happiness),2));
    
    %bivariate income and health
    xFSDy_income_health = F_income_health_x<=F_income_health_y_perm;
    yFSDx_income_health = F_income_health_y_perm<=F_income_health_x;
    xFSDy_income_health = double(xFSDy_income_health);
    yFSDx_income_health = double(yFSDx_income_health);
    xFSDy_income_health_intersection = xFSDy_income_health(:,1:poverty_line.health,1:poverty_line.income);
    yFSDx_income_health_intersection = yFSDx_income_health(:,1:poverty_line.health,1:poverty_line.income);
    prop_xFSDy_income_health(j,:,:) = mean(xFSDy_income_health);
    prop_yFSDx_income_health(j,:,:) = mean(yFSDx_income_health);
    overall_xFSDy_income_health(j,1) = mean(prod(prod(xFSDy_income_health,2),3));
    overall_yFSDx_income_health(j,1) = mean(prod(prod(yFSDx_income_health,2),3));
    overall_xFSDy_income_health_intersection(j,1) = mean(prod(prod(xFSDy_income_health_intersection,2),3));
    overall_yFSDx_income_health_intersection(j,1) = mean(prod(prod(yFSDx_income_health_intersection,2),3));
    
    %bivariate income and health
    xSSDy_income_health = H_income_health_x<=H_income_health_y_perm;
    ySSDx_income_health = H_income_health_y_perm<=H_income_health_x;
    xSSDy_income_health = double(xSSDy_income_health);
    ySSDx_income_health = double(ySSDx_income_health);
    xSSDy_income_health_intersection = xSSDy_income_health(:,1:poverty_line.health,1:poverty_line.income);
    ySSDx_income_health_intersection = ySSDx_income_health(:,1:poverty_line.health,1:poverty_line.income);
    prop_xSSDy_income_health(j,:,:) = mean(xSSDy_income_health);
    prop_ySSDx_income_health(j,:,:) = mean(ySSDx_income_health);
    overall_xSSDy_income_health(j,1) = mean(prod(prod(xSSDy_income_health,2),3));
    overall_ySSDx_income_health(j,1) = mean(prod(prod(ySSDx_income_health,2),3));
    overall_xSSDy_income_health_intersection(j,1) = mean(prod(prod(xSSDy_income_health_intersection,2),3));
    overall_ySSDx_income_health_intersection(j,1) = mean(prod(prod(ySSDx_income_health_intersection,2),3));
    
    
    %bivariate income and education
    xFSDy_income_educ = F_income_educ_x<=F_income_educ_y_perm;
    yFSDx_income_educ = F_income_educ_y_perm<=F_income_educ_x;
    xFSDy_income_educ = double(xFSDy_income_educ);
    yFSDx_income_educ = double(yFSDx_income_educ);
    xFSDy_income_educ_intersection = xFSDy_income_educ(:,1:poverty_line.educ,1:poverty_line.income);
    yFSDx_income_educ_intersection = yFSDx_income_educ(:,1:poverty_line.educ,1:poverty_line.income);
    prop_xFSDy_income_educ(j,:,:) = mean(xFSDy_income_educ);
    prop_yFSDx_income_educ(j,:,:) = mean(yFSDx_income_educ);
    overall_xFSDy_income_educ(j,1) = mean(prod(prod(xFSDy_income_educ,2),3));
    overall_yFSDx_income_educ(j,1) = mean(prod(prod(yFSDx_income_educ,2),3));
    overall_xFSDy_income_educ_intersection(j,1) = mean(prod(prod(xFSDy_income_educ_intersection,2),3));
    overall_yFSDx_income_educ_intersection(j,1) = mean(prod(prod(yFSDx_income_educ_intersection,2),3));
    
    %bivariate income and education
    xSSDy_income_educ = H_income_educ_x<=H_income_educ_y_perm;
    ySSDx_income_educ = H_income_educ_y_perm<=H_income_educ_x;
    xSSDy_income_educ = double(xSSDy_income_educ);
    ySSDx_income_educ = double(ySSDx_income_educ);
    xSSDy_income_educ_intersection = xSSDy_income_educ(:,1:poverty_line.educ,1:poverty_line.income);
    ySSDx_income_educ_intersection = ySSDx_income_educ(:,1:poverty_line.educ,1:poverty_line.income);
    prop_xSSDy_income_educ(j,:,:) = mean(xSSDy_income_educ);
    prop_ySSDx_income_educ(j,:,:) = mean(ySSDx_income_educ);
    overall_xSSDy_income_educ(j,1) = mean(prod(prod(xSSDy_income_educ,2),3));
    overall_ySSDx_income_educ(j,1) = mean(prod(prod(ySSDx_income_educ,2),3));
    overall_xSSDy_income_educ_intersection(j,1) = mean(prod(prod(xSSDy_income_educ_intersection,2),3));
    overall_ySSDx_income_educ_intersection(j,1) = mean(prod(prod(ySSDx_income_educ_intersection,2),3));

    %bivariate income and happiness
    xFSDy_income_happiness = F_income_happiness_x<=F_income_happiness_y_perm;
    yFSDx_income_happiness = F_income_happiness_y_perm<=F_income_happiness_x;
    xFSDy_income_happiness = double(xFSDy_income_happiness);
    yFSDx_income_happiness = double(yFSDx_income_happiness);
    xFSDy_income_happiness_intersection = xFSDy_income_happiness(:,1:poverty_line.happiness,1:poverty_line.income);
    yFSDx_income_happiness_intersection = yFSDx_income_happiness(:,1:poverty_line.happiness,1:poverty_line.income);
    prop_xFSDy_income_happiness(j,:,:) = mean(xFSDy_income_happiness);
    prop_yFSDx_income_happiness(j,:,:) = mean(yFSDx_income_happiness);
    overall_xFSDy_income_happiness(j,1) = mean(prod(prod(xFSDy_income_happiness,2),3));
    overall_yFSDx_income_happiness(j,1) = mean(prod(prod(yFSDx_income_happiness,2),3));
    overall_xFSDy_income_happiness_intersection(j,1) = mean(prod(prod(xFSDy_income_happiness_intersection,2),3));
    overall_yFSDx_income_happiness_intersection(j,1) = mean(prod(prod(yFSDx_income_happiness_intersection,2),3));
    
    %bivariate income and happiness
    xSSDy_income_happiness = H_income_happiness_x<=H_income_happiness_y_perm;
    ySSDx_income_happiness = H_income_happiness_y_perm<=H_income_happiness_x;
    xSSDy_income_happiness = double(xSSDy_income_happiness);
    ySSDx_income_happiness = double(ySSDx_income_happiness);
    xSSDy_income_happiness_intersection = xSSDy_income_happiness(:,1:poverty_line.happiness,1:poverty_line.income);
    ySSDx_income_happiness_intersection = ySSDx_income_happiness(:,1:poverty_line.happiness,1:poverty_line.income);
    prop_xSSDy_income_happiness(j,:,:) = mean(xSSDy_income_happiness);
    prop_ySSDx_income_happiness(j,:,:) = mean(ySSDx_income_happiness);
    overall_xSSDy_income_happiness(j,1) = mean(prod(prod(xSSDy_income_happiness,2),3));
    overall_ySSDx_income_happiness(j,1) = mean(prod(prod(ySSDx_income_happiness,2),3));
    overall_xSSDy_income_happiness_intersection(j,1) = mean(prod(prod(xSSDy_income_happiness_intersection,2),3));
    overall_ySSDx_income_happiness_intersection(j,1) = mean(prod(prod(ySSDx_income_happiness_intersection,2),3));

    % U2 dominance result
    overall_xFSDy_U2(j,1) = mean(prod(prod(xFSDy_income_health,2),3).*prod(prod(xFSDy_income_educ,2),3).*prod(prod(xFSDy_income_happiness,2),3));
    overall_yFSDx_U2(j,1) = mean(prod(prod(yFSDx_income_health,2),3).*prod(prod(yFSDx_income_educ,2),3).*prod(prod(yFSDx_income_happiness,2),3));
    overall_xFSDy_U2_intersection(j,1) = mean(prod(prod(xFSDy_income_health_intersection,2),3).*prod(prod(xFSDy_income_educ_intersection,2),3).*prod(prod(xFSDy_income_happiness_intersection,2),3));
    overall_yFSDx_U2_intersection(j,1) = mean(prod(prod(yFSDx_income_health_intersection,2),3).*prod(prod(yFSDx_income_educ_intersection,2),3).*prod(prod(yFSDx_income_happiness_intersection,2),3));
    
    %SSD only consider income and health compensated pronciple
      
    overall_xSSDy_income_health_compensated(j,1) = mean(prod(xSSDy_health,2).*prod(prod(xSSDy_income_health,2),3));
    overall_ySSDx_income_health_compensated(j,1) = mean(prod(ySSDx_health,2).*prod(prod(ySSDx_income_health,2),3));
    overall_xSSDy_income_health_compensated_intersection(j,1) = mean(prod(xSSDy_health(:,1:poverty_line.health),2).*prod(prod(xSSDy_income_health_intersection,2),3));
    overall_ySSDx_income_health_compensated_intersection(j,1) = mean(prod(ySSDx_health(:,1:poverty_line.health),2).*prod(prod(ySSDx_income_health_intersection,2),3));

      %SSD only consider income and educ compensated principle
      overall_xSSDy_income_educ_compensated(j,1) = mean(prod(xSSDy_educ(:,1:end-1),2).*prod(prod(xSSDy_income_educ,2),3));
      overall_ySSDx_income_educ_compensated(j,1) = mean(prod(ySSDx_educ(:,1:end-1),2).*prod(prod(ySSDx_income_educ,2),3));
      overall_xSSDy_income_educ_compensated_intersection(j,1) = mean(prod(xSSDy_educ(:,1:poverty_line.educ),2).*prod(prod(xSSDy_income_educ_intersection,2),3));
      overall_ySSDx_income_educ_compensated_intersection(j,1) = mean(prod(ySSDx_educ(:,1:poverty_line.educ),2).*prod(prod(ySSDx_income_educ_intersection,2),3));
       
       
      %SSD only consider income and happiness compensated principle
      overall_xSSDy_income_happiness_compensated(j,1) = mean(prod(xSSDy_happiness(:,1:end-1),2).*prod(prod(xSSDy_income_happiness,2),3));
      overall_ySSDx_income_happiness_compensated(j,1) = mean(prod(ySSDx_happiness(:,1:end-1),2).*prod(prod(ySSDx_income_happiness,2),3));
      overall_xSSDy_income_happiness_compensated_intersection(j,1) = mean(prod(xSSDy_happiness(:,1:poverty_line.happiness),2).*prod(prod(xSSDy_income_happiness_intersection,2),3));
      overall_ySSDx_income_happiness_compensated_intersection(j,1) = mean(prod(ySSDx_happiness(:,1:poverty_line.happiness),2).*prod(prod(ySSDx_income_happiness_intersection,2),3));
                    
      %SSD consider income educ health happiness compensated principle
      overall_xSSDy_income_health_educ_happiness_compensated(j,1) = mean(prod(xSSDy_health,2).*prod(xSSDy_educ(:,1:end-1),2).*prod(xSSDy_happiness(:,1:end-1),2).*...
                                                                prod(prod(xSSDy_income_health,2),3).*prod(prod(xSSDy_income_educ,2),3).*prod(prod(xSSDy_income_happiness,2),3));
      overall_ySSDx_income_health_educ_happiness_compensated(j,1) = mean(prod(ySSDx_health,2).*prod(ySSDx_educ(:,1:end-1),2).*prod(ySSDx_happiness(:,1:end-1),2).*...
                                                               prod(prod(ySSDx_income_health,2),3).*prod(prod(ySSDx_income_educ,2),3).*prod(prod(ySSDx_income_happiness,2),3));
      overall_xSSDy_income_health_educ_happiness_compensated_lowest(j,1) = mean(prod(xSSDy_health(:,1:poverty_line.health),2).*prod(xSSDy_educ(:,1:poverty_line.educ),2).*prod(xSSDy_happiness(:,1:poverty_line.happiness),2).*...
          prod(prod(xSSDy_income_health_intersection,2),3).*prod(prod(xSSDy_income_educ_intersection,2),3).*prod(prod(xSSDy_income_happiness_intersection,2),3));
      overall_ySSDx_income_health_educ_happiness_compensated_lowest(j,1) = mean(prod(ySSDx_health(:,1:poverty_line.health),2).*prod(ySSDx_educ(:,1:poverty_line.educ),2).*prod(ySSDx_happiness(:,1:poverty_line.happiness),2).*...
          prod(prod(ySSDx_income_health_intersection,2),3).*prod(prod(ySSDx_income_educ_intersection,2),3).*prod(prod(ySSDx_income_happiness_intersection,2),3));
    
    
end
 

save(['indices_result_',num2str(year_x),'_',num2str(year_y),'.mat'],'mean_income_x','mean_health_x','headcount_income_x','headcount_health_x','headcount_educ_x','headcount_happiness_x',...
                                            'multi_headcount_k1_x','multi_headcount_k2_x','multi_headcount_k3_x','multi_headcount_k4_x',...
                                            'multi_M0_k1_x','multi_M0_k2_x','multi_M0_k3_x','multi_M0_k4_x',...
                                            'FGT1_income_x','FGT2_income_x','FGT1_health_x','FGT2_health_x',...
                                            'CF01_education_x','CF09_education_x','CF01_happiness_x','CF09_happiness_x',...
                                            'mean_income_y','mean_health_y','headcount_income_y','headcount_health_y','headcount_educ_y','headcount_happiness_y',...
                                            'multi_headcount_k1_y','multi_headcount_k2_y','multi_headcount_k3_y','multi_headcount_k4_y',...
                                            'multi_M0_k1_y','multi_M0_k2_y','multi_M0_k3_y','multi_M0_k4_y',...
                                            'FGT1_income_y','FGT2_income_y','FGT1_health_y','FGT2_health_y',...
                                            'CF01_education_y','CF09_education_y','CF01_happiness_y','CF09_happiness_y',...
                                            'Gini_income_x','Gini_income_y','Gini_health_x','Gini_health_y');
                                        
save(['dominance_result_',num2str(year_x),'_',num2str(year_y),'.mat'],...
                                            'overall_xFSDy_income','overall_yFSDx_income','overall_xFSDy_income_lowest','overall_yFSDx_income_lowest',...
                                            'overall_xSSDy_income','overall_ySSDx_income','overall_xSSDy_income_lowest','overall_ySSDx_income_lowest',...
                                            'overall_xFSDy_health','overall_yFSDx_health','overall_xFSDy_health_lowest','overall_yFSDx_health_lowest',...
                                            'overall_xSSDy_health','overall_ySSDx_health','overall_xSSDy_health_lowest','overall_ySSDx_health_lowest',...
                                            'overall_xFSDy_educ','overall_yFSDx_educ','overall_xFSDy_educ_lowest','overall_yFSDx_educ_lowest',...
                                            'overall_xSSDy_educ','overall_ySSDx_educ','overall_xSSDy_educ_lowest','overall_ySSDx_educ_lowest',...
                                            'overall_xFSDy_happiness','overall_yFSDx_happiness','overall_xFSDy_happiness_lowest','overall_yFSDx_happiness_lowest',...
                                            'overall_xSSDy_happiness','overall_ySSDx_happiness','overall_xSSDy_happiness_lowest','overall_ySSDx_happiness_lowest',...
                                            'overall_xFSDy_income_health','overall_yFSDx_income_health','overall_xFSDy_income_health_intersection','overall_yFSDx_income_health_intersection',...
                                            'overall_xFSDy_income_educ','overall_yFSDx_income_educ','overall_xFSDy_income_educ_intersection','overall_yFSDx_income_educ_intersection',...
                                            'overall_xFSDy_income_happiness','overall_yFSDx_income_happiness','overall_xFSDy_income_happiness_intersection','overall_yFSDx_income_happiness_intersection',...
                                            'overall_xSSDy_income_health','overall_ySSDx_income_health','overall_xSSDy_income_health_intersection','overall_ySSDx_income_health_intersection',...
                                            'overall_xSSDy_income_educ','overall_ySSDx_income_educ','overall_xSSDy_income_educ_intersection','overall_ySSDx_income_educ_intersection',...
                                            'overall_xSSDy_income_happiness','overall_ySSDx_income_happiness','overall_xSSDy_income_happiness_intersection','overall_ySSDx_income_happiness_intersection',...
                                            'prop_xFSDy_income','prop_yFSDx_income',...
                                            'prop_xSSDy_income','prop_ySSDx_income',...
                                            'prop_xFSDy_health','prop_yFSDx_health',...
                                            'prop_xSSDy_health','prop_ySSDx_health',...
                                            'prop_xFSDy_educ','prop_yFSDx_educ',...
                                            'prop_xSSDy_educ','prop_ySSDx_educ',...
                                            'prop_xFSDy_happiness','prop_yFSDx_happiness',...
                                            'prop_xSSDy_happiness','prop_ySSDx_happiness',...
                                            'prop_xFSDy_income_health','prop_yFSDx_income_health',...
                                            'prop_xSSDy_income_health','prop_ySSDx_income_health',...
                                            'prop_xFSDy_income_educ','prop_yFSDx_income_educ',...
                                            'prop_xSSDy_income_educ','prop_ySSDx_income_educ',...
                                            'prop_xFSDy_income_happiness','prop_yFSDx_income_happiness',...
                                            'prop_xSSDy_income_happiness','prop_ySSDx_income_happiness',...
                                            'overall_xSSDy_income_health_compensated','overall_ySSDx_income_health_compensated','overall_xSSDy_income_health_compensated_intersection','overall_ySSDx_income_health_compensated_intersection',...
                                            'overall_xSSDy_income_educ_compensated','overall_ySSDx_income_educ_compensated','overall_xSSDy_income_educ_compensated_intersection','overall_ySSDx_income_educ_compensated_intersection',...
                                            'overall_xSSDy_income_happiness_compensated','overall_ySSDx_income_happiness_compensated','overall_xSSDy_income_happiness_compensated_intersection','overall_ySSDx_income_happiness_compensated_intersection',...
                                            'overall_xSSDy_income_health_educ_happiness_compensated','overall_ySSDx_income_health_educ_happiness_compensated',...
                                            'overall_xSSDy_income_health_educ_happiness_compensated_lowest','overall_ySSDx_income_health_educ_happiness_compensated_lowest',...
                                            'overall_xFSDy_U2','overall_yFSDx_U2',...
                                            'overall_xFSDy_U2_intersection','overall_yFSDx_U2_intersection');
                                                                                
                                        
                                        