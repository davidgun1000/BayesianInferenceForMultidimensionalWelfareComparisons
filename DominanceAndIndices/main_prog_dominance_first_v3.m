%mainprogfourdimension
num_FPBB=200;
year_x='2001';
year_y='2010';

num_component_gamma_x = 3;
num_component_beta_x = 3;

num_component_gamma_y = 2;
num_component_beta_y = 3;



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



grid_income = (4:0.0335:7.3132)';
grid_health = 0.01:0.01:0.99;
grid_income = (exp(grid_income));
grid_health = grid_health';
grid_educ = [1;2;3;4;5];
grid_happiness = [1;2;3;4;5];

poverty_line.income=40;
poverty_line.health=50;
poverty_line.educ=2;
poverty_line.happiness=2;

num_param = length(w_gamma_x);
F_income_health_educ_happiness_x = [];
F_income_health_educ_happiness_y = [];

parpool(48)
parfor j=1:num_param

    [output_x]=obtain_dominance_matrix_first(w_gamma_x(j,:),m_gamma_x(j,:),v_gamma_x(j,:),w_beta_x(j,:),m_beta_x(j,:),s_beta_x(j,:),theta_gauss_x(j,:),prop_educ_x(j,:),prop_happiness_x(j,:),poverty_line,grid_income,grid_health,grid_educ,grid_happiness,...
                num_component_gamma_x,num_component_beta_x);
    F_income_health_educ_happiness_x = [output_x.income_health_educ_happiness];
    
    [output_y]=obtain_dominance_matrix_first(w_gamma_y(j,:),m_gamma_y(j,:),v_gamma_y(j,:),w_beta_y(j,:),m_beta_y(j,:),s_beta_y(j,:),theta_gauss_y(j,:),prop_educ_y(j,:),prop_happiness_y(j,:),poverty_line,grid_income,grid_health,grid_educ,grid_happiness,...
                num_component_gamma_y,num_component_beta_y);
    
    F_income_health_educ_happiness_y = [output_y.income_health_educ_happiness];
    
    xFSDy_U1 =  F_income_health_educ_happiness_x<=F_income_health_educ_happiness_y;
    yFSDx_U1 =  F_income_health_educ_happiness_y<=F_income_health_educ_happiness_x;
    xFSDy_U1 = double(xFSDy_U1);
    yFSDx_U1 = double(yFSDx_U1);
    xFSDy_U1_intersection = xFSDy_U1(:,1:poverty_line.happiness,1:poverty_line.educ,1:poverty_line.health,1:poverty_line.income);
    yFSDx_U1_intersection = yFSDx_U1(:,1:poverty_line.happiness,1:poverty_line.educ,1:poverty_line.health,1:poverty_line.income);
    
    overall_xFSDy_U1_temp(j,1) = (prod(prod(prod(prod(xFSDy_U1,2),3),4),5));
    overall_yFSDx_U1_temp(j,1) = (prod(prod(prod(prod(yFSDx_U1,2),3),4),5));
    overall_xFSDy_U1_intersection_temp(j,1) = (prod(prod(prod(prod(xFSDy_U1_intersection,2),3),4),5));
    overall_yFSDx_U1_intersection_temp(j,1) = (prod(prod(prod(prod(yFSDx_U1_intersection,2),3),4),5));
    
    
end

overall_xFSDy_U1 = mean(overall_xFSDy_U1_temp);
overall_yFSDx_U1 = mean(overall_yFSDx_U1_temp);
overall_xFSDy_U1_intersection = mean(overall_xFSDy_U1_intersection_temp);
overall_yFSDx_U1_intersection = mean(overall_yFSDx_U1_intersection_temp);



% m=size(w_gamma_x,1);
% 
% for j=1:1000
%     %j
%     R = (randperm(m))';
%     for i=1:m
%         
%         F_income_health_educ_happiness_y_perm(i,:,:,:,:) = F_income_health_educ_happiness_y(R(i,1),:,:,:,:);
%         
%     end
%     
%     xFSDy_U1 =  F_income_health_educ_happiness_x<=F_income_health_educ_happiness_y_perm;
%     yFSDx_U1 =  F_income_health_educ_happiness_y_perm<=F_income_health_educ_happiness_x;
%     xFSDy_U1 = double(xFSDy_U1);
%     yFSDx_U1 = double(yFSDx_U1);
%     xFSDy_U1_intersection = xFSDy_U1(:,1:poverty_line.happiness,1:poverty_line.educ,1:poverty_line.health,1:poverty_line.income);
%     yFSDx_U1_intersection = yFSDx_U1(:,1:poverty_line.happiness,1:poverty_line.educ,1:poverty_line.health,1:poverty_line.income);
%     
%     overall_xFSDy_U1(j,1) = mean(prod(prod(prod(prod(xFSDy_U1,2),3),4),5));
%     overall_yFSDx_U1(j,1) = mean(prod(prod(prod(prod(yFSDx_U1,2),3),4),5));
%     overall_xFSDy_U1_intersection(j,1) = mean(prod(prod(prod(prod(xFSDy_U1_intersection,2),3),4),5));
%     overall_yFSDx_U1_intersection(j,1) = mean(prod(prod(prod(prod(yFSDx_U1_intersection,2),3),4),5));
%     
%     
%     
% end
 

                                        
save(['first_dominance_result_',num2str(year_x),'_',num2str(year_y),'.mat'],...
                                            'overall_xFSDy_U1','overall_yFSDx_U1',...
                                            'overall_xFSDy_U1_intersection','overall_yFSDx_U1_intersection');
                                                                                
                                        
                                        