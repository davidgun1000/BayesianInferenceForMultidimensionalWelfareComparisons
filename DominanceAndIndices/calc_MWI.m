%calculate bivariate headcount, bivariate MWI and overall MWI, for 2001, 2010, 2015, and 2019.
year_x='2001';

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

num_param = length(w_gamma_x);
num_component_gamma_x = size(m_gamma_x,2);
num_component_beta_x = size(m_beta_x,2);

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

parpool(48)
parfor j=1:num_param
   
   [output] = gen_obs(w_gamma_x(j,:),m_gamma_x(j,:),v_gamma_x(j,:),w_beta_x(j,:),m_beta_x(j,:),s_beta_x(j,:),theta_gauss_x(j,:),prop_educ_x(j,:),prop_happiness_x(j,:),...
                  poverty_line, grid_income,grid_health,grid_educ,grid_happiness, num_component_gamma_x,num_component_beta_x);
   
   biv_headcount_income_health(j,1) = output.biv_headcount_income_health;
   biv_headcount_income_education(j,1) = output.biv_headcount_income_education;
   biv_headcount_income_happiness(j,1) = output.biv_headcount_income_happiness;
   
   MWI_income_health(j,1) = output.MWI_income_health;
   MWI_income_education(j,1) = output.MWI_income_education;
   MWI_income_happiness(j,1) = output.MWI_income_happiness;
   MWI_4dim(j,1) = output.MWI_4dim;
              
              
          
end

save(['MWI_',num2str(year_x),'.mat'],'biv_headcount_income_health','biv_headcount_income_education','biv_headcount_income_happiness',...
                                    'MWI_income_health','MWI_income_education','MWI_income_happiness','MWI_4dim');

mean(biv_headcount_income_health)
std(biv_headcount_income_health)

mean(biv_headcount_income_education)
std(biv_headcount_income_education)
                                
mean(biv_headcount_income_happiness)
std(biv_headcount_income_happiness)

mean(MWI_income_health)
std(MWI_income_health)

mean(MWI_income_education)
std(MWI_income_education)

mean(MWI_income_happiness)
std(MWI_income_happiness)

mean(MWI_4dim)
std(MWI_4dim)




                                