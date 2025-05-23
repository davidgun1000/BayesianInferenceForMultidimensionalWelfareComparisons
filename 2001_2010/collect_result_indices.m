load('indices_result_2001_2010_v2.mat');

load('param2001.mat');
mean_education_x = prop_educ(:,1)*1 + prop_educ(:,2)*2 + prop_educ(:,3)*3 + prop_educ(:,4)*4 + prop_educ(:,5)*5;
mean_happiness_x = prop_happiness(:,1)*1 + prop_happiness(:,2)*2 + prop_happiness(:,3)*3 + prop_happiness(:,4)*4 + prop_happiness(:,5)*5;

load('param2010.mat');
mean_education_y = prop_educ(:,1)*1 + prop_educ(:,2)*2 + prop_educ(:,3)*3 + prop_educ(:,4)*4 + prop_educ(:,5)*5;
mean_happiness_y = prop_happiness(:,1)*1 + prop_happiness(:,2)*2 + prop_happiness(:,3)*3 + prop_happiness(:,4)*4 + prop_happiness(:,5)*5;



mean(mean_income_y) - mean(mean_income_x)%Table 4
mean(mean_health_y) - mean(mean_health_x)%Table 4
mean(mean_education_y) - mean(mean_education_x)%Table 4
mean(mean_happiness_y) - mean(mean_happiness_x)%Table 4


load('MWI_2001_v2.mat');
MWI_income_health_x = MWI_income_health;
MWI_income_education_x = MWI_income_education;
MWI_income_happiness_x = MWI_income_happiness; 
biv_headcount_income_health_x = biv_headcount_income_health;
biv_headcount_income_education_x = biv_headcount_income_education;
biv_headcount_income_happiness_x = biv_headcount_income_happiness;


load('MWI_2010_v2.mat');
MWI_income_health_y = MWI_income_health;
MWI_income_education_y = MWI_income_education;
MWI_income_happiness_y = MWI_income_happiness; 
biv_headcount_income_health_y = biv_headcount_income_health;
biv_headcount_income_education_y = biv_headcount_income_education;
biv_headcount_income_happiness_y = biv_headcount_income_happiness;



mean(MWI_income_health_y) - mean(MWI_income_health_x)%Table 5
mean(MWI_income_education_y) - mean(MWI_income_education_x)%Table 5
mean(MWI_income_happiness_y) - mean(MWI_income_happiness_x)%Table 5


mean(headcount_income_x) - mean(headcount_income_y)%Table 6
mean(headcount_health_x) - mean(headcount_health_y)%Table 6
mean(headcount_educ_x) - mean(headcount_educ_y)%Table 6
mean(headcount_happiness_x) - mean(headcount_happiness_y)%Table 6

mean(biv_headcount_income_health_x) - mean(biv_headcount_income_health_y)%Table7
mean(biv_headcount_income_education_x) - mean(biv_headcount_income_education_y)%Table7
mean(biv_headcount_income_happiness_x) - mean(biv_headcount_income_happiness_y)%Table7

mean(multi_headcount_k1_x) - mean(multi_headcount_k1_y)%Table8
mean(multi_headcount_k2_x) - mean(multi_headcount_k2_y)%Table8
mean(multi_headcount_k3_x) - mean(multi_headcount_k3_y)%Table8
mean(multi_headcount_k4_x) - mean(multi_headcount_k4_y)%Table8

mean(multi_M0_k1_x) - mean(multi_M0_k1_y)%Table8
mean(multi_M0_k2_x) - mean(multi_M0_k2_y)%Table8
mean(multi_M0_k3_x) - mean(multi_M0_k3_y)%Table8
mean(multi_M0_k4_x) - mean(multi_M0_k4_y)%Table8



