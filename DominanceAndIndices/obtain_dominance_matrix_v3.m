function [output]=obtain_dominance_matrix_v3(w_gamma,m_gamma,v_gamma,w_beta,m_beta,s_beta,theta_gauss,prop_educ,prop_happiness,poverty_line,grid_income,grid_health,grid_educ,grid_happiness,...
                  num_component_gamma,num_component_beta)

   num_rep=100000;%the number of observations generated from the copula model for estimating indices and dominance probabilities
   dim_cop=4; 
   num_category_educ = length(prop_educ'); 
   num_category_happiness = length(prop_happiness');
   
   %construct the correlation matrix
   chol_gauss=diag([ones(dim_cop,1)]);
   number_gauss = 1;
   for i=1:dim_cop-1
       for j=i:dim_cop-1
           chol_gauss(i,j+1) = theta_gauss(1,number_gauss); 
           number_gauss=number_gauss+1;
       end
   end
   sigma = chol_gauss'*chol_gauss;
   diag_mat = diag([1./sqrt(diag(sigma))]);
   covmat = diag_mat*sigma*diag_mat;  

   %generate observations from Gaussian copula
   cop_draws = mvnrnd(zeros(1,dim_cop),covmat,num_rep); 
   cop_draws_uniform = normcdf(cop_draws); 
   num_rep_income = 100000;
   u1_income = rand(num_rep_income,1);
   obs_income = [];
   for k=1:num_component_gamma
       z_income = (u1_income>sum(w_gamma(1:(k-1)))) & (u1_income<=sum(w_gamma(1:k)));
       n_income = sum(z_income);
       obs_temp_income = gamrnd(v_gamma(k),m_gamma(k)/v_gamma(k),n_income,1);
       obs_income = [obs_income;obs_temp_income];
   end
   obs_income = sort(obs_income,1);
   cdf_income = sum(repmat(w_gamma(1:num_component_gamma), num_rep_income, 1) .* ...
       gamcdf(repmat(obs_income, 1, num_component_gamma), repmat(v_gamma(1:num_component_gamma), num_rep_income, 1),repmat(m_gamma(1:num_component_gamma), num_rep_income, 1)./repmat(v_gamma(1:num_component_gamma), num_rep_income, 1)), 2);
   cdf_income(end,1) = 1;
   table_income = [obs_income,cdf_income];
   %generate observations for income
   for k=1:num_rep
        indx_income = find(cop_draws_uniform(k,1)<=table_income(:,2),1,'first');
        if indx_income == 1
           data_dominance(k,1) = table_income(indx_income(1,1),1); 
        else
           data_dominance(k,1) = unifrnd(table_income(indx_income(1,1)-1,1),table_income(indx_income(1,1),1)); 
        end
   end
   
   num_rep_health=100000;
   u1_health = rand(num_rep_health,1);
   obs_health = [];
   for k=1:num_component_beta
       z_health = (u1_health>sum(w_beta(1:(k-1)))) & (u1_health<=sum(w_beta(1:k)));
       n_health = sum(z_health);
       obs_temp_health = betarnd(s_beta(k)*m_beta(k),s_beta(k)*(1-m_beta(k)),n_health,1); 
       obs_health = [obs_health;obs_temp_health];       
   end
   obs_health = sort(obs_health,1);
   
   cdf_health = sum(repmat(w_beta(1:num_component_beta), num_rep_health, 1) .* ...
       betacdf(repmat(obs_health, 1, num_component_beta), repmat(s_beta(1:num_component_beta), num_rep_income, 1).*repmat(m_beta(1:num_component_beta), num_rep_income, 1),...
       repmat(s_beta(1:num_component_beta), num_rep_income, 1).*(1-repmat(m_beta(1:num_component_beta), num_rep_income, 1))), 2);
   cdf_health(end,1) = 1;
   table_health = [obs_health,cdf_health];
   %generate observations for mental health
   for k=1:num_rep
       indx_health = find(cop_draws_uniform(k,2)<=table_health(:,2),1,'first');
       if indx_health ==1
          data_dominance(k,2) = table_health(indx_health(1,1),1);           
       else
          data_dominance(k,2) = unifrnd(table_health(indx_health(1,1)-1,1),table_health(indx_health(1,1),1));   
       end
   end
   
   %generate observations for education
   cum_marginal_educ = cumsum(prop_educ);
   cum_marginal_educ(1,end)=1;
   for k=1:num_category_educ
       if k==1
          id = cop_draws_uniform(:,3)<=cum_marginal_educ(k);
          data_dominance(id,3) = k; 
        else
          id = cum_marginal_educ(k-1)<=cop_draws_uniform(:,3) & cop_draws_uniform(:,3)<=cum_marginal_educ(k);
           data_dominance(id,3) = k;
        end
   end
   
   %generate observations for mental health
   cum_marginal_happiness = cumsum(prop_happiness);
   cum_marginal_happiness(1,end)=1;
   for k=1:num_category_happiness
       if k==1
          id = cop_draws_uniform(:,4)<=cum_marginal_happiness(k);
          data_dominance(id,4) = k; 
        else
          id = cum_marginal_happiness(k-1)<=cop_draws_uniform(:,4) & cop_draws_uniform(:,4)<=cum_marginal_happiness(k);
           data_dominance(id,4) = k;
        end
   end
   
   %quantities of interest
   
   output.mean_income = sum(w_gamma(1:num_component_gamma).*m_gamma(1:num_component_gamma));%estimate mean income
   output.mean_health = sum(w_beta(1:num_component_beta).*m_beta(1:num_component_beta));%estimate mean mental health score
   
   output.headcount_income = mean(data_dominance(:,1)<=grid_income(poverty_line.income,1));%estimate headcount for income
   output.headcount_health = mean(data_dominance(:,2)<=grid_health(poverty_line.health,1)); %estimate headcount for health
   output.headcount_educ = mean(data_dominance(:,3)<=grid_educ(poverty_line.educ,1));%estimate headcount for education
   output.headcount_happiness = mean(data_dominance(:,4)<=grid_happiness(poverty_line.happiness,1));%estimate headcount for happiness
   
   output.FGT1_income = mean((data_dominance(:,1)<=grid_income(poverty_line.income,1)).*((grid_income(poverty_line.income,1) - data_dominance(:,1))./grid_income(poverty_line.income,1)));%estimate FGT1 index for incomr
   output.FGT2_income = mean((data_dominance(:,1)<=grid_income(poverty_line.income,1)).*(((grid_income(poverty_line.income,1) - data_dominance(:,1))./grid_income(poverty_line.income,1)).^2));%estimate FGT2 index for income
   
   output.FGT1_health = mean((data_dominance(:,2)<=grid_health(poverty_line.health,1)).*((grid_health(poverty_line.health,1) - data_dominance(:,2))./grid_health(poverty_line.health,1)));%estimate FGT1 index for health
   output.FGT2_health = mean((data_dominance(:,2)<=grid_health(poverty_line.health,1)).*(((grid_health(poverty_line.health,1) - data_dominance(:,2))./grid_health(poverty_line.health,1)).^2));%estimate FGT2 index for health
   
   output.CF01_education = computing_Cowell_idx(prop_educ(1,:),0.1,num_category_educ);%estimate Cowell-Flachaire index for education
   output.CF09_education = computing_Cowell_idx(prop_educ(1,:),0.9,num_category_educ);
   
   output.CF01_happiness = computing_Cowell_idx(prop_happiness(1,:),0.1,num_category_happiness);%estimate Cowell-Flachaire for happiness
   output.CF09_happiness = computing_Cowell_idx(prop_happiness(1,:),0.9,num_category_happiness);
   
   output.Gini_income = samplegini(data_dominance(:,1));%estimate sample Gini index
   output.Gini_health = samplegini(data_dominance(:,2));
   
   cutoff_k1=1;
   cutoff_k2=2;
   cutoff_k3=3;
   cutoff_k4=4;

   id_poor_income = data_dominance(:,1)<=grid_income(poverty_line.income,1);
   id_poor_health = data_dominance(:,2)<=grid_health(poverty_line.health,1);
   id_poor_educ = data_dominance(:,3)<=grid_educ(poverty_line.educ,1);
   id_poor_happiness = data_dominance(:,4)<=grid_happiness(poverty_line.happiness,1);
   
   %estimate multidimensional headcount index
   output.multi_headcount_k1 = mean((id_poor_income+id_poor_health+id_poor_educ+id_poor_happiness)>=cutoff_k1);
   output.multi_headcount_k2 = mean((id_poor_income+id_poor_health+id_poor_educ+id_poor_happiness)>=cutoff_k2);
   output.multi_headcount_k3 = mean((id_poor_income+id_poor_health+id_poor_educ+id_poor_happiness)>=cutoff_k3);
   output.multi_headcount_k4 = mean((id_poor_income+id_poor_health+id_poor_educ+id_poor_happiness)>=cutoff_k4);
   
   depr_matrix = [id_poor_income,id_poor_health,id_poor_educ,id_poor_happiness];
   depr_count = sum(depr_matrix,2);
   
   %estimate M0 index
   
   id_mult_poor_k1 = depr_count>=cutoff_k1;
   output.multi_M0_k1 = (sum(sum(depr_matrix.*id_mult_poor_k1,2)))/(length(depr_matrix)*4);
   
   id_mult_poor_k2 = depr_count>=cutoff_k2;
   output.multi_M0_k2 = (sum(sum(depr_matrix.*id_mult_poor_k2,2)))/(length(depr_matrix)*4);
   
   id_mult_poor_k3 = depr_count>=cutoff_k3;
   output.multi_M0_k3 = (sum(sum(depr_matrix.*id_mult_poor_k3,2)))/(length(depr_matrix)*4);
   
   id_mult_poor_k4 = depr_count>=cutoff_k4;
   output.multi_M0_k4 = (sum(sum(depr_matrix.*id_mult_poor_k4,2)))/(length(depr_matrix)*4);
   
   length_grid_income = length(grid_income);
   length_grid_health = length(grid_health);
   
   %count_income_health = 0;
   %count_income_educ = 0;
   %count_income_happiness=0;
      
   %H_income_health_y(count1,1) = mean((pov_line_income - income_y_draw).*ind_poor_income_y.*ind_poor_health_y) + ...
   %              mean((grid_income(1,1) - income_y_draw).*ind_poor_income_y_first.*ind_poor_health_y_first) - ...
   %              mean((grid_income(1,1) - income_y_draw).*ind_poor_income_y_first) - ...
   %              mean((grid_income(1,1) - income_y_draw).*ind_poor_health_y_first);
   
   %estimate the required quantities for computing dominance probabilities
   
   for ss1 = 1:length_grid_income
       ind_poor_income = (data_dominance(:,1)<=grid_income(ss1,1));
       ind_poor_income_y_first = (data_dominance(:,1)<=grid_income(1,1));
       
       output.F_income(ss1,1) = mean(ind_poor_income);   
       output.H_income(ss1,1) = mean(((grid_income(ss1,1) - data_dominance(:,1))).*(data_dominance(:,1)<=grid_income(ss1,1)))-...
           mean((grid_income(1,1) - data_dominance(:,1)).*ind_poor_income_y_first);
       
       
   end
   
   
   for ss2 = 1:length_grid_health
       ind_poor_health = data_dominance(:,2)<=grid_health(ss2,1);
       ind_poor_health_y_first = (data_dominance(:,2)<=grid_health(1,1));
       output.F_health(ss2,1) = mean(ind_poor_health);
       output.H_health(ss2,1) = mean(((grid_health(ss2,1) - data_dominance(:,2))).*(data_dominance(:,2)<=grid_health(ss2,1)))-...
           mean((grid_health(1,1) - data_dominance(:,2)).*ind_poor_health_y_first);       
          
   end
   
   output.F_educ = cumsum(prop_educ,2);
   output.H_educ = cumsum(output.F_educ,2);
   output.F_happiness = cumsum(prop_happiness,2);
   output.H_happiness = cumsum(output.F_happiness,2);
   
   for ss1 = 1:length_grid_income
       for ss2 = 1:length_grid_health 
           ind_poor_income = data_dominance(:,1)<=grid_income(ss1,1);
           ind_poor_health = data_dominance(:,2)<=grid_health(ss2,1);  
           ind_poor_income_y_first = (data_dominance(:,1) <=grid_income(1,1));
          
           output.F_income_health(1,ss2,ss1) = mean(ind_poor_income.*ind_poor_health);
           output.H_income_health(1,ss2,ss1) = mean((grid_income(ss1,1) - data_dominance(:,1)).*(data_dominance(:,1)<=grid_income(ss1,1)).*(data_dominance(:,2)<=grid_health(ss2,1))) - ...                          
                    mean((grid_income(1,1) - data_dominance(:,1)).*ind_poor_income_y_first.*ind_poor_health); 
           
       end
       
       for ss3=1:num_category_educ
           ind_poor_income = data_dominance(:,1)<=grid_income(ss1,1);
           ind_poor_educ = data_dominance(:,3)<=grid_educ(ss3,1);
           ind_poor_income_y_first = (data_dominance(:,1)<=grid_income(1,1));
           output.F_income_educ(1,ss3,ss1) = mean(ind_poor_income.*ind_poor_educ);
           output.H_income_educ(1,ss3,ss1) = mean((grid_income(ss1,1) - data_dominance(:,1)).*(data_dominance(:,1)<=grid_income(ss1,1)).*(data_dominance(:,3)<=grid_educ(ss3,1))) - ...
               mean((grid_income(1,1) - data_dominance(:,1)).*ind_poor_income_y_first.*ind_poor_educ);
            
       end

       for ss4=1:num_category_happiness
           ind_poor_income = data_dominance(:,1)<=grid_income(ss1,1);
           ind_poor_happiness = data_dominance(:,4)<=grid_happiness(ss4,1);
           ind_poor_income_y_first = (data_dominance(:,1)<=grid_income(1,1));
           output.F_income_happiness(1,ss4,ss1) = mean(ind_poor_income.*ind_poor_happiness);   
           output.H_income_happiness(1,ss4,ss1) = mean((grid_income(ss1,1) - data_dominance(:,1)).*(data_dominance(:,1)<=grid_income(ss1,1)).*(data_dominance(:,4)<=grid_happiness(ss4,1))) - ...
               mean((grid_income(1,1) - data_dominance(:,1)).*ind_poor_income_y_first.*ind_poor_happiness);
            
       end
   end
   
   
   
end

%(w_gamma_x,m_gamma_x,v_gamma_x,w_beta_x,m_beta_x,s_beta_x,theta_gauss_x,prop_educ_x,prop_happiness_x);