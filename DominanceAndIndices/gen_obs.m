function [output] = gen_obs(w_gamma,m_gamma,v_gamma,w_beta,m_beta,s_beta,theta_gauss,prop_educ,prop_happiness,poverty_line,grid_income,grid_health,grid_educ,grid_happiness,...
                  num_component_gamma,num_component_beta)


   num_rep=100000;
   dim_cop=4;
   %num_component=3;
   num_category_educ = length(prop_educ');
   num_category_happiness = length(prop_happiness');
   
   
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
   
   for k=1:num_rep
       indx_health = find(cop_draws_uniform(k,2)<=table_health(:,2),1,'first');
       if indx_health ==1
          data_dominance(k,2) = table_health(indx_health(1,1),1);           
       else
          data_dominance(k,2) = unifrnd(table_health(indx_health(1,1)-1,1),table_health(indx_health(1,1),1));   
       end
   end
   
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

   output.biv_headcount_income_health = mean((data_dominance(:,1)<=grid_income(poverty_line.income,1)).*(data_dominance(:,2)<=grid_health(poverty_line.health,1)));
   output.biv_headcount_income_education = mean((data_dominance(:,1)<=grid_income(poverty_line.income,1)).*(data_dominance(:,3)<=grid_educ(poverty_line.educ,1)));
   output.biv_headcount_income_happiness = mean((data_dominance(:,1)<=grid_income(poverty_line.income,1)).*(data_dominance(:,4)<=grid_happiness(poverty_line.happiness,1)));

   normalised_data_dominance = data_dominance;
   normalised_data_dominance(:,1) = (log(normalised_data_dominance(:,1)) - log(min(normalised_data_dominance(:,1))))./(log(max(normalised_data_dominance(:,1))) - log(min(normalised_data_dominance(:,1))));
   normalised_data_dominance(:,3) = (normalised_data_dominance(:,3) - min(normalised_data_dominance(:,3)))./(max(normalised_data_dominance(:,3)) - min(normalised_data_dominance(:,3)));
   normalised_data_dominance(:,4) = (normalised_data_dominance(:,4) - min(normalised_data_dominance(:,4)))./(max(normalised_data_dominance(:,4)) - min(normalised_data_dominance(:,4)));
   
   output.MWI_income_health = mean(0.5*normalised_data_dominance(:,1) + 0.5*normalised_data_dominance(:,2));
   output.MWI_income_education = mean(0.5*normalised_data_dominance(:,1) + 0.5*normalised_data_dominance(:,3));
   output.MWI_income_happiness = mean(0.5*normalised_data_dominance(:,1) + 0.5*normalised_data_dominance(:,4));
   output.MWI_4dim = mean(0.25*normalised_data_dominance(:,1) + 0.25*normalised_data_dominance(:,2) + 0.25*normalised_data_dominance(:,3) + 0.25*normalised_data_dominance(:,4));
   

end