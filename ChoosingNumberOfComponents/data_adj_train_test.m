%split the data to training and test sets.

year_x='2019';
load(['wellbeing_data_',num2str(year_x),'.mat']);
data_used = data_overall_2019;
weight_used = weight_overall_2019;
thres = 0.8;
unif = rand(length(data_used(:,1)),1);
id_train = (unif<0.8);
id_test = (unif>0.8);

data_train = data_used(id_train,:);
weight_train = weight_used(id_train,:);

data_test = data_used(id_test,:);
weight_test = weight_used(id_test,:);

save(['wellbeing_data_train_test_',num2str(year_x),'.mat'],'data_train','weight_train','data_test','weight_test');

