
%This is the code to calculate posterior means of education and happiness
%for years 2001, 2010, 2015, 2019
load('param2001.mat');
m = 10000;
educ = [1;2;3;4;5];
happiness = [1;2;3;4;5];

cumprop_educ = cumsum(prop_educ,2);
cumprop_happiness = cumsum(prop_happiness,2);
for i=1:m
    mean_educ_2001(i,1) = sum(prop_educ(i,:)'.*educ);
    mean_happiness_2001(i,1) = sum(prop_happiness(i,:)'.*happiness);
    median_educ_2001(i,1) = find(cumprop_educ(i,:) >= 0.5,1,'first');
    median_happiness_2001(i,1) = find(cumprop_happiness(i,:) >= 0.5,1,'first');
end

posmean_educ_2001 = mean(mean_educ_2001);
posstd_educ_2001 = std(mean_educ_2001);

posmean_happiness_2001 = mean(mean_happiness_2001);
posstd_happiness_2001 = std(mean_happiness_2001);

load('param2010.mat');
cumprop_educ = cumsum(prop_educ,2);
cumprop_happiness = cumsum(prop_happiness,2);
for i=1:m
    mean_educ_2010(i,1) = sum(prop_educ(i,:)'.*educ);
    mean_happiness_2010(i,1) = sum(prop_happiness(i,:)'.*happiness);
    median_educ_2010(i,1) = find(cumprop_educ(i,:) >= 0.5,1,'first');
    median_happiness_2010(i,1) = find(cumprop_happiness(i,:) >= 0.5,1,'first');
end

posmean_educ_2010 = mean(mean_educ_2010);
posstd_educ_2010 = std(mean_educ_2010);

posmean_happiness_2010 = mean(mean_happiness_2010);
posstd_happiness_2010 = std(mean_happiness_2010);




load('param2015.mat');
cumprop_educ = cumsum(prop_educ,2);
cumprop_happiness = cumsum(prop_happiness,2);
for i=1:m
    mean_educ_2015(i,1) = sum(prop_educ(i,:)'.*educ);
    mean_happiness_2015(i,1) = sum(prop_happiness(i,:)'.*happiness);
    median_educ_2015(i,1) = find(cumprop_educ(i,:) >= 0.5,1,'first');
    median_happiness_2015(i,1) = find(cumprop_happiness(i,:) >= 0.5,1,'first');
end

posmean_educ_2015 = mean(mean_educ_2015);
posstd_educ_2015 = std(mean_educ_2015);

posmean_happiness_2015 = mean(mean_happiness_2015);
posstd_happiness_2015 = std(mean_happiness_2015);



load('param2019.mat');
cumprop_educ = cumsum(prop_educ,2);
cumprop_happiness = cumsum(prop_happiness,2);
for i=1:m
    mean_educ_2019(i,1) = sum(prop_educ(i,:)'.*educ);
    mean_happiness_2019(i,1) = sum(prop_happiness(i,:)'.*happiness);
    median_educ_2019(i,1) = find(cumprop_educ(i,:) >= 0.5,1,'first');
    median_happiness_2019(i,1) = find(cumprop_happiness(i,:) >= 0.5,1,'first');

end

posmean_educ_2019 = mean(mean_educ_2019);
posstd_educ_2019 = std(mean_educ_2019);

posmean_happiness_2019 = mean(mean_happiness_2019);
posstd_happiness_2019 = std(mean_happiness_2019);


