%this function generates the pseudo representative sample. 
%the weights are provided by HILDA survey.

function [obs_FPBB]=generate_pseudo_representative(income_data,sampling_weight,n)

    l_boot=zeros(n,1);
    
    weight_selected_sample = sampling_weight;
    N=sum(weight_selected_sample);
    Nnn = (N-n)/n;
    
    com = [income_data weight_selected_sample];
    length_FPBB= round(N-n);
    for k=1:length_FPBB
         
         newweights_num=com(:,5)-1+l_boot.*Nnn;
         newweights_den=(N-n)+(k-1)*Nnn;
         newweights=newweights_num./newweights_den;
         idx=newweights<0;
         newweights(idx,1)=0;
         [y_select(k,:),idx]=datasample(com(:,1:4),1,'replace',true,'weights',newweights);
         lk=zeros(n,1);
         lk(idx,1)=1;
         l_boot=l_boot+lk;
     end
     
     FPBB_population=[com(:,1:4);y_select];
     obs_FPBB=datasample(FPBB_population(:,1:4),n,'replace',true);


end