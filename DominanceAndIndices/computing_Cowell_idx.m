function [Cowell_idx]=computing_Cowell_idx(prop,alpha,numcat)
   
   length_alpha=length(alpha); 
   for i=1:length_alpha 
   if alpha(1,i)==0 
      for k=1:numcat 
          temp_Cowell_idx(1,k) = -prop(1,k).*log(sum(prop(1,1:k)));  
      end
      Cowell_idx(1,i)=sum(temp_Cowell_idx,2);       
   else
      for k=1:numcat 
          temp_Cowell_idx(1,k) = prop(1,k).*(sum(prop(1,1:k)).^alpha(1,i));  
      end
      Cowell_idx(1,i)=(1/(alpha(1,i)*(alpha(1,i)-1)))*(sum(temp_Cowell_idx,2)-1); 
   end
   end
    



end