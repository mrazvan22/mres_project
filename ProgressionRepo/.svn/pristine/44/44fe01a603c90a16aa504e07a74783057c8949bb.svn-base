function [order_indx]=findNmostLikelyOrderings(posterior,n)

if n==1
[order_indx]=find(posterior==max(posterior));
else
    [order_indx]=find(posterior==max(posterior));
    posterior(order_indx)=realmin;
    [order_indx]=[order_indx findNmostLikelyOrderings(posterior,n-1)];
end
