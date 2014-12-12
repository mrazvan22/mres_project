function [logProbXGivenCase]=findLikelihoodCurrentOrdering(likelihood,event_order_current)

[nr_events nr_patients dummy]=size(likelihood);
pr_k=ones(1,nr_events)/nr_events;

I_select=[1:nr_patients];
probXjGivenCase=zeros(1,length(I_select));

likelihood(likelihood==0)=min(likelihood(likelihood~=0));

jCount=0;
%for each patient
for j=1:nr_patients
    jCount=jCount+1;
    probXjGivenSK=[];
    % for each stage
    for k=1: nr_events
        mCount=0;
        nCount=0;
        probXjGivenSK_event=[];
        probXjGivenSK_no_event=[];
        
        for m=1:k
            mCount=mCount+1;
            probXjGivenSK_event(mCount)=likelihood(event_order_current(m),j,2);
        end
        for n=k+1:nr_events
            nCount=nCount+1;
            probXjGivenSK_no_event(nCount)=likelihood(event_order_current(n),j,1);
        end
        probXjGivenSK(k)=sum(log(probXjGivenSK_event))+sum(log(probXjGivenSK_no_event))+log(pr_k(k));
        
        
    end % end stage k
    
    probXjGivenCase(jCount)=sum(exp(probXjGivenSK));
end % end patient j
temp=probXjGivenCase;
logProbXGivenCase=sum(log(temp));


