function est_order=findNEventOrdering(likelihood)

[nr_events nr_patients dummy]=size(likelihood);

S=perms([1:nr_events]);
nr_cases=size(S,1);
pr_k=ones(1,nr_events)/nr_events;

I_select=[1:nr_patients];
probXjGivenCase=zeros(nr_cases,length(I_select));



%for each ordering
for c_case=1:nr_cases
    jCount=0;
    %for each patient
    for j=I_select
        jCount=jCount+1;    
        % for each stage
        for k=1: nr_events
            mCount=0;
            nCount=0;
            probXjGivenSK_event=[];
            probXjGivenSK_no_event=[];
            
            for m=1:k
                mCount=mCount+1;
                probXjGivenSK_event(mCount)=likelihood(S(c_case,m),j,2);
            end
            for n=k+1:nr_events
                nCount=nCount+1;
                probXjGivenSK_no_event(nCount)=likelihood(S(c_case,n),j,1);
            end
            probXjGivenSK(k)=sum(log(probXjGivenSK_event))+sum(log(probXjGivenSK_no_event));
            
            
        end % end stage k
        
        probXjGivenCase(c_case,jCount)=sum(exp(probXjGivenSK));
    end % end patient j
    temp=probXjGivenCase(c_case,:);
    temp(temp==0)=realmin;
    logProbXGivenCase(c_case)=sum(log(temp)); 
end % end c_case

%scale the logs to avoid exp=0;
mx=max(logProbXGivenCase);
logProbXGivenCaseScale=logProbXGivenCase-mx;

PrSGivenX=exp(logProbXGivenCaseScale)/sum(exp(logProbXGivenCaseScale));

%Choose the orgering with the highest probability
order_indx=find(PrSGivenX==max(PrSGivenX))
est_order=S(order_indx,:);
est_order_log_like=logProbXGivenCase(order_indx);
