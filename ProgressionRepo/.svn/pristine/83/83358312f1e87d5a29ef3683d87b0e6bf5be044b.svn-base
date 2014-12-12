function logProbXGivenCaseMuSigma=findLikelihoodCurrentOrderingAndParamsZeroK(likelihood,event_order_current,prior_thisMean_control,prior_thisMean_patient,prior_thisCov_control,prior_thisCov_patient)

[nr_events nr_patients dummy]=size(likelihood);
pr_k=ones(1,nr_events+1)/(nr_events+1);

I_select=[1:nr_patients];
probXjGivenCase=zeros(1,length(I_select));
jCount=0;

likelihood(likelihood==0)=min(likelihood(likelihood~=0));

%for this ordering

%for each patient
for j=I_select
    jCount=jCount+1;
    probXjGivenSK=[];
    % for each stage
    for k=0: nr_events
        mCount=0;
        nCount=0;
        probXjGivenSK_event=[];
        probXjGivenSK_no_event=[];
        
        for m=1:k
            mCount=mCount+1;
            probXjGivenSK_event(mCount)=likelihood(event_order_current(m),j,2)*prior_thisMean_patient(event_order_current(m))*prior_thisCov_patient(event_order_current(m));
        end
        for n=k+1:nr_events
            nCount=nCount+1;
            probXjGivenSK_no_event(nCount)=likelihood(event_order_current(n),j,1)*prior_thisMean_control(event_order_current(n))*prior_thisCov_control(event_order_current(n));
        end
        probXjGivenSKMuSigma(k+1)=sum(log(probXjGivenSK_event))+sum(log(probXjGivenSK_no_event));
        
        
    end % end stage k
    
    probXjGivenCaseMuSigma(jCount)=sum(exp(probXjGivenSKMuSigma));
end % end patient j
temp=probXjGivenCaseMuSigma;
logProbXGivenCaseMuSigma=sum(log(temp));


% %scale the logs to avoid exp=0;
% mx=max(logProbXGivenCaseMuSigma);
% logProbXGivenCaseScale=logProbXGivenCaseMuSigma-mx;
% 
% PrSGivenXMuSigma=exp(logProbXGivenCaseScale)/sum(exp(logProbXGivenCaseScale));
% est_order_log_like=logProbXGivenCaseMuSigma;
