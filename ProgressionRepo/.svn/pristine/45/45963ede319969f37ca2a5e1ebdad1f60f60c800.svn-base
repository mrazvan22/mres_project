function findTwoEventOrdering(likelihood)

%% Problem 1 : Noticed that the problematic patients are the ouliers i.e. close to controls
% th reason is that the likelihood values for these patients using only the
% patient gaussian component is very cllose to zero hence the likelihood of
% event having occured ends up being close to zero and this highly
% influcensed the position of that event in the ordering 

%% Solution 1 to problem 1:  a hack which hubert used in his imaging_biomarker experiment
%which is to set the values to some small number
%likelihood(likelihood<1e-3)=1e-3;

%% solution 2 to problem 1: Fit a distribution to the data which acount for outlier e.g.
% a T-distribution


%% Problem 2: some biomarkers such as hippovol have a very low covariance e.g. 1e-7
%therefore they end of having higher likelihood than that of abeta142

%% solution 1 to problem 2:  set all low covariance to 1e-3
%this partially solved the problem but bring another low-cov biomarker to the top

%% initialize
[nr_events nr_patients dummy]=size(likelihood);
pr_k1=0.5;
pr_k2=0.5;

%% case 1 is when we assume S=(1,2)
for j=1:nr_patients
probXjGivenCase1(j)=exp(log(likelihood(1,j,2)) + log(likelihood(2,j,1))) + ...  
    exp(log(likelihood(1,j,2))+ log(likelihood(2,j,2)));
end

%% case 2 is when we assume S=(2,1)
for j=1:nr_patients
probXjGivenCase2(j)=exp(log(likelihood(2,j,2)) + log(likelihood(1,j,1))) + ...
    exp(log(likelihood(2,j,2)) + log(likelihood(1,j,2)));
end

probXjGivenCase1(probXjGivenCase1==0)=realmin;
probXjGivenCase2(probXjGivenCase2==0)=realmin;


I_select=[1:nr_patients];
logProbXGivenCase1=sum(log(probXjGivenCase1(I_select)));
logProbXGivenCase2=sum(log(probXjGivenCase2(I_select)));

%scale the logs to avoid exp = 0
mx=max([logProbXGivenCase1, logProbXGivenCase2]);
logProbXGivenCase1=logProbXGivenCase1-mx;
logProbXGivenCase2=logProbXGivenCase2-mx;

PrS1GivenX=exp(logProbXGivenCase1)/(exp(logProbXGivenCase1)+exp(logProbXGivenCase2))
PrS2GivenX=exp(logProbXGivenCase2)/(exp(logProbXGivenCase1)+exp(logProbXGivenCase2))

i

