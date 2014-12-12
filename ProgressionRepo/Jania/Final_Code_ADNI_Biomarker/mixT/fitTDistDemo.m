function r = fitTDistDemo

%close all previous figures
close all;


%generate some t-ish data
N_DATA = 1000;
actNu = 1.5;actMu = 12; actSigma = 0.8;
actU = gamrnd(actNu/2,actNu/2,1,N_DATA);
actSD = sqrt(actSigma./actU);
data = randn(1,N_DATA).*actSD+actMu;


figure;set(gcf,'Color',[1 1 1]);
x = 1:0.5:20;
h=hist(data,x);
bar(x,h/(N_DATA*0.5),1.0);hold on;set(gca,'Box','Off');

%fit t-distribution
[estMu estSigma estNu] = fitTDist(data)


x2 = 1:0.1:20;
plot(x2,tpdf((x2-estMu)/sqrt(estSigma),estNu)/sqrt(estSigma),'r-','LineWidth',2);

%===============================================
%main loop for EM algorithm

function [estMu estSigma estNu] = fitTDist(data);

N_ITER = 30;
N_DATA = size(data,2);
%set initial mean, variance, dof
estMu = mean(data,2);
estSigma = std(data,[],2);
estNu = rand*10+1;

expectedU = zeros(1,N_DATA);
expectedLogU = zeros(1,N_DATA);
for (cIter = 1:N_ITER)
    x = 1:0.5:20;hold off;
    h=hist(data,x);
    bar(x,h/(N_DATA*0.5),1.0);hold on;set(gca,'Box','Off');

    x2 = 1:0.1:20;
    plot(x2,tpdf((x2-estMu)/sqrt(estSigma),estNu)/sqrt(estSigma),'r-','LineWidth',2);
    drawnow;   
    
    %E-Step - calculate expected U and expected Log U
    a = repmat(estNu/2+1/2,1,N_DATA);
    b = estNu/2+0.5*(((data-estMu).^2)/estSigma);
    expectedU = a./b;
    expectedLogU = psi(a)-log(b);
        
    %M-Step - calculate mean and variance 
    estMu = sum(expectedU.*data)/sum(expectedU);
    estSigma = mean(((data-estMu).^2).*expectedU);
    
    %minimization function here
    estNu= fminbnd(@(x) optCrit(x,expectedLogU,expectedU),-1,5);
    estNu=exp(estNu);

    fprintf('End of Iter %d, Mean = %4.3f, Var = %4.3f, DOF = %4.3f\n',cIter,estMu,estSigma,estNu);

end;

%==================================
function r=optCrit(logNu,expectedLogU, expectedU)

nu = exp(logNu);
N = length(expectedU);
r = N*nu*log(nu/2)/2-N*log(gamma(nu/2))+(nu/2-1)*sum(expectedLogU)-nu*sum(expectedU)/2;
r = r*-1;

