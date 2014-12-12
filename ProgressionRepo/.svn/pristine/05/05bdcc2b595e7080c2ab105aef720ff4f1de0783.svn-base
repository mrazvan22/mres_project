function demo

close all
%load('data.mat');
%DEFINE TRUE PARAMETERS FOR MIXTURE OF TWO GAUSSIANS
gaussTrue.nGauss = 2;
gaussTrue.p = [0.3 0.7];
gaussTrue.mean = [-1 0.5];
gaussTrue.std = [0.5 0.25];

%GENERATE DATA FROM MIXTURE OF GAUSSIANS
nData = 400;
data = zeros(nData,1);
%for each data point
for (cData = 1:nData)
    %draw random number between zero and one
    rand01 = rand(1);
    %if less than pre-defined value then the data comes from first gaussian
    if (rand01<gaussTrue.p(1))
        data(cData) = randn*gaussTrue.std(1)+gaussTrue.mean(1);
        true_label(cData)=1;
    else%otherwise it comes from second gaussian
        data(cData) = randn*gaussTrue.std(2)+gaussTrue.mean(2);
        true_label(cData)=2;
    end;
end;

[label, model, llh] = emgm(data',2); 
%spread(data',label);


% draw the true and the estimated distributions
class1_trueProb=getGaussProb(data(true_label==1),gaussTrue.mean(1),gaussTrue.std(1));
class2_trueProb=getGaussProb(data(true_label==2),gaussTrue.mean(2),gaussTrue.std(2));

class1_estProb=getGaussProb(data(label==1),model.mu(1),model.Sigma(1));
class2_estProb=getGaussProb(data(label==2),model.mu(2),model.Sigma(2));

subplot(121), hold on
plot(data(true_label==1),class1_trueProb,'r.');
plot(data(true_label==2),class2_trueProb,'g.');

subplot(122), hold on
plot(data(label==1),class1_estProb,'r.');
plot(data(label==2),class2_estProb,'g.');






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%subroutine to return gaussian probabilities
function prob = getGaussProb(x,mean,sd)

prob = exp(-0.5*((x-mean).^2)/(sd*sd));
prob = prob/ sqrt(2*pi*sd*sd);