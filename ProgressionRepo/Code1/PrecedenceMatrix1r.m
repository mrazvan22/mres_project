function p = PrecedenceMatrix1r(jacobians)
% Estimates a probability of precedence matrix p(Ri<Rj) from the atrophy
% scores of each subject in each region.  This is a very simple version
% in which each probability derives from the relative amount of atrophy in
% the two regions.
% This version negates the ordering of PrecedenceMatrix1, which is more
% logical when the atrophy scores are jacobians, because they are lower
% when atrophy is higher.
[numMeas numEvents] = size(jacobians);

% Compute probability matrix that event i occurs before event j.
p = ones(numEvents, numEvents);
ct = zeros(numEvents, numEvents);
for i=1:(numEvents-1)
    for j=(i+1):numEvents
        ct(i,j)=sum(jacobians(:,i)<=jacobians(:,j));
        ct(j,i) = numMeas - ct(i,j);
        dinds = find(jacobians(:,i)<=jacobians(:,j));
        dis(i,j) = sum(jacobians(dinds,j)-jacobians(dinds,i));
        dinds = find(jacobians(:,i)>jacobians(:,j));
        dis(j,i) = sum(jacobians(dinds,i)-jacobians(dinds,j));
        p(i,j) = ct(i,j)*dis(i,j)/(ct(i,j)*dis(i,j) + ct(j,i)*dis(j,i));
        p(j,i) = 1-p(i,j);
    end
end

% Replace zero probabilities with a small number
p(find(p==0)) = min(p(find(p>0)))/10;
