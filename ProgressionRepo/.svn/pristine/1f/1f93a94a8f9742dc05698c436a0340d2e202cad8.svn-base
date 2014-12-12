function p = PrecedenceMatrix1(atrophy)
% Estimates a probability of precedence matrix p(Ri<Rj) from the atrophy
% scores of each subject in each region.  This is a very simple version
% in which each probability derives from the relative amount of atrophy in
% the two regions.
[numMeas numEvents] = size(atrophy);

% Compute probability matrix that event i occurs before event j.
p = ones(numEvents, numEvents);
% reg = 0.2;
% for i=1:(numEvents-1)
%     for j=(i+1):numEvents
%         p(i,j) = mean((atrophy(:,i)+reg)./(atrophy(:,i)+atrophy(:,j)+2*reg));
%         p(j,i) = 1-p(i,j);
%     end
% end
% Alternative p from counting number of times atrophy is greater.
% for i=1:numEvents
%     for j=1:numEvents
%         p(i,j)=sum(atrophy(:,i)>=atrophy(:,j))/numMeas;
%     end
% end
% Now one that accounts for how big the differences are.
ct = zeros(numEvents, numEvents);
for i=1:(numEvents-1)
    for j=(i+1):numEvents
        ct(i,j)=sum(atrophy(:,i)>=atrophy(:,j));
        ct(j,i) = numMeas - ct(i,j);
        dinds = find(atrophy(:,i)>=atrophy(:,j));
        dis(i,j) = sum(atrophy(dinds,i)-atrophy(dinds,j));
        dinds = find(atrophy(:,i)<atrophy(:,j));
        dis(j,i) = sum(atrophy(dinds,j)-atrophy(dinds,i));
        p(i,j) = ct(i,j)*dis(i,j)/(ct(i,j)*dis(i,j) + ct(j,i)*dis(j,i));
        p(j,i) = 1-p(i,j);
    end
end

% Replace zero probabilities with a small number
p(find(p==0)) = min(p(find(p>0)))/10;
