function p = PrecedenceMatrix2(ControlJacs, PatientJacs)
% Estimates a probability of precedence matrix p(Ri<Rj) from the atrophy
% scores of each patient in each region by comparison with the distribution
% of atrophy scores in controls.

mat = mean(ControlJacs);
sat = std(ControlJacs);
mmat = mean(ControlJacs(:));
msat = std(ControlJacs(:));

[numMeas numEvents] = size(PatientJacs);

% Probability of atrophy for each region in each patient data set assuming
% uniform distribution for jacobian given atrophy and using the normal
% model of unatrophied jacobians from controls.
pat = zeros(size(PatientJacs));
for i=1:numMeas
    for j=1:numEvents
        % Regional means and standard deviations
%        tp = exp(-(PatientJacs(i,j) - mat(j))^2/(2*sat(j)^2))/(sqrt(2*pi*sat(j)^2));
        % Pooled mean and standard deviation
        tp = exp(-(PatientJacs(i,j) - mmat)^2/(2*msat^2))/(sqrt(2*pi*msat^2));
        pat(i,j) = 1/(1+tp);
    end
end


% Probability of atrophy in the first of a pair, but not the second in each
% patient data set.
par1nar2 = zeros(numMeas, numEvents, numEvents);
for i=1:numMeas
    par1nar2(i,:,:) = pat(i,:)'*(1-pat(i,:));
end


% Now produce the p(R1<R2) matrix.
% May need to be careful about entries for which probability is close to
% zero everywhere, but this does not seem to be a problem at the moment.
% If so, add a regularizing factor to the denominator.
p = squeeze(sum(par1nar2.^2)./sum(par1nar2));

% Probabilities of simultaneous atrophy.  Here we use a simple but naive
% definition, but it should really account better for the distribution of
% samples along the L-curve.
%peq = (1 - p).*(1 - p');

