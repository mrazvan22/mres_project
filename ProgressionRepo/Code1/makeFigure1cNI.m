% Making figure 1c of the neuroimage revision2, making sure that this time,
% p(x | E) and p(x | noE) are in the right order and there's only a cutoff 
% in p(x | E)

clear all
close all

mu = [1.4 2.4 3.4 5.5];
sig = [1 1.2 1.4 1.1];
x = [1:6];

p_E = zeros(4, 6);
for e = 1:4,
    
    p_E(e, :) = normpdf(x, mu(e), sig(e));
    
end
p_noE = zeros(4, 6);
p_noE(1, round(mu(1)):end) = 0.71;
p_noE(2, round(mu(2)):end) = 0.6;
p_noE(3, round(mu(3)):end) = 0.55;
p_noE(4, round(mu(4)):end) = 0.65;

figure(1), clf
subplot(121), imagesc(p_E)
axis equal, axis off
subplot(122), imagesc(p_noE, [0 1])
axis equal, axis off
    