function rt=getDoubleTerm(X,Xp,L)

%function to calcualte the joint probability of x,xp,h i,e
%Pr(xi,xp,hi)=Pr(xi|hi)Pr(xp|hi)Pr(hi)

n_hidden=size(L,2);

for(j=1:n_hidden)
         term(j)=getLogTProb(X,L(j).mean,L(j).cov,L(j).dof)+getLogTProb(Xp,L(j).mean,L(j).cov,L(j).dof)+log(L(j).prior); 
 end

 %rescale log probabilities by their maximim
maxTerm=max(term(:)); 
term=term-maxTerm;

%do summation
summedTerm=sum(exp(term));

%convert back to log probabilities
term=log(summedTerm)+maxTerm;


rt=term;