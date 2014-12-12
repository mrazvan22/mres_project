function rt=getSingleTermTDist(x,L)

%function to calcualte the joint probability of x,h i,e
%Pr(xi,hi)=Pr(xi|hi)Pr(hi)

n_hidden=size(L,2);


for(j=1:n_hidden)
         term(j)=getLogTProb(x,L(j).mean,L(j).cov,L(j).dof)+log(L(j).prior);   
end
    
%rescale log probabilities by their maximim
maxTerm=max(term(:)); 
term=term-maxTerm;

%do summation
summedTerm=sum(exp(term));

%convert back to log probabilities
term=log(summedTerm)+maxTerm;

rt=term; 