function prob=getExpMixTProb(data,tDist)

prob=0;
N_DIST=length(tDist);
for i=1:N_DIST
prob=prob+ exp(getLogTProb(data,tDist(i).mean,tDist(i).cov,tDist(i).dof))*tDist(i).prior;
end