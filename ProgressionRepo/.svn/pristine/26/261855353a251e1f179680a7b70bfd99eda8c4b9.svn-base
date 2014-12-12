function rt=getProbModel(data,Xp,L,S,model_num)

% a function to find the probability of a given model


D=getDoubleTermTDist(data(:,model_num),Xp,L);



pr_model=(sum(S)+D)-S(model_num);
rt=pr_model;