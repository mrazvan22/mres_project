function rt=version3GetProbModel(data, XP,library,singleTerms,cModel)


%find doubleTerm
    doubleTerm=version3getDoubleTermTDist(data(:,cModel), XP,library);

    %add doubleTerm to correct model
    logProbModel=(sum(singleTerms)+doubleTerm)-singleTerms(cModel);
    
 rt= logProbModel;