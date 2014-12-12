function   tDist= mixGaussFitAppearanceTest2(data,N_DIST,N_APP,covType,nAttempts,dataIDMatrix);
%fits a Gaussian mixture model with N_GAUSS or fewer components
%
%covType = 0;  %unconstrained - individual covariance matrices
%covType = 1;  %unconstrained - shared covariance matrices
%covType = 2;  %diagonal - individual covariance matrices
%covType = 3;  %diagonal - shared covariance matrices
%covType = 4;  %spherical - individual covariance matrices
%covType = 5;  %spherical - shared covariance matrix



if (covType<2)
    DIAG_FLAG = 0;
else
    DIAG_FLAG = 1;
end

%retrive number of different individuals
[N_DATA N_PERSON] = size(dataIDMatrix);
sampleSizes=sum(dataIDMatrix);

%measure covariance of data.
%find smallest tolerable SD for fitted gauss - sd must be at last 10000 of
%the Sd in the smallest direction

[N_DIM, N_DATA] = size(data);
prDist = 1/N_DIST;
 
prHiddenVector=zeros(N_DIST*N_APP,N_DATA);

%initialize responsibilities
[clusters,meanVec] = kMeans(data,N_DIST*N_APP);
prHidden = zeros(N_DIST,N_APP,N_DATA);
for (cDist = 1:N_DIST*N_APP)
      prHiddenVector(cDist,find(clusters==cDist)) =1;
end;

prHidden=reshape(prHiddenVector,N_DIST,N_APP,N_DATA);
%initialize U values
a = 1.6/2;b = 1.6/2;
expectedU = ones(N_DIST,N_APP,N_DATA);
expectedLogU = (psi(a)-log(b))*ones(N_DIST,N_APP,N_DATA);

MAX_ITER = 20;
MIN_ITER = 20;
DRAW_FLAG = 1;

if (DRAW_FLAG&N_DIM==2)
    figure;
    randColors = 0.3+0.7*rand(N_DIST*N_APP,3);
end;

%initialize prior
alpha0 = 1;


cIter = 0;
while(1)
    cIter = cIter+1;
     if (cIter>MAX_ITER) 
        break;
    end;
   
    fprintf('%d...',cIter);
    
    
    %=======================================================================
    %============================= M-STEP ==================================
    %=======================================================================

    
    %calculate prior
    N = sum(prHidden,3);   
    %N=sum(extPrHidden);
    alphaAndN=alpha0+N;
    prior = (alphaAndN)/sum(alphaAndN(:));  %dirichlet prior on weights  
    
    
    for (cDist =1:N_DIST)
        for(cApp=1:N_APP)
            tDist(cDist,cApp).prior = prior(cDist,cApp);
        end
    end;
    
    %calculate mean
     
    for (cDist = 1:N_DIST)
        for(cApp=1:N_APP)
            tDist(cDist,cApp).mean = sum(repmat(squeeze(expectedU(cDist,cApp,:).*prHidden(cDist,cApp,:))',N_DIM,1).*data,2);
            tDist(cDist,cApp).mean = tDist(cDist,cApp).mean/sum(expectedU(cDist,cApp,:).*prHidden(cDist,cApp,:),3);
        end
    end;
     
    
    
    %calculate covariance 
    if (covType==0|covType==2|covType==4)%    - if individual
        for (cDist = 1:N_DIST)
            for(cApp=1:N_APP)
                xMinusMean = data-repmat(tDist(cDist,cApp).mean,1,N_DATA);
                if (covType==0) %full covariance
                    tDist(cDist,cApp).cov = zeros(N_DIM);
                    for(cData = 1:N_DATA)
                        tDist(cDist,cApp).cov = tDist(cDist,cApp).cov+prHidden(cDist,cApp,cData)*expectedU(cDist,cApp,cData)*xMinusMean(:,cData)*xMinusMean(:,cData)';
                    end;
                elseif(covType==2) % diagonal covariance
                    prHidden(prHidden==0)=realmin;
                    %tDist(cDist).cov = sum((repmat(extPrHidden(cDist,:).*expectedU(cDist,:),N_DIM,1).*xMinusMean).*xMinusMean,2);
                    cc=sum((repmat(prHidden(cDist,cApp,:).*expectedU(cDist,cApp,:),N_DIM,1).*xMinusMean).*xMinusMean,2);
                    cc(cc==0)=realmin;
                    tDist(cDist,cApp).cov=cc;
                else %(covType==4) % spherical covariance
                    tDist(cDist,cApp).cov = sum(sum((repmat(prHidden(cDist,cApp,:).*expectedU(cDist,cApp,:),N_DIM,1).*xMinusMean).*xMinusMean,2));
                    tDist(cDist,cApp).cov = tDist(cDist,cApp).cov*ones(N_DIM,1)/N_DIM;
                end;
                tDist(cDist,cApp).cov = tDist(cDist,cApp).cov/sum(prHidden(cDist,cApp,:));

            end
        end;
        
    else%pooled covariance
        if (covType==1)
            tDist(1).cov = zeros(N_DIM);
        else
            tDist(1).cov = zeros(N_DIM,1);
        end;
        
        for (cDist = 1:N_DIST)
            for(cApp=1:N_APP)
                xMinusMean = data-repmat(tDist(cDist,cApp).mean,1,N_DATA);
                if (covType==1) %full covariance
                    for(cData = 1:N_DATA)
                        tDist(1).cov = tDist(1).cov+prHidden(cDist,cApp,cData)*expectedU(cDist,cApp,cData)*xMinusMean(:,cData)*xMinusMean(:,cData)';
                    end;
                elseif(covType==3) % diagonal covariance
                    prHidden(prHidden==0)=realmin;
                    %                tDist(1).cov = tDist(1).cov+sum((repmat(extPrHidden(cDist,:).*expectedU(cDist,:),N_DIM,1).*xMinusMean).*xMinusMean,2);
                    cc=sum((repmat(squeeze(expectedU(cDist,cApp,:).*prHidden(cDist,cApp,:))',N_DIM,1).*xMinusMean).*xMinusMean,2);
                    cc(cc==0)=realmin;
                    tDist(1).cov = tDist(1).cov+cc;
                else %(covType==5) % spherical covariance
                    tDist(1).cov = tDist(1).cov+sum(sum((repmat(prHidden(cDist,cApp,:).*expectedU(cDist,cApp,:),N_DIM,1).*xMinusMean).*xMinusMean,2))/N_DIM;
                end;
            end
        end;
        %normalize
        tDist(1).cov = tDist(1).cov/N_DATA;
        
        %copy to remaining points
        for (cDist = 1:N_DIST)
            for(cApp=1:N_APP)
                tDist(cDist,cApp).cov = tDist(1).cov;
            end
        end;
    end;
    
     %update degrees of freedom
    for (cDist = 1:N_DIST)
        for(cApp=1:N_APP)
            estNu= fminbnd(@(x) optCrit(x,expectedLogU(cDist,cApp,:),expectedU(cDist,cApp,:),prHidden(cDist,cApp,:)),-1,5);
            tDist(cDist,cApp).dof=exp(estNu);
            %      tDist(cDist).dof=1.6;
        end

    end;
    
    %=======================================================================
    %========================== DRAW DATA ==================================
    %=======================================================================
    
    if (DRAW_FLAG&N_DIM==2)
        plot(data(1,:),data(2,:),'k.');
        for (c1 = 1:N_DIST)
            for (c2 = 1:N_APP)
                set(gcf,'Color',[1 1 1]);
                drawGaussianOutline(gauss(c1,c2).mean,gauss(c1,c2).cov,randColors(((c1-1)*N_APP)+c2,:));
                hold on;
            end
        end;
        plot(data(1,:),data(2,:),'k.');
        axis square;axis equal;
        axis off;
        hold off;drawnow;
    end;
%     %=======================================================================
%     %========================== DRAW MEANS =================================
%     %=======================================================================
%     
%     if (MEAN_DRAW_FLAG&N_DIM>2)
%         mCount=0;
%         close
%         for (c1 = 1:4)
%             set(gcf,'Color',[1 1 1]);
%             mCount=mCount+1;
%             subplot(4,2,mCount); imshow([reshape(gauss(c1,1).mean,80,60)]);
%             colorbar
%             mCount=mCount+1;
%             subplot(4,2,mCount); imshow(reshape(gauss(c1,2).mean,80,60));
%             colorbar
%         end
%         drawnow;
%         fileName=['Mat data/shared cov/iter',num2str(cIter) '.jpg'];
%         print ('-djpeg',fileName);
%         
%         %resclaing the range of the intensities so it wont come out black
%         im=reshape(gauss(c1,1).cov,80,60);
%         mx=max(max(im));
%         mn=min(min(im));
%         im2=im-mn;
%         im2=im2./(mx-mn);
%         close;
%         imshow(im2);
%         colorbar;
%         drawnow;
%         fileName=['Mat data/shared cov/cov iter',num2str(cIter) '.jpg'];
%         print ('-djpeg',fileName);
%         %imwrite(im2,['Mat data/shared cov/cov iter',num2str(cIter) '.jpg'],'jpg');
% 
%     end
    %=======================================================================
    %============================= E-STEP ==================================
    %=======================================================================
    
    %SIMON'S COMMENTS
    %For the i'th group of 4 training images we will have 1 hidden variable
    %hi which can take values 1 or 2 and we will have four hidden variables
    %aij which can take 3 values.
    
    %E-STEP
    totalLikelihood = 0;
    dataLikelihood=[];
    prHidden=[];
    
    %for each group of 4 people
    for(cIndiv=1:size(dataIDMatrix,2)) %training person
        indivIndx=find(dataIDMatrix(:,cIndiv));
        indivSamples=data(:,indivIndx);
        logProbAijAndHiGivenXi=[];
        for(cSamples=1:size(indivSamples,2)) %sample of training person
            for(cInd=1:N_DIST)%library person
                for(cApp=1:N_APP) %appearance of library person
                    %to get term1
                    logProbAijGivenHi = log(tDist(cInd,cApp).prior/sum([tDist(cInd,:).prior]));
                    logProbTermMatrix(cInd,cApp,cSamples) = getLogTProb(indivSamples(:,cSamples),tDist(cInd,cApp).mean,tDist(cInd,cApp).cov,tDist(cInd,cApp).dof);
                    logProbTermMatrix(cInd,cApp,cSamples) = logProbTermMatrix(cInd,cApp,cSamples)+logProbAijGivenHi;
                end
                %rescale
                logProbTermMatrix(cInd,:,cSamples) = logProbTermMatrix(cInd,:,cSamples)-max(logProbTermMatrix(cInd,:,cSamples));
                %convert to linear probs
                linMatrix = exp(logProbTermMatrix(cInd,:,cSamples));
                %to avoind log of zero, replace the zeros wit realMin
                linMatrix(linMatrix==0)=realmin;
                %normalize
                linMatrix = linMatrix/sum(linMatrix(:));
                %store back in log prob term matrix
                logProbTermMatrix(cInd,:,cSamples) = log(linMatrix);
            end
        end
        
        %to get term2
        %for each sample
        %collect the 2x3 probability table of pr(a_ij,h_i)
        for(cInd=1:N_DIST)%library person
            logLikeDataGivenHi(cInd) = 0;
            for(cSamples=1:size(indivSamples,2))
                
                for(cApp=1:N_APP) %appearance of library person
                    logProbAijAndHi = log(tDist(cInd,cApp).prior/sum([tDist(cInd,:).prior]));
                    %to get term1
                    logProbTerm2Matrix(cInd,cApp,cSamples) = getLogTProb(indivSamples(:,cSamples),tDist(cInd,cApp).mean,tDist(cInd,cApp).cov,tDist(cInd,cApp).dof);
                    logProbTerm2Matrix(cInd,cApp,cSamples) = logProbTerm2Matrix(cInd,cApp,cSamples)+logProbAijAndHi;
                end
                %summation of logs equals product over samples
                logLikeDataGivenHi(cInd) = logLikeDataGivenHi(cInd)+sumLog(logProbTerm2Matrix(cInd,:,cSamples));
            end;
        end
        %rescale
        logLikeDataGivenHi = logLikeDataGivenHi-max(logLikeDataGivenHi);
        %convert to linear
        likeDataGivenHi = exp(logLikeDataGivenHi);
        %to avoind log of zero, replace the zeros wit realMin
        likeDataGivenHi(likeDataGivenHi==0)=realmin;
        %normalize 
        probHiGivenAllSamples = likeDataGivenHi/sum(likeDataGivenHi);
        %convert back to logs
        logProbHiGivenAllSamples = log(probHiGivenAllSamples);
%         
        %log post probability
        for(cSamples = 1:size(indivSamples,2))
            for (cInd = 1:N_DIST)
                for(cApp = 1:N_APP)
                    logProbAijAndHiGivenXi(cInd,cApp,cSamples) = logProbTermMatrix(cInd,cApp,cSamples)+logProbHiGivenAllSamples(cInd);


                end;
            end;

        end;
        
        
        %convert back to linear
        linProbAijAndHiGivenXi=exp(logProbAijAndHiGivenXi);
        prHidden(:,:,cIndiv*4-3:cIndiv*4)=linProbAijAndHiGivenXi;
        
        %find data linear likelihood
       %for each library person
       for(cInd=1:N_DIST)
           %for each library appearance
           for(cApp=1:N_APP)
               %for each sample
               for (cSample=1:size(indivSamples,2))
                   %find the likelihood of each sample given hidden variables
                   dataLikeMatrix(cSample,cInd,cApp)=getLogTProb(indivSamples(:,cSamples),tDist(cInd,cApp).mean,tDist(cInd,cApp).cov,tDist(cInd,cApp).dof)...
                                    +log(tDist(cInd,cApp).prior);
               end
           end
       end
       

       for(cInd=1:N_DIST)
           likelihoodGivenHi=1;
           sampleLikelihood=[];
           %for each sample
           for (cSample=1:size(indivSamples,2))
               likelihoodGivenHi=likelihoodGivenHi+sumLog(dataLikeMatrix(cSample,cInd,:));
               %likelihood of each data sample to be used in check means
               sampleLikelihood(cSample)=likelihoodGivenHi;
           end
           totalLikelihood=totalLikelihood+likelihoodGivenHi;
       end
       
       dataLikelihood=[dataLikelihood sampleLikelihood];
    end
        
%        %find likelihood of each sample 
%        for(cSamples = 1:size(indivSamples,2))
%             sampleLogLike(cSamples)=sumLog(sampleLogLikeMatrix(cSamples,:,:));
%        end
%        %storing the likelihood of each datapoint in a matrix later to be
%        %used in check means
%        dataLogLikeMatrix=[dataLogLikeMatrix sampleLogLike];
%        
%        %find total likelihood of the data
%        totalDataLogeLikelihood=totalDataLogLikelihood+sum(sampleLogLike);
%        
       
    
    
    
    
        
     
    %for each group of 4 people we should wind up with 1 Pr(hi) (2x1)
    %vector
    %and 8 Pr(aij) vectors - 1 for each person = 4   x2 for the two possible values of hi 
           
   
%     %quit iterating if we are not improving
%     %[totalLogDataLikelihood prevLogLike totalLogDataLikelihood-prevLogLike]
%     if (((totalLikelihood-prevLogLike)<0.001&cIter>MIN_ITER)|cIter>MAX_ITER)
%         fprintf('\n%d Iterations, Likelihood %f\n',cIter,totalLikelihood');
%         break;
%     end;
%     prevLogLike = totalLikelihood;
%     %fprintf('Total Likelihood = %f\n',totalLikelihood);
%     
  
    
%     %======================================================================
%     %find log likelihood of the data after the M step
%     %calculate the likelihood of this group of samples
%     for(cSamples = 1:size(indivSamples,2))
%         for (cInd = 1:N_DIST)
%             for(cApp = 1:N_APP)
%                 sampleLogLikeMatrix(cSamples,cInd,cApp)=getLogLikelihoodGivenHidden(indivSamples(:,cSamples),tDist(cInd,cApp).mean,tDist(cInd,cApp).cov)...
%                     +log(tDist(cInd,cApp).prior);
%             end
%         end
%     end
% 
%     %find likelihood of each sample
%     for(cSamples = 1:size(indivSamples,2))
%         sampleLogLike(cSamples)=sumLog(sampleLogLikeMatrix(cSamples,:,:));
%     end
%     %storing the likelihood of each datapoint in a matrix later to be
%     %used in check means
%     dataLogLikeMatrix=[dataLogLikeMatrix sampleLogLike];
% 
%     %find total likelihood of the data
%     totalDataLogeLikelihood=totalDataLogLikelihood+sum(sampleLogLike);
%        
%     fprintf('Total Likelihood after M= %f\n',totalDataLogLikelihood)   
    
%     %===========================================================================    
%     %check that we haven't dissapeared into a local extrema thingy
%     for (c1 = 1:N_DIST)
%         for(c2=1:N_APP)
%             if (DIAG_FLAG)
%                 minCov = min(tDist(c1,c2).cov);
%             else
%                 [U L V] =svd(tDist(c1,c2).cov);
%                 minCov = min(diag(L));
%             end;
%             if minCov<minAcceptableCov
%                 failedFlag = 1;
%                 break;
%             end
%         end
%     end;
%     %if the process has failed then break out of the loop and return
%     %nothing
%     if (failedFlag)
%         totalDataLogLikelihood = -realmax;
%         tDist = [];
%         break;
%     end;
end %while    

 
prData=totalLikelihood;

%===================================================================

function logProb= getLogMVarGaussProb(data,mean,covar)

nDim = size(data,1);
covar(covar==0)=0.001;

if (size(covar,2)==1)%if cov is diagonal
    logProb = -0.5*nDim*log(2*pi)-0.5*sum(covar)-0.5*sum(((data-mean).^2)./covar,1);
else %if it is general (slower)
    logProb= -0.5*nDim*log(2*pi)-0.5*log(det(covar))-0.5*(data-mean)'*inv(covar)*(data-mean);
end;
    
%===================================================================
function r= drawGaussian(m,s,c,priority)
hold on;
angleInc = 0.1;

for (cAngle = 0:angleInc:2*pi)
    angle1 = cAngle;
    angle2 = cAngle+angleInc;
    [x1 y1] = getGaussian2SD(m,s,angle1);
    [x2 y2] = getGaussian2SD(m,s,angle2);
    h = patch([m(1) x1 x2],[m(2) y1 y2],[priority priority priority],c);
    set(h,'EdgeColor',c);
    plot([x1 x2],[y1 y2],'k-','LineWidth',2);
end

    
%===================================================================

    
    
    
function r= drawGaussianNoBorder(m,s,c,priority)
hold on;
angleInc = 0.1;

for (cAngle = 0:angleInc:2*pi)
    angle1 = cAngle;
    angle2 = cAngle+angleInc;
    [x1 y1] = getGaussian2SD(m,s,angle1);
    [x2 y2] = getGaussian2SD(m,s,angle2);
    h = patch([m(1) x1 x2],[m(2) y1 y2],[priority priority priority],c);
    set(h,'EdgeColor',c);
%    plot([x1 x2],[y1 y2],'k-','LineWidth',2);
end


%===================================================================

%draw 2DGaussian
function r= drawGaussianOutline(m,s,c)

hold on;
angleInc = 0.1;

for (cAngle = 0:angleInc:2*pi)
    angle1 = cAngle;
    angle2 = cAngle+angleInc;
    [x1 y1] = getGaussian2SD(m,s,angle1);
    [x2 y2] = getGaussian2SD(m,s,angle2);
    plot([x1 x2],[y1 y2],'k-','LineWidth',2,'Color',c);
end
%===================================================================

%find position of in xy co-ordinates at 2SD out for a certain angle
function [x,y]= getGaussian2SD(m,s,angle1)

if (size(s,2)==1)
    s = diag(s);
end;

vec = [cos(angle1) sin(angle1)];
factor = 4/(vec*inv(s)*vec');

x = cos(angle1) *sqrt(factor);
y = sin(angle1) *sqrt(factor);

x = x+m(1);
y = y+m(2);

%===================================================================

%return whether insied 2D Gaussian
function r= insideGauss2SD(m,s,x)

r= (((x-m)*inv(s)*(x-m)')<4);


% 
% ===================================================================
% function r=getLogLikelihoodGivenHidden(samples,mean,covar)
% 
% sampleDim=size(samples,2);
% prSamples=zeros(sampleDim,1);
% 
% for(i=1:sampleDim)
%     prSamples(i)=getLogMVarGaussProb(samples(:,i),mean,covar);
% end
% 
% logPrSamples=sum(prSamples);
% 
% r=logPrSamples;

% %=====================================================================
% function r=getProbHidden(samples,mean,covar)
% 
% sampleDim=size(samples,2);
% prSamples=zeros(sampleDim,1);
% 
% for(i=1:sampleDim)
%     prSamples(i)=exp(getLogProbHidden(samples(:,i),mean,covar));
% end
% 
% expPrSamples=sum(prSamples);
% 
% r=expPrSamples;
% 
%========================================================================
function r=checkMeanDist(means,data,likelihood);
% a function to check whether an two gaussians have the same mean or very
% close means. to avoid local minima
 

M_DIM=size(means,2);
N_DATA=size(data,2);

k1=kron(1:M_DIM,ones(1,M_DIM));
k2=kron(ones(1,M_DIM),1:M_DIM);

%find the distance between the means of the gaussians
dist=(means(:,k1)-means(:,k2)).^2;
dist=sqrt(sum(dist,1));

%remove the distance between a mean and itself
dist=reshape(dist,[M_DIM,M_DIM]);
djj=diag([ones(1,M_DIM)*100]);
dist=dist+djj;

%check if there is a very small distance between means
[x,y]=find(dist<0.1);



%determing rand order point
randPointOrder = randperm(N_DATA);

if(size(x,1)==0)
    r=means;
else
    %find the datapoint with the lowest likelihood
    l=likelihood;
    for(i=1:100)
    mn=find(l==min(l));
    lowLikelihoodIndx(i)=mn(1);
    l(lowLikelihoodIndx(i))=1000;
    end

    %avoid jumping to the same mean again
    %by taking a random index amongthe least likeliy data points  
    indx=round(rand*100);
    if(indx==0)
        indx=1;
    end
    dataIndx=lowLikelihoodIndx(indx);
    if(dataIndx>size(data,2))
        dataIndx
    end
    for(i=1:size(x,1))
        %checking if the mean has been randomized befor
        if(x(i)~=0)
            %setting the mean of the gaussian to a random point
            %means(:,x(i))= data(:,randPointOrder(i));
            means(:,x(i))= data(:,dataIndx);
            %removing the overlapping gaussian from the list
            x(find(x==y(i)))=0;
        end
    end


end
r=means;
%========================================================================
function r=getLikelihoodGivenHidden(samples,mean,covar)

sampleDim=size(samples,2);
prSamples=zeros(sampleDim,1);

for(i=1:sampleDim)
    prSamples(i)=exp(getLogProbHidden(samples(:,i),mean,covar));
end

expPrSamples=sum(prSamples);

r=expPrSamples;
%===================================================================
function r=getLogProbHidden(samples,mean,covar)

sampleDim=size(samples,2);
prSamples=zeros(sampleDim,1);
covar(covar==0)=0.1;

for(i=1:sampleDim)
    prSamples(i)=getLogMVarGaussProb(samples(:,i),mean,covar);
end

logPrSamples=sum(prSamples);

r=logPrSamples;

%========================================================================
function rt=sumLog(lst)
% a function that takes a list of terms, rescales the logs and return the
dim=size(lst,2);
logSumTerms=0;
sumTerms=0;
%find the rescale factor
rescaleFactor=max(lst(:));
%subtract the rescale factor
scaledTerm=lst-rescaleFactor;
%takes the exp
scaledTerm=exp(scaledTerm);
%add them together
sumTerms=sum(scaledTerm(:));
%convert back to log  probabilities and add the rescale factor
logSumTerms=log(sumTerms)+rescaleFactor;
 
rt=logSumTerms;

%===================================================================
%===================================================================

function logProb= getLogTProb(data,meanVal,covar,dof)

 [N_DIM N_DATA] = size(data);
 xMinusMean =data-repmat(meanVal,1,N_DATA);
 
 
 if (size(covar,2)==1)%if cov is diagonal
     deltaVal = sum(xMinusMean.*(repmat(1./covar,1,N_DATA)).*xMinusMean,1);
     detCov = prod(covar);
 else %if it is general (slower)
     deltaVal = sum((xMinusMean'*inv(covar))'.*xMinusMean,1);
     detCov = det(covar);
 end;
 
 if(detCov==inf)
     detCov=realmax;
 end

 
 if(detCov==0)
     detCov=realmin;
 end
 
 logProb = log(gamma(dof/2+N_DIM/2))-log(gamma(dof/2));
 logProb = logProb-0.5*log(detCov);
 logProb = logProb-(N_DIM/2)*log(dof*pi);
 logProb = logProb-(dof/2+N_DIM/2)*log(1+deltaVal/dof);
 
 %===================================================================
 %===================================================================
function [clusters,meanVec] = kMeans(data,N_CLUST)



[nDim nData] = size(data);
if (nData<N_CLUST)
    fprintf('Less datapoints than required clusters\n');
    return;
end

if (N_CLUST==1)
    clusters = ones(1,nData);
    meanVec = mean(data,2); 
    return;
end;



%if conditions are suitable, initialize using principal components
if (N_CLUST<nDim*2)
    dataMean = mean(data,2);
    dataMinusMean = data-repmat(dataMean,1,size(data,2));
    %Initialize means by placing at +/- 2Sd of each 
    %principal component in turn;
    meanVec = zeros(nDim,N_CLUST);
    thisF =  trainPCA(dataMinusMean);
    thisF = thisF(:,1:ceil(N_CLUST/2));
    L = (thisF'*dataMinusMean);
    L =L*L'/size(dataMinusMean,2);
    LSqrt = diag(sqrt(diag(L)));
    thisF = thisF*LSqrt;
    
    for (cClust = 1:N_CLUST)
        thisComp = ceil(cClust/2);
        if (rem(cClust,2))
            meanVec(:,cClust)=dataMean+2*thisF(:,thisComp); 
        else
            meanVec(:,cClust)=dataMean-2*thisF(:,thisComp);         
        end;
    end;
else
    %else initialize with random datapoints
    randOrder = randperm(nData);
    for (cClust = 1:N_CLUST)
        meanVec(:,cClust) = data(:,randOrder(cClust));
    end;
end;


N_ITER =6;
clusters = zeros(1,nData);
distToMeans = zeros(N_CLUST,nData);
cIter = 1;
while(cIter<=N_ITER)
    fprintf('K Means - Iteration %d\n',cIter);
    %classify each point as belonging to one of the means
    for (cClust =1:N_CLUST)
        distToMeans(cClust,:) = sum((data-repmat(meanVec(:,cClust),1,nData)).^2,1);
    end;
    
    [minDist clusters] = min(distToMeans);
    clusters = clusters(1,:);
    
    %check at least one point in each cluster
    for (cClust = 1:N_CLUST)
        N_PTS_PER_CLUST(cClust) = sum(clusters==cClust);
    end;
    if (min(N_PTS_PER_CLUST)==0)
        randOrder = randperm(nData);
        for (cClust = 1:N_CLUST)
            meanVec(:,cClust) = data(:,randOrder(cClust));
        end;
        cIter = 1;
        fprintf('Local Minimum - Re-initializing\n');
        continue;
  end;
    
    %update the means accordingly
    for (cClust =1:N_CLUST)
        meanVec(:,cClust) = mean(data(:,find(clusters==cClust)),2);
    end;
       cIter=cIter+1;
 
end;


%=======================================================
function [princComp,meanVec] = trainPCA(data)
%assumes that data is in columns


[nDim nData] = size(data);
meanVec= mean(data,2);
data = data-repmat(meanVec,1,nData);


if (nDim>nData)
    XXT = data'*data;
    [dummy LSq V]= svd(XXT);
    LInv = 1./sqrt(diag(LSq));
    princComp  = data * V * diag(LInv);
else
    [princComp L V] = svd(data*data');
end;

function r=optCrit(logNu,expectedLogU, expectedU,prHidden)


prHidden = length(prHidden)*prHidden/sum(prHidden);
totalWeight = sum(prHidden);
expectedLogU = (expectedLogU.*prHidden);
expectedU = (expectedU.*prHidden);

nu = exp(logNu);
N = length(expectedU);
r = N*nu*log(nu/2)/2-N*log(gamma(nu/2))+(nu/2-1)*sum(expectedLogU)-nu*sum(expectedU)/(2);
r = r*-1;
