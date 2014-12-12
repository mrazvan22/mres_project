function tDist = mixTFitFixedComp(data,N_DIST,covType,fix_comp_mean,control_prior,min_prior_flag)
%fits a Gaussian mixture model with N_DIST or fewer components
%
%covType = 0;  %unconstrained - individual covariance matrices
%covType = 1;  %unconstrained - shared covariance matrices
%covType = 2;  %diagonal - individual covariance matrices
%covType = 3;  %diagonal - shared covariance matrices
%covType = 4;  %spherical - individual covariance matrices
%covType = 5;  %spherical - shared covariance matrix

rand('seed',1)
randn('seed',1)

if (covType<2)
    DIAG_FLAG = 0;
else
    DIAG_FLAG = 1;
end

%measure covariance of data.
%find smallest tolerable SD for fitted gauss - sd must be at last 10000 of
%the Sd in the smallest direction

[N_DIM, N_DATA] = size(data);
prDist = 1/N_DIST;
 
%initialize responsibilities
[clusters,meanVec] = kMeans(data,N_DIST);
prHidden = zeros(N_DIST,N_DATA);
for (cDist = 1:N_DIST)
    prHidden(cDist,find(clusters==cDist)) =1;
end;

%initialize U values
a = 1.6/2;b = 1.6/2;
expectedU = ones(N_DIST,N_DATA);
expectedLogU = (psi(a)-log(b))*ones(N_DIST,N_DATA);

MAX_ITER = 20;
MIN_ITER = 20;
DRAW_FLAG = 1;

if (DRAW_FLAG&N_DIM==2)
    figure;
    randColors = 0.3+0.7*rand(N_DIST,3);
end;

%initialize prior
alpha0 = 1;



cIter = 0;
while(1)
    tic
    cIter = cIter+1;
    if (cIter>MAX_ITER) 
        break;
    end;
    fprintf('%d...',cIter);

    %=======================================================================
    %============================= M-STEP ==================================
    %=======================================================================
    

    
    %calculate prior
    N = sum(prHidden,2)';   
    prior = (alpha0+N)/sum(alpha0+N);  %dirichlet prior on weights  
    for (cDist =1:N_DIST)
        tDist(cDist).prior = prior(cDist);    
    end;
   
    %force a minimum prior for the unaccfected component equal to control_prior
    if min_prior_flag
        if tDist(1).prior<control_prior
            prior_adjust_value=(tDist(1).prior-control_prior)/(N_DIST-1);
            tDist(1).prior=control_prior;
            
            for cDist=2:N_DIST
                tDist(cDist).prior=tDist(cDist).prior+prior_adjust_value;
            end
        end
    end
    %keep the first mean as teh fixed comp mean
    tDist(1).mean=fix_comp_mean;
    %update the remaining means
    for (cDist = 2:N_DIST)
        tDist(cDist).mean = sum(repmat(expectedU(cDist,:).*prHidden(cDist,:),N_DIM,1).*data,2);
        tDist(cDist).mean = tDist(cDist).mean/sum(expectedU(cDist,:).*prHidden(cDist,:),2);
    end;
    
    %calculate covariance 
    if (covType==0|covType==2|covType==4)%    - if individual
        for (cDist = 1:N_DIST)
            xMinusMean = data-repmat(tDist(cDist).mean,1,N_DATA);
            if (covType==0) %full covariance
                tDist(cDist).cov = zeros(N_DIM);
                for(cData = 1:N_DATA)
                    tDist(cDist).cov = tDist(cDist).cov+prHidden(cDist,cData)*expectedU(cDist,cData)*xMinusMean(:,cData)*xMinusMean(:,cData)';
                end;
            elseif(covType==2) % diagonal covariance
                tDist(cDist).cov = sum((repmat(prHidden(cDist,:).*expectedU(cDist,:),N_DIM,1).*xMinusMean).*xMinusMean,2);
            else %(covType==4) % spherical covariance
                tDist(cDist).cov = sum(sum((repmat(prHidden(cDist,:).*expectedU(cDist,:),N_DIM,1).*xMinusMean).*xMinusMean,2));
                tDist(cDist).cov = tDist(cDist).cov*ones(N_DIM,1)/N_DIM;
            end;
            tDist(cDist).cov = tDist(cDist).cov/sum(prHidden(cDist,:));
        end;
        
    else%pooled covariance
        if (covType==1)
            tDist(1).cov = zeros(N_DIM);
        else
            tDist(1).cov = zeros(N_DIM,1);
        end;
        
        for (cDist = 1:N_DIST)
            xMinusMean = data-repmat(tDist(cDist).mean,1,N_DATA);
            if (covType==1) %full covariance
                for(cData = 1:N_DATA)
                    tDist(1).cov = tDist(1).cov+prHidden(cDist,cData)*expectedU(cDist,cData)*xMinusMean(:,cData)*xMinusMean(:,cData)';
                end;
            elseif(covType==3) % diagonal covariance
                tDist(1).cov = tDist(1).cov+sum((repmat(prHidden(cDist,:).*expectedU(cDist,:),N_DIM,1).*xMinusMean).*xMinusMean,2);
            else %(covType==5) % spherical covariance
                tDist(1).cov = tDist(1).cov+sum(sum((repmat(prHidden(cDist,:).*expectedU(cDist,:),N_DIM,1).*xMinusMean).*xMinusMean,2))/N_DIM;
            end;
        end;
        %normalize
        tDist(1).cov = tDist(1).cov/N_DATA;
        
        %copy to remaining points
        for (cDist = 2:N_DIST)
            tDist(cDist).cov = tDist(1).cov;
        end;
    end;
    
    %update degrees of freedom
    for (cDist = 1:N_DIST)
        estNu= fminbnd(@(x) optCrit(x,expectedLogU(cDist,:),expectedU(cDist,:),prHidden(cDist,:)),-1,5);
        tDist(cDist).dof=exp(estNu);
       %      tDist(cDist).dof=1.6;
  
 
    end;
    
    %=======================================================================
    %========================== DRAW DATA ==================================
    %=======================================================================
    
    if (DRAW_FLAG&N_DIM==2)
        minX = min(data(1,:)); maxX = max(data(1,:));
        minY = min(data(2,:)); maxY = max(data(2,:));
        X = minX:(maxX-minX)/99:maxX;
        Y = minY:(maxY-minY)/99:maxY;
        nX = length(X);nY = length(Y);
        X = repmat(X,[nY 1]);
        Y = repmat(Y',[1 nX]);
        X2 = X(:);Y2 = Y(:);
        plot(data(1,:),data(2,:),'k.');hold on;

        for (cDist = 1:N_DIST)
            Z = exp(log(tDist(cDist).prior)+getLogTProb([X2';Y2'],tDist(cDist).mean,tDist(cDist).cov,tDist(cDist).dof));
            Z = reshape(Z,size(X));
            contour(X,Y,Z); colormap hsv;hold on;

        end;
        axis square; axis equal; hold off; drawnow;

              
        %plot(data(1,:),data(2,:),'k.');
        %hold on;
        %for (cDist = 1:N_DIST)
        %    set(gcf,'Color',[1 1 1]);
        %    %drawGaussian(tDist(cDist).mean,tDist(cDist).cov,randColors(cDist,:),0.01+0.01*cDist);
        %    drawGaussianOutline(tDist(cDist).mean,tDist(cDist).cov,randColors(cDist,:));
        %    hold on;
        %end;
        %axis square;axis equal;
        %axis off;
        %hold off;drawnow;
    end;
    if (DRAW_FLAG&N_DIM==1)
        [h x] = hist(data,N_DATA/4);
        xWidth = x(2)-x(1);
        bar(x,h/(xWidth*N_DATA),1);hold on;
        x2 = min(x):(max(x)-min(x))/1000:max(x);
        prob=zeros(size(x2));
        for (cDist = 1:N_DIST);
            estMu = tDist(cDist).mean;
            estSigma = tDist(cDist).cov;
            estNu = tDist(cDist).dof;
            estPrior = tDist(cDist).prior;
            prob = prob+estPrior*tpdf((x2-estMu)/sqrt(estSigma),estNu)/sqrt(estSigma);
        end;
        plot(x2,prob,'r-','LineWidth',2);set(gca,'Box','off');
        hold off;drawnow;
    end;
    
    
    %=======================================================================
    %============================= E-STEP ==================================
    %=======================================================================
    
    %calculate distribution over hidden variables - do in log probs to
    %avoid numerical problems.
    totalLogDataLikelihood = 0;
    for (cDist = 1:N_DIST)
        prHidden(cDist,:) = log(tDist(cDist).prior)+getLogTProb(data,tDist(cDist).mean,tDist(cDist).cov,tDist(cDist).dof);
    end;
    rescaleFactor = max(prHidden,[],1);
    prHidden = prHidden-repmat(rescaleFactor,N_DIST,1);
    %normalize to create posterior probability for each point
    prHidden = exp(prHidden);
    prHidden = prHidden./repmat(sum(prHidden,1),N_DIST,1);
    
    %calcualate expected U values
    a = zeros(N_DIST,N_DATA);
    b = zeros(N_DIST,N_DATA);
    for (cDist = 1:N_DIST)
        a(cDist,:) = repmat(tDist(cDist).dof/2+N_DIM/2,1,N_DATA);
        xMinusMean =data-repmat(tDist(cDist).mean,1,N_DATA);
        if (covType<2) % if full covariance
            deltaVal = sum((xMinusMean'*inv(tDist(cDist).cov))'.*xMinusMean,1);
        else
            deltaVal = sum(xMinusMean.*(repmat(1./tDist(cDist).cov,1,N_DATA)).*xMinusMean,1);
        end;
        b(cDist,:) = tDist(cDist).dof/2+0.5*deltaVal;
    end;
    expectedU = a./b;
    expectedLogU = psi(a)-log(b);
    toc
end %while    




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
 
 logProb = log(gamma(dof/2+N_DIM/2))-log(gamma(dof/2));
 logProb = logProb-0.5*log(detCov);
 logProb = logProb-(N_DIM/2)*log(dof*pi);
 logProb = logProb-(dof/2+N_DIM/2)*log(1+deltaVal/dof);
 

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
