function [bestGauss,bestPrData,label]  = mixGaussFit_fixedcomp(data,fix_mean,N_GAUSS_IN,covType,nAttempts);
%fits a Gaussian mixture model with N_GAUSS or fewer components
%
%covType = 0;  %unconstrained - individual covariance matrices
%covType = 1;  %unconstrained - shared covariance matrices
%covType = 2;  %diagonal - individual covariance matrices
%covType = 3;  %diagonal - shared covariance matrices
%covType = 4;  %spherical - individual covariance matrices
%covType = 5;  %spherical - shared covariance matrix

N_GAUSS= N_GAUSS_IN;

N_ITER = nAttempts;
failedFlag = 0;
bestGauss = []; bestPrData = -realMax;
while (1)
    %try fitting mixture of Gaussian N_ITER times and choose the best  If it keeps failing, reduce number of Gaussians 
    fprintf('Fitting %d Gaussians\n',N_GAUSS);
    for (cIter = 1:N_ITER)
        fprintf('ATTEMPT %d\n',cIter);
        [allGaussEst{cIter} allPrData(cIter),label]  = mixGaussFitInternal(data,fix_mean,N_GAUSS,covType);
        gaussEst=allGaussEst{cIter};
        prData=allPrData(cIter);
        if (prData==-realMax)
            fprintf('FITTING FAILED\n');
        end;
        if (prData>bestPrData)
            bestPrData = prData;
            bestGauss = gaussEst;
        end;   
    end;
    %if succeeded then return
    if ~isempty(bestGauss)
        break;
    end;
        
    %if failed all 10 times then reduce number of Gaussians by 1
    N_GAUSS = N_GAUSS-1;
    
end;

%if we have fit less Gaussians then we were asked to then add some dummy
%ones at the end with weight zero so the output has the desired
%dimensionality
for (cGauss = N_GAUSS+1:N_GAUSS_IN)
    bestGauss(cGauss).prior = 0;
    bestGauss(cGauss).mean = bestGauss(N_GAUSS).mean;
    bestGauss(cGauss).cov = bestGauss(N_GAUSS).cov;
end;

gaussEst = bestGauss;


function [gaussEst,prData,label]  = mixGaussFitInternal(data,fix_mean,N_GAUSS,covType);
%function gaussEst = mixGaussFit(data,N_GAUSS,covType);
%fits N_GAUSS gaussians to the passed data (1 pt per column)  
%
%covType = 0;  %unconstrained - individual covariance matrices
%covType = 1;  %unconstrained - shared covariance matrices
%covType = 2;  %diagonal - individual covariance matrices
%covType = 3;  %diagonal - shared covariance matrices
%covType = 4;  %spherical - individual covariance matrices
%covType = 5;  %spherical - shared covariance matrix


[N_DIM, N_DATA] = size(data);
prGauss = 1/N_GAUSS;
covData = cov(data');

%find smallest tolerable SD for fitted gauss - sd must be at last 1000 of
%the Sd in the smallest direction
minAcceptableCov = min(diag(covData))/10000;

%determing rand order point
randPointOrder = randperm(N_DATA);

%set up initial values for Gaussian
for (c1 = 1:N_GAUSS)
    gauss(c1).mean = data(:,randPointOrder(c1));
    gauss(c1).prior = prGauss;
    gauss(c1).cov = N_DIM*diag(diag(covData));
end;
gauss(1).mean=fix_mean;

MAX_ITER = 50;
MIN_ITER = 20;
DRAW_FLAG = 1;

if (DRAW_FLAG&N_DIM==2)
    figure;
    randColors = 0.3+0.7*rand(N_GAUSS,3);
end;

%previous log likelihood (initially v. small!)
prevLogLike = -realmax;
cIter = 0;
prHidden = zeros(N_GAUSS,N_DATA);
while(1)
    cIter = cIter+1;
    fprintf('%d...',cIter);
    %=======================================================================
    %========================== DRAW DATA ==================================
    %=======================================================================
    
    if (DRAW_FLAG&N_DIM==2)
        plot(data(1,:),data(2,:),'k.');
        for (cGauss = 1:N_GAUSS)
            set(gcf,'Color',[1 1 1]);
            %drawGaussian(gauss(cGauss).mean,gauss(cGauss).cov,randColors(cGauss,:),0.01+0.01*cGauss);
            drawGaussianOutline(gauss(cGauss).mean,gauss(cGauss).cov,randColors(cGauss,:));
            hold on;
        end;
        plot(data(1,:),data(2,:),'k.');
        axis square;axis equal;
        axis off;
        hold off;drawnow;
    end;
    
    %=======================================================================
    %============================= E-STEP ==================================
    %=======================================================================
    
    %calculate distribution over hidden variables - do in log probs to
    %avoid numerical problems.
    totalLogDataLikelihood = 0;
    for (cData = 1:N_DATA)
        if (covType>=2) % if diagonal then only pass diagonal which is fast to invert
            for (cGauss = 1:N_GAUSS)
                 prHidden(cGauss,cData) = log(gauss(cGauss).prior)+getLogMVarGaussProb(data(:,cData),gauss(cGauss).mean,diag(gauss(cGauss).cov));
            end;
        else %pass full covariance matrix - must do inversion, which is slower
            for (cGauss = 1:N_GAUSS)
                 prHidden(cGauss,cData) = log(gauss(cGauss).prior)+getLogMVarGaussProb(data(:,cData),gauss(cGauss).mean,gauss(cGauss).cov);
            end;            
        end;
        
        rescaleFactor = max(prHidden(:,cData));
        prHidden(:,cData) = prHidden(:,cData)-rescaleFactor;
        
        totalPr = sum(exp(prHidden(:,cData)));
        
        totalLogPr = log(totalPr)+rescaleFactor;
        totalLogDataLikelihood = totalLogDataLikelihood+totalLogPr;
        
        %normalize to create posterior probability for each point
        prHidden(:,cData) = exp(prHidden(:,cData));
        prHidden(:,cData) = prHidden(:,cData)/sum(prHidden(:,cData));
    end;
    %quit iterating if we are not improving
    %[totalLogDataLikelihood prevLogLike totalLogDataLikelihood-prevLogLike]
    if (((totalLogDataLikelihood-prevLogLike)<0.001&cIter>MIN_ITER)|cIter>MAX_ITER)
        fprintf('\n%d Iterations, Log Likelihood %f\n',cIter,totalLogDataLikelihood');
        break;
    end;
    prevLogLike = totalLogDataLikelihood;
    [~,label(1,:)] = max(prHidden,[],1);
    %=======================================================================
    %============================= M-STEP ==================================
    %=======================================================================
    
    failedFlag = 0;
    
    for (cGauss = 1:N_GAUSS)
        %check that this Gaussian does have some pts associated with it
        if (sum(prHidden(cGauss,:))==0)
            fprintf('FAILED - One gaussian with no weight assigned\n');
            failedFlag = 1;
            break;
        end;
        %calculate new prior
        gauss(cGauss).prior = sum(prHidden(cGauss,:))/sum(prHidden(:));
        %calculate new mean
        weightedMean = sum(repmat(prHidden(cGauss,:),N_DIM,1).*data,2);
        gauss(cGauss).mean = weightedMean/sum(prHidden(cGauss,:));
    end;
    gauss(1).mean=fix_mean;
    
    if failedFlag
        break;
    end;
    
    
    %calculate new covariance
    switch covType
        case 0 %covType = 0;  %unconstrained - individual covariance matrices
            for (cGauss = 1:N_GAUSS)
                weightedOuterProduct = zeros(N_DIM,N_DIM);
                for (cData = 1:N_DATA)
                    weightedOuterProduct =weightedOuterProduct+prHidden(cGauss,cData)*(data(:,cData)-gauss(cGauss).mean)*(data(:,cData)-gauss(cGauss).mean)';
                end;
                gauss(cGauss).cov = weightedOuterProduct/sum(prHidden(cGauss,:));
            end;
            
        case 1%covType = 1;  %unconstrained - shared covariance matrices
            
            weightedOuterProduct = zeros(N_DIM,N_DIM);
            for (cGauss = 1:N_GAUSS)
                for (cData = 1:N_DATA)
                    weightedOuterProduct =weightedOuterProduct+prHidden(cGauss,cData)*(data(:,cData)-gauss(cGauss).mean)*(data(:,cData)-gauss(cGauss).mean)';
                end;
            end;
            for (cGauss = 1:N_GAUSS)
                gauss(cGauss).cov = weightedOuterProduct/sum(prHidden(:));
            end;
            
        case 2%covType = 2;  %diagonal - individual covariance matrices
            for (cGauss = 1:N_GAUSS)
                gauss(cGauss).cov = zeros(N_DIM,N_DIM);
                for (cDim = 1:N_DIM)
                    gauss(cGauss).cov(cDim,cDim) = sum(prHidden(cGauss,:).*(data(cDim,:)-gauss(cGauss).mean(cDim)).*(data(cDim,:)-gauss(cGauss).mean(cDim)))/sum(prHidden(cGauss,:));    
                end;
            end;
        case 3%covType = 3;  %diagonal - shared covariance matrices
            gauss(1).cov = zeros(N_DIM,N_DIM);
            
            for (cGauss = 1:N_GAUSS)
                for (cDim = 1:N_DIM)
                    gauss(1).cov(cDim,cDim) = gauss(1).cov(cDim,cDim)+sum(prHidden(cGauss,:).*(data(cDim,:)-gauss(cGauss).mean(cDim)).*(data(cDim,:)-gauss(cGauss).mean(cDim)));    
                end;
            end;
            gauss(1).cov = gauss(1).cov/sum(prHidden(:));
            for (cGauss = 1:N_GAUSS)
                gauss(cGauss).cov = gauss(1).cov;
            end;
            
        case 4%covType = 4;  %spherical - individual covariance matrices
            for (cGauss = 1:N_GAUSS)
                weightedSqDiff = 0;
                for (cDim = 1:N_DIM)
                    weightedSqDiff = weightedSqDiff+prHidden(cGauss,:).*(data(cDim,:)-gauss(cGauss).mean(cDim)).*(data(cDim,:)-gauss(cGauss).mean(cDim));
                end;
                gauss(cGauss).cov = eye(N_DIM) *sum(weightedSqDiff,2)/(N_DIM*sum(prHidden(cGauss,:)));
            end;
            
        case 5%covType = 5;  %spherical - shared covariance matrix
            weightedSqDiff = 0;
            for (cGauss = 1:N_GAUSS)
                for (cDim = 1:N_DIM)
                    weightedSqDiff = weightedSqDiff+prHidden(cGauss,:).*(data(cDim,:)-gauss(cGauss).mean(cDim)).*(data(cDim,:)-gauss(cGauss).mean(cDim));
                end;
            end;
            for (cGauss = 1:N_GAUSS)
                gauss(cGauss).cov = eye(N_DIM) *sum(weightedSqDiff,2)/(N_DIM*sum(prHidden(:)));
            end;
            
    end;
    
    %===========================================================================    
    %check that we haven't dissapeared into a local extrema thingy
    for (cGauss=1:N_GAUSS)
        [U L V] =svd(gauss(cGauss).cov);
        minCov = min(diag(L));
        if minCov<minAcceptableCov
            failedFlag = 1;
            break;
        end
    end;
    %if the process has failed then break out of the loop and return
    %nothing
    if (failedFlag)
        totalLogDataLikelihood = -realmax;
        gauss = [];
        break;
    end;
end %while    


gaussEst = gauss;
prData=totalLogDataLikelihood;


%===================================================================

function logProb= getLogMVarGaussProb(data,mean,covar)

nDim = size(data,1);

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



