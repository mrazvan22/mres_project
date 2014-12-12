function [likelihood, gmix_struct]=computeMixTLikelihood(data_controls,data_patients,threshold_flag,min_prior_flag)

data=[data_controls,data_patients];
[nr_events, nr_subj]=size(data);
%prior_control=size(data_controls,2)/nr_subj;
prior_control=0.25;

%cov type parameter
%covType = 0;  %unconstrained - individual covariance matrices
%covType = 1;  %unconstrained - shared covariance matrices
%covType = 2;  %diagonal - individual covariance matrices
%covType = 3;  %diagonal - shared covariance matrices
covType = 4;  %spherical - individual covariance matrices
%covType = 5;  %spherical - shared covariance matrices


for roi=1:nr_events
    figure
    tEst_control = mixTFit(data_controls(roi,:),1,covType);
    figure
    tEst_tot = mixTFitFixedCompVersion2(data(roi,:),2,covType,tEst_control(1).mean,prior_control,min_prior_flag);
    tEst_patient=tEst_tot(1);
    tEst_patient.prior=1;
    
    gmix_struct(roi).tEst_tot=tEst_tot;
    gmix_struct(roi).tEst_control=tEst_control;
    gmix_struct(roi).tEst_patient=tEst_patient;

        
    %% find likelihood event_occurred and event_not_occured
    likelihood(roi, :, 1) = exp(getLogTProb(data_patients(roi,:),tEst_control.mean,tEst_control.cov,tEst_control.dof));
    likelihood(roi, :, 2) = exp(getLogTProb(data_patients(roi,:),tEst_patient.mean,tEst_patient.cov,tEst_patient.dof));
%  
%     %% Visualize how well the T-distributions fit the data
%  
%     p_controls=exp(getLogTProb(data_controls(roi,:),tEst_control.mean,tEst_control.cov,tEst_control.dof));
%     p_patients=exp(getLogTProb(data_patients(roi,:),tEst_patient.mean,tEst_patient.cov,tEst_patient.dof));
%     p_all_MoT=getExpMixTProb(data(roi,:),tEst_tot);
%     
%     [hist_c, x_c] = ksdensity(data_controls(roi,:));
%     [hist_p, x_p] = ksdensity(data_patients(roi,:));
%     
%     
%     figure
%     subplot(121), hold on
%     plot(x_c, hist_c, 'g');
%     plot(x_p, hist_p, 'r');
%     
%     subplot(122), hold on
%     plot(data_controls(roi,:), p_controls, 'g.')
%     plot(data_patients(roi,:), p_patients, 'r.');
%     plot(data(roi,:), p_all_MoT, 'b.');
%     
%     set(gcf,'Color',[1 1 1])    
%     pause


   
if threshold_flag
    % set a maximumn likelilhood ratio
    T=[3.2];
    ratio=abs(log(likelihood(roi, :, 1))-log(likelihood(roi, :, 2)));
    indx=find(ratio>log(T));
    
    for cIndx=1:length(indx)
        thisLike=likelihood(roi,indx(cIndx),:);
        %if you are in the far left of the distributrion
        if data_patients(roi,indx(cIndx))<tEst_control.mean
            minIndx=find(thisLike==min(thisLike));
            mxIndx=find(thisLike==max(thisLike));
            thisLike(minIndx)=thisLike(mxIndx)/T;
        else
            % if you are in the far right of the distribution
            minIndx=find(thisLike==min(thisLike));
            mxIndx=find(thisLike==max(thisLike));
            thisLike(mxIndx)=thisLike(minIndx)*T;
        end
        likelihood(roi,indx(cIndx),:)=thisLike;
    end
end

indx=find(data_patients(roi,:)>tEst_control.mean);
like2=likelihood(roi,:,2);
like2(indx)=0;
likelihood(roi,:,2)=like2;
    
end