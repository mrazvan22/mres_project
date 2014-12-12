function [gmix_out, BIC] = gmmm(gmix_in, data1, data2)

MIN_COVAR = eps;
check_covars = 0;

gmix_out = gmix_in;

ndata = size(data1, 1);
post = gmmpost(gmix_in, data1);
% Adjust the new estimates for the parameters
new_pr = sum(post, 1);
new_c = post' * data2;

% Now move new estimates to old parameter vectors
nr_dim = size(data2, 2);
gmix_out.nin = nr_dim;
gmix_out.priors = new_pr ./ ndata;
gmix_out.centres = new_c ./ (new_pr' * ones(1, gmix_out.nin));
gmix_dummy = gmm(gmix_out.nin, gmix_out.ncentres, gmix_out.covar_type);
init_covars = gmix_dummy.covars;

switch gmix_in.covar_type
    case 'spherical'
        n2 = dist2(data2, gmix_out.centres);
        for j = 1:gmix_out.ncentres
            v(j) = (post(:,j)'*n2(:,j));
        end
        gmix_out.covars = ((v./new_pr))./gmix_out.nin;
        if check_covars
            % Ensure that no covariance is too small
            for j = 1:gmix_out.ncentres
                if gmix_out.covars(j) < MIN_COVAR
                    gmix_out.covars(j) = init_covars(j);
                end
            end
        end
    case 'diag'
        gmix_out.covars = zeros(gmix_out.ncentres, gmix_out.nin);
        for j = 1:gmix_out.ncentres
            diffs = data2 - (ones(ndata, 1) * gmix_out.centres(j,:));
            gmix_out.covars(j,:) = sum((diffs.*diffs).*(post(:,j)*ones(1, ...
                gmix_out.nin)), 1)./new_pr(j);
        end
        if check_covars
            % Ensure that no covariance is too small
            for j = 1:gmix_out.ncentres
                if min(gmix_out.covars(j,:)) < MIN_COVAR
                    gmix_out.covars(j,:) = init_covars(j,:);
                end
            end
        end
    case 'full'
        gmix_out.covars = zeros(gmix_out.nin, gmix_out.nin, gmix_out.ncentres);
        for j = 1:gmix_out.ncentres
            diffs = data2 - (ones(ndata, 1) * gmix_out.centres(j,:));
            diffs = diffs.*(sqrt(post(:,j))*ones(1, gmix_out.nin));
            gmix_out.covars(:,:,j) = (diffs'*diffs)/new_pr(j);
        end
        if check_covars
            % Ensure that no covariance is too small
            for j = 1:gmix_out.ncentres
                if min(svd(gmix_out.covars(:,:,j))) < MIN_COVAR
                    gmix_out.covars(:,:,j) = init_covars(:,:,j);
                    disp(j)
                end
            end
        end
    case 'ppca'
        gmix_out.U = zeros(gmix_out.nin, gmix_out.ppca_dim, gmix_out.ncentres); 
        for j = 1:gmix_out.ncentres,
            diffs = data2 - (ones(ndata, 1) * gmix_out.centres(j,:));
            diffs = diffs.*(sqrt(post(:,j))*ones(1, gmix_out.nin));
            [tempcovars, tempU, templambda] = ...
                ppca((diffs'*diffs)/new_pr(j), gmix_out.ppca_dim);
            if length(templambda) ~= gmix_out.ppca_dim
                error('Unable to extract enough components');
            else
                gmix_out.covars(j) = tempcovars;
                gmix_out.U(:, :, j) = tempU;
                gmix_out.lambda(j, :) = templambda;
            end
        end
        if check_covars
            if gmix_out.covars(j) < MIN_COVAR
                gmix_out.covars(j) = init_covars(j);
            end
        end
    otherwise
        error(['Unknown covariance type ', gmix_out.covar_type]);
end
prob = gmmprob(gmix_out, data2);
if strcmp(gmix_out.covar_type, 'full'),

    nr_parm = gmix_out.ncentres + gmix_out.ncentres*nr_dim + gmix_out.ncentres*sum([1:1:nr_dim]);

else

    nr_parm = gmix_out.nwts;

end
logLik = sum(log(prob));
BIC = -2*logLik + nr_parm*log(ndata);

