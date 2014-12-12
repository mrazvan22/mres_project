function [mix, options, post, errlog] = gmmem2(mix, x, options)
%GMMEM	EM algorithm for Gaussian mixture model.
%
%	Description
%	[MIX, OPTIONS, ERRLOG] = GMMEM(MIX, X, OPTIONS) uses the Expectation
%	Maximization algorithm of Dempster et al. to estimate the parameters
%	of a Gaussian mixture model defined by a data structure MIX. The
%	matrix X represents the data whose expectation is maximized, with
%	each row corresponding to a vector.    The optional parameters have
%	the following interpretations.
%
%	OPTIONS(1) is set to 1 to display error values; also logs error
%	values in the return argument ERRLOG. If OPTIONS(1) is set to 0, then
%	only warning messages are displayed.  If OPTIONS(1) is -1, then
%	nothing is displayed.
%
%	OPTIONS(3) is a measure of the absolute precision required of the
%	error function at the solution. If the change in log likelihood
%	between two steps of the EM algorithm is less than this value, then
%	the function terminates.
%
%	OPTIONS(5) is set to 1 if a covariance matrix is reset to its
%	original value when any of its singular values are too small (less
%	than MIN_COVAR which has the value eps).   With the default value of
%	0 no action is taken.
%
%	OPTIONS(14) is the maximum number of iterations; default 100.
%
%	The optional return value OPTIONS contains the final error value
%	(i.e. data log likelihood) in OPTIONS(8).
%
%	See also
%	GMM, GMMINIT
%

%	Copyright (c) Ian T Nabney (1996-2001)

% Check that inputs are consistent

[ndata, xdim, nr_sub] = size(x);

% Sort out the options
if (options(14))
  niters = options(14);
else
  niters = 100;
end

display = options(1);
store = 0;
if (nargout > 2)
  store = 1;	% Store the error values to return them
  errlog = zeros(1, niters);
end
test = 0;
if options(3) > 0.0
  test = 1;	% Test log likelihood for termination
end

check_covars = 0;
if options(5) >= 1
    if display >= 0
        %     disp('check_covars is on');
    end
    check_covars = 1;	% Ensure that covariances don't collapse
    MIN_COVAR = 2*eps;	% Minimum singular value of covariance matrix
    init_covars = mix{1}.covars;
end

% Main loop of algorithm
for n = 1:niters

    % Calculate posteriors based on old parameters

    post = zeros(ndata, mix{1}.ncentres, nr_sub);
    for sub = 1:nr_sub,
               
        [post(:, :, sub), act(:, :, sub)] = ...
            gmmpost(mix{sub}, x(:, :, sub));
        
    end
    post = mean(post, 3);
    act = mean(act, 3);
    % Calculate error value if needed
    if (display | store | test)
        
        for sub = 1:nr_sub,
            
            prob = act*(mix{sub}.priors)';
            % Error value is negative log likelihood of data
            e(sub) = - sum(log(prob));
            
        end
        e = sum(e);
        if store
            errlog(n) = e;
        end
        if display > 0
            fprintf(1, 'Cycle %4d  Error %11.6f\n', n, e);
        end
        if test
            if (n > 1 & abs(e - eold) < options(3))
                options(8) = e;
                return;
            else
                eold = e;
            end
        end
    end

    % Adjust the new estimates for the parameters
    new_pr = sum(post, 1);

    for sub = 1:nr_sub,        
            
        new_c = post' * x(:, :, sub);

        % Now move new estimates to old parameter vectors
        mix{sub}.priors = new_pr ./ ndata;
        mix{sub}.centres = new_c ./ (new_pr' * ones(1, mix{sub}.nin));

        switch mix{sub}.covar_type
            case 'spherical'
                n2 = dist2(x(:, :, sub), mix{sub}.centres);
                for j = 1:mix.ncentres
                    v(j) = (post(:,j)'*n2(:,j));
                end
                mix{sub}.covars = ((v./new_pr))./mix{sub}.nin;
                if check_covars
                    % Ensure that no covariance is too small
                    for j = 1:mix{sub}.ncentres
                        if mix{sub}.covars(j) < MIN_COVAR
                            mix{sub}.covars(j) = init_covars(j);
                        end
                    end
                end
            case 'diag'
                for j = 1:mix{sub}.ncentres
                    diffs = x(:, :, sub) - (ones(ndata, 1) * mix{sub}.centres(j,:));
                    mix{sub}.covars(j,:) = sum((diffs.*diffs).*(post(:,j)*ones(1, ...
                        mix{sub}.nin)), 1)./new_pr(j);
                end
                if check_covars
                    % Ensure that no covariance is too small
                    for j = 1:mix{sub}.ncentres
                        if min(mix{sub}.covars(j,:)) < MIN_COVAR
                            mix{sub}.covars(j,:) = init_covars(j,:);
                        end
                    end
                end
            case 'full'
                for j = 1:mix{sub}.ncentres
                    diffs = x(:, :, sub) - (ones(ndata, 1) * mix{sub}.centres(j,:));
                    diffs = diffs.*(sqrt(post(:,j))*ones(1, mix{sub}.nin));
                    mix{sub}.covars(:,:,j) = (diffs'*diffs)/new_pr(j);
                end
                if check_covars
                    % Ensure that no covariance is too small
                    for j = 1:mix{sub}.ncentres
                        if min(svd(mix{sub}.covars(:,:,j))) < MIN_COVAR
                            mix{sub}.covars(:,:,j) = init_covars(:,:,j);
                        end
                    end
                end
            case 'ppca'
                for j = 1:mix{sub}.ncentres
                    diffs = x(:, :, sub) - (ones(ndata, 1) * mix{sub}.centres(j,:));
                    diffs = diffs.*(sqrt(post(:,j))*ones(1, mix{sub}.nin));
                    [tempcovars, tempU, templambda] = ...
                        ppca((diffs'*diffs)/new_pr(j), mix{sub}.ppca_dim);
                    if length(templambda) ~= mix{sub}.ppca_dim
                        error('Unable to extract enough components');
                    else
                        mix{sub}.covars(j) = tempcovars;
                        mix{sub}.U(:, :, j) = tempU;
                        mix{sub}.lambda(j, :) = templambda;
                    end
                end
                if check_covars
                    if mix{sub}.covars(j) < MIN_COVAR
                        mix{sub}.covars(j) = init_covars(j);
                    end
                end
            otherwise
                error(['Unknown covariance type ', mix{sub}.covar_type]);
        end
    end
end

for sub = 1:nr_sub,
    
    opt(sub) = -sum(log(gmmprob(mix{sub}, x(:, :, sub))));
    
end
options(8) = sum(opt);
if (display >= 0)
    disp(maxitmess);
end
