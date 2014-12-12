function mix = gmminit2(mix, x, options)
%GMMINIT Initialises Gaussian mixture model from data
%
%	Description
%	MIX = GMMINIT(MIX, X, OPTIONS) uses a dataset X to initialise the
%	parameters of a Gaussian mixture model defined by the data structure
%	MIX.  The k-means algorithm is used to determine the centres. The
%	priors are computed from the proportion of examples belonging to each
%	cluster. The covariance matrices are calculated as the sample
%	covariance of the points associated with (i.e. closest to) the
%	corresponding centres. For a mixture of PPCA model, the PPCA
%	decomposition is calculated for the points closest to a given centre.
%	This initialisation can be used as the starting point for training
%	the model using the EM algorithm.
%
%	See also
%	GMM
%

%	Copyright (c) Ian T Nabney (1996-2001)

[ndata, xdim, nr_sub] = size(x);

% Arbitrary width used if variance collapses to zero: make it 'large' so
% that centre is responsible for a reasonable number of points.
GMM_WIDTH = 1.0;

% Use kmeans algorithm to set centres
options(5) = 1;	
[mix{1}.centres, options, post_mix{1}] = ...
    kmeans_netlab(mix{1}.centres, x(:, :, 1), options);
for sub = 2:nr_sub,
    
    [mix{sub}.centres, options, post_mix{sub}] = ...
        kmeans_netlab(mix{sub}.centres, x(:, :, sub), ...
        options, post_mix{sub-1});
    
end

% Set priors depending on number of points in each cluster

for sub = 1:nr_sub,
    
    cluster_sizes = max(sum(post_mix{sub}, 1), 1);  % Make sure that no prior is zero
    mix{sub}.priors = cluster_sizes/sum(cluster_sizes); % Normalise priors
    
end

for sub = 1:nr_sub,
    
    post = post_mix{sub};
    switch mix{sub}.covar_type

        case 'spherical'

            if mix{sub}.ncentres > 1
                % Determine widths as distance to nearest centre
                % (or a constant if this is zero)
                cdist = dist2(mix{sub}.centres, mix{sub}.centres);
                cdist = cdist + diag(ones(mix{sub}.ncentres, 1)*realmax);
                mix{sub}.covars = min(cdist);
                mix{sub}.covars = mix{sub}.covars + GMM_WIDTH*(mix{sub}.covars < eps);
            else
                % Just use variance of all data points averaged over all
                % dimensions
                mix{sub}.covars = mean(diag(cov(x(:, :, sub))));
            end
        case 'diag'
            for j = 1:mix{sub}.ncentres
                % Pick out data points belonging to this centre
                c = x(find(post(:, j)), :, sub);
                diffs = c - (ones(size(c, 1), 1) * mix{sub}.centres(j, :));
                mix{sub}.covars(j, :) = sum((diffs.*diffs), 1)/size(c, 1);
                % Replace small entries by GMM_WIDTH value
                mix{sub}.covars(j, :) = mix{sub}.covars(j, :) + GMM_WIDTH.*(mix{sub}.covars(j, :)<eps);
            end
        case 'full'
            for j = 1:mix{sub}.ncentres
                % Pick out data points belonging to this centre
                c = x(find(post(:, j)),:, sub);
                diffs = c - (ones(size(c, 1), 1) * mix{sub}.centres(j, :));
                mix{sub}.covars(:,:,j) = (diffs'*diffs)/(size(c, 1));
                % Add GMM_WIDTH*Identity to rank-deficient covariance matrices
                if rank(mix{sub}.covars(:,:,j)) < mix{sub}.nin
                    mix{sub}.covars(:,:,j) = mix{sub}.covars(:,:,j) + GMM_WIDTH.*eye(mix{sub}.nin);
                end
            end
        case 'ppca'
            for j = 1:mix{sub}.ncentres
                % Pick out data points belonging to this centre
                c = x(find(post(:,j)),:, sub);
                diffs = c - (ones(size(c, 1), 1) * mix{sub}.centres(j, :));
                [tempcovars, tempU, templambda] = ...
                    ppca((diffs'*diffs)/size(c, 1), mix{sub}.ppca_dim);
                if length(templambda) ~= mix{sub}.ppca_dim
                    error('Unable to extract enough components');
                else
                    mix{sub}.covars(j) = tempcovars;
                    mix{sub}.U(:, :, j) = tempU;
                    mix{sub}.lambda(j, :) = templambda;
                end
            end
        otherwise
            error(['Unknown covariance type ', mix{sub}.covar_type]);
    end
    
end

