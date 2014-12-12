% [gmix_opt, BIC_c, gmix_c] = netlab_wrapper(data, Cmin, Cmax, vec_it,
% covar_type, verbose, varargin)
function [gmix_opt, BIC_c, gmix_c] = netlab_wrapper(data, Cmin, Cmax, ...
    vec_it, covar_type, verbose, varargin)

nr_c = length(vec_it);
c_vec = Cmin:Cmax;
[nr_dp, nr_dim] = size(data);
if covar_type == 1,
    
    covar_type = 'diag';
    
elseif covar_type == 2,
    
    covar_type = 'full';
    
elseif covar_type == 3,
    
    covar_type = 'ppca';
    ppca_dim = varargin{1};
    
end

opt_std = foptions;
if verbose,
    
    opt_std(1) = 1;
    
else
    
    opt_std(1) = -1;
    
end
opt_std(3) = 1e-3;
opt_std(5) = 1;
opt_std(14) = 1e3;
for it_c = 1:nr_c,

    opt = opt_std;
    try

        if verbose,
            
            fprintf('kmeans initialization\n')
            
        end
        if strcmp(covar_type, 'full'),

            gmix_init = gmm(nr_dim, c_vec(it_c), 'full');
            nr_parm = c_vec(it_c) + c_vec(it_c)*nr_dim + c_vec(it_c)*sum([1:1:nr_dim]);

        elseif strcmp(covar_type,  'ppca'),

            gmix_init = gmm(nr_dim, c_vec(it_c), 'ppca', ppca_dim);
            nr_parm = gmix_init.nwts;

        else

            gmix_init = gmm(nr_dim, c_vec(it_c), covar_type);
            nr_parm = gmix_init.nwts;

        end
        gmix_init = gmminit(gmix_init, data, opt);
        [gmix_it{1}, opt] = gmmem(gmix_init, data, opt);
        logLik = -opt(8);
        gmix_it{1}.logLik = logLik;
        BIC_it(1) = -2*logLik + nr_parm*log(nr_dp);

    catch

        BIC_it(1) = Inf;

    end

    mult_fact = ceil(vec_it(it_c)/nr_dp);
    rand_vec = [];
    for it_rand = 1:mult_fact,

        rand_vec_temp = zeros(c_vec(it_c), nr_dp);
        for cc = 1:c_vec(it_c),

            rand_vec_temp(cc, :)  = randperm(nr_dp);

        end
        rand_vec = [rand_vec rand_vec_temp];

    end

    cnt_it = [1:1:vec_it(it_c)];
    for it = 2:vec_it(it_c),

        if strcmp(covar_type, 'full'),

            gmix_init = gmm(nr_dim, c_vec(it_c), 'full');
            nr_parm = c_vec(it_c) + c_vec(it_c)*nr_dim + c_vec(it_c)*sum([1:1:nr_dim]);

        elseif strcmp(covar_type, 'ppca'),

            gmix_init = gmm(nr_dim, c_vec(it_c), 'ppca', ppca_dim);
            nr_parm = gmix_init.nwts;

        else

            gmix_init = gmm(nr_dim, c_vec(it_c), covar_type);
            nr_parm = gmix_init.nwts;

        end
        gmix_init.centres = data(rand_vec(:, it), :);

        opt = opt_std;
        try

            [gmix_it{it}, opt] = gmmem(gmix_init, data, opt);
            logLik = -opt(8);
            gmix_it{it}.logLik = logLik;
            BIC_it(it) = -2*logLik + nr_parm*log(nr_dp);

        catch

            if strcmp(covar_type, 'full'),

                gmix_init = gmm(nr_dim, c_vec(it_c), 'full');
                nr_parm = c_vec(it_c) + c_vec(it_c)*nr_dim + c_vec(it_c)*sum([1:1:nr_dim]);

            elseif strcmp(covar_type, 'ppca'),

                gmix_init = gmm(nr_dim, c_vec(it_c), 'ppca', ppca_dim);
                nr_parm = gmix_init.nwts;

            else

                gmix_init = gmm(nr_dim, c_vec(it_c), covar_type);
                nr_parm = gmix_init.nwts;

            end
            gmix_init.centres = data(rand_vec(:, it), :);
            [gmix_it{it}, opt] = gmmem(gmix_init, data, opt);
            logLik = -opt(8);
            gmix_it{it}.logLik = logLik;
            BIC_it(it) = -2*logLik + nr_parm*log(nr_dp);

        end
        if verbose
            if find(cnt_it == it),
                
                fprintf('Iteration %d in %d clusters\n', it, c_vec(it_c));
                
            end
        end

    end
    it_opt = find(BIC_it == min(BIC_it));
    if verbose
        if it_opt == 1,
            
            fprintf('Fuzzy clustering won...\n');
            
        end
    end
    BIC_c(it_c) = BIC_it(it_opt(1));
    gmix_c{it_c} = gmix_it{it_opt(1)};
    gmix_c{it_c}.nr_parm = nr_parm;
    if it_opt == 1,

        gmix_c{it_c}.flg_fcm = 1;


    else

        gmix_c{it_c}.flg_fcm = 0;

    end

end

c_opt = find(BIC_c == min(BIC_c));
gmix_opt = gmix_c{c_opt(1)};

% -------------------------------------------------------------------------------
% -------------------------------------------------------------------------------

function gmix = fcm_init(gmix, data)

GMM_WIDTH = 1.0;

[nr_dp, nr_dim] = size(data);
parm_fcm.e = 1e-9;
parm_fcm.m = 7;
parm_fcm.c = gmix.ncentres;
data_fcm.X = data;

results_fcm = GKclust(data_fcm, parm_fcm);
centres_fcm = results_fcm.cluster.v;
post_fcm = results_fcm.data.f;
gmix.centres = centres_fcm;
for k = 1:parm_fcm.c,

    gmix.priors(k) = mean(post_fcm(:, k));
    % Pick out data points belonging to this centre
    data_c = data(post_fcm(:, k) == max(post_fcm, [], 2), :);
    diffs = data_c - (ones(size(data_c, 1), 1) * gmix.centres(k, :));
    gmix.covars(:,:,k) = (diffs'*diffs)/(size(data_c, 1));
    % Add GMM_WIDTH*Identity to rank-deficient covariance matrices
    if rank(gmix.covars(:,:,k)) < gmix.nin,

        gmix.covars(:,:,k) = gmix.covars(:,:,k) + GMM_WIDTH.*eye(gmix.nin);

    end

end
