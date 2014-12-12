% Function to calculate BIC from a mixture structure and data
% BIC = gmmem_BIC(gmix, data)

function BIC = gmmem_BIC(gmix, data)

nr_clust = gmix.ncentres;
covar_type = gmix.covar_type;
[nr_dp, nr_dim] = size(data);
if strcmp(covar_type, 'full'),

    %     nr_parm = nr_clust + nr_clust*nr_dim + nr_clust*sum([1:1:nr_dim]);
    nr_parm = gmix.nwts;

else

    nr_parm = gmix.nwts;

end
logLik = -sum(log(gmmprob(gmix, data)));
BIC = 2*logLik + nr_parm*log(nr_dp);
