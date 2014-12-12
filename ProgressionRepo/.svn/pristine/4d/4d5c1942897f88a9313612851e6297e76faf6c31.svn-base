% Wrapper function that fits a mixture of fixed and free Gaussian
% components to the data, using gmmem_fixcomp.m, which is an adjusted 
% version of netlab's gmmem.m 
% [gmix_fixcomp, BIC_fixcomp, gmix_struct] = ...
%     gmmem_fixcomp_wrapper(gmix_fix, data, nr_extracomp_max, nr_it, ...
%     opts, flag_direction_extracomponents)
% flag_direction_extracomponents indicates in which direction the extra
% components are placed
% flag_direction_extracomponents == 1: only extra components to the smaller 
% than existing comopnents are placed
% flag_direction_extracomponents == 2: only extra components to the larger 
% than existing comopnents are placed
% flag_direction_extracomponents == 3: extra components both smaller and
% larger than existing comopnents are placed

function [gmix_fixcomp, BIC_fixcomp, gmix_struct] = ...
    gmmem_fixcomp_wrapper(gmix_fix, data, nr_extracomp_max, nr_it, ...
    opts, flag_direction_extracomponents)

nr_centres_start = gmix_fix.ncentres;
BIC_fixcomp = zeros(nr_extracomp_max + 1, 1);
gmix_struct{1} = gmix_fix;
BIC_fixcomp(1) = gmmem_BIC(gmix_fix, data);
for nr_extracomp = 1:nr_extracomp_max,
    
    BIC_it = zeros(nr_it, 1);
    for it = 1:nr_it,
        
        gmix_it{it} = gmm(1, nr_centres_start + nr_extracomp, 'full');
        if flag_direction_extracomponents == 1,
            
            minplus = -1*ones(nr_extracomp, 1);
            
        elseif flag_direction_extracomponents == 2,
            
            minplus = ones(nr_extracomp, 1);
            
        elseif flag_direction_extracomponents == 3,
            
            minplus = rand(nr_extracomp, 1);
            minplus(minplus > 0.5) = 1;
            minplus(minplus < 0.5) = -1;
            
        end
        for extracomp = 1:nr_extracomp,
            
            if minplus == -1,
                
                gmix_it{it}.centres(nr_centres_start + extracomp) = ...
                    min(gmix_fix.centres) + (min(data) - min(gmix_fix.centres))*rand;
                
            elseif minplus == 1,
                
                gmix_it{it}.centres(nr_centres_start + extracomp) = ...
                    max(gmix_fix.centres) + (max(data) - max(gmix_fix.centres))*rand;
                
            end
            
        end
        try
            gmix_it{it} = gmmem_fixcomp(gmix_it{it}, gmix_fix, ...
                data, opts);
            BIC_it(it) = gmmem_BIC(gmix_it{it}, data);
        catch
            BIC_it(it) = Inf;           
        end
        
    end
    if length(find(isinf(BIC_it))) == nr_it,
        
        BIC_fixcomp(nr_extracomp + 1) = Inf;
        gmix_struct{nr_extracomp + 1} = [];
        
    else
        
        it_opt = find(BIC_it == min(BIC_it));
        gmix_struct{nr_extracomp + 1} = gmix_it{it_opt(1)};
        BIC_fixcomp(nr_extracomp + 1) = BIC_it(it_opt(1));
        
    end
    
end
it_opt = find(BIC_fixcomp == min(BIC_fixcomp));
gmix_fixcomp = gmix_struct{it_opt(1)};        