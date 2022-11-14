function [V,Vc,Ji_now] = NI_function(o,p,i)
% compute Nikaido-Isoda function
% W. Ananduta
% 10/11/2021
    
    Ji_now = cost_compute(o,p,i);
    
    % compute aggregative terms
    o.sigma_mg = p.en.sumPd(1:p.h);
    o.sigma_gu = p.gn.sumGd(1:p.h);
    for j=1:p.n
        o.sigma_mg = o.sigma_mg + o.p_mg{j};
        o.sigma_gu = o.sigma_gu + o.d_gu{j};
    end
    
    
    V = zeros(p.n,1);
    Vc = zeros(p.n,1);
    
    
    % compute best cost for MI problem
    Ji_best = loc_opt_NI_MI(o,p,i);

    V = Ji_now - Ji_best;

    % compute best cost for convexified problem
    Ji_best_c = loc_opt_NI_c(o,p,i);
    Vc = Ji_now - Ji_best_c;
    
end