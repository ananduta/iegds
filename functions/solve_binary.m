function o = solve_binary(o,p)
% Obtain binary variables for PWA model of gas-flow equations
% W. Ananduta
% 06/10/2021
eps = 1e-8;
    
    for i = 1:p.n
        for jj = 1:p.gn.noN(i)
            j=p.gn.N{i}(jj);
            
            % define the parameters of the affine function at each region
            nf = @(y) y^2/p.gn.cf(i,j);
            pwaf = pwa_approx_nf(p.r,-p.gn.phi_max(i,j),p.gn.phi_max(i,j),nf);
            
            for h = 1:p.h
                
                % \delta_(i,j)^\psi
                if o.phi{i,j}(h,1) <= 0
                    o.delta_psi_s{i,j}(h,1) = 1;
                else
                    o.delta_psi_s{i,j}(h,1) = 0;
                end
            
            
                % \delta_(i,j)^m            
            
                for m = 1:p.r

                    if o.phi{i,j}(h,1) >= pwaf.m(m)-eps && o.phi{i,j}(h,1) <= pwaf.M(m)+eps
                        o.delta_s{i,j}{m}(h,1) = 1;
                    else
                        o.delta_s{i,j}{m}(h,1) = 0;
                    end
                end
            end
            
        end
        
    end
    
end