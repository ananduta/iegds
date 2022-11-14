function o = solve_binarySOCP(o,p)
% Obtain binary variables for SOCP model of gas-flow equations
% W. Ananduta
% 27/09/2022
eps = 0*1e-7;
    
    for i = 1:p.n
        for jj = 1:p.gn.noN(i)
            j=p.gn.N{i}(jj);
                       
            for h = 1:p.h
                
                % flow direction alpha_{i,j}
                if o.phi{i,j}(h,1) <= 0
                    o.alpha_s{i,j}(h,1) = 0;
                else
                    o.alpha_s{i,j}(h,1) = 1;
                end
            
            end
            
        end
        
    end
    
end