function er = gasFlow_error(p,o)
% Compute error of gas flow non-linear equations
% W. Ananduta
% 01/03/2022
c = 1;
for i = 1:p.n
    for jj = 1:length(p.gn.N{i})
        j = p.gn.N{i}(jj);
        
        for h = 1:p.h
        % gas flow equation
        
            if o.psi{i} >= o.psi{j}
                sig = -1;
            else
                sig = 1;
            end

            if abs(o.phi{i,j}(h)) < 1e-6
                phi = 0;
            else
                phi = o.phi{i,j}(h);
            end


            flow(c,h) = sqrt(abs(o.psi{i}(h)-o.psi{j}(h)));

            if flow < 1e-6
                flow(c,h) = 0;
                gf(c,h) = 0;
                gf_abs(c,h) = 0;
                g_f(i,j,h) = gf(c,h);
            else
                gf(c,h) = ((1/sqrt(p.gn.cf(i,j))*(phi)) - sig*sqrt(abs(o.psi{i}(h)-o.psi{j}(h))))/flow(c,h);
                gf_abs(c,h) = (1/sqrt(p.gn.cf(i,j))*(phi)) - sig*sqrt(abs(o.psi{i}(h)-o.psi{j}(h)));
                g_f(i,j,h) = gf(c,h);
            end

            if abs(gf(c,h)) < 1e-5
                gf_tight(c,h) = 1;
            else
                gf_tight(c,h) = 0;
            end
        
        
        end
        c = c+1;
    end
    
    
end
er.max = max(max(abs(gf)));
er.mean = mean(mean(abs(gf)));
er.gf = gf;
er.gf_abs = gf_abs;
er.gf_tight = gf_tight;
er.flow = flow;
end


