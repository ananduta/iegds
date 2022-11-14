%% check feasibility of solution to the second stage\
c = 1;
for i = 1:p.n
    for jj = 1:p.gn.noN(i)
        j= p.gn.N{i}(jj);
        
        for h = 1:p.h
            cond(c,h) = -(1-2*o.alpha_s{i,j}(h))*(o.psi{i}(h) - o.psi{j}(h)) -o.tau{i,j}(h) + (1/p.gn.cf(i,j))*(o.phi{i,j}(h)^2); 
            cond1(c,h) = -(1-2*o.alpha_s{i,j}(h))*(o.psi{i}(h) - o.psi{j}(h)) + (1/p.gn.cf(i,j))*(o.phi{i,j}(h)^2); 
        end
        c = c+1;
    end
end