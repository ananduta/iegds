function [p,o,q] = algSOCP_ce(p)
% with convex envelope of concave constraints    


    
    for r = 1:1
        tic
        p.solverFlag = 2; %1 = GUROBI, 2= fmincon
        p.scpFlag = 0;
        [p,o]= solveCentralized_SOCP_ce(p);
        disp('SOCP method with CE done')
        % compute cost
        
        for i=1:p.n
            [q.J(i),q.P(i)] = cost_compute(o,p,i);
        end

        q.Jt = sum(q.J);
        q.Pt = sum(q.P);

        % compute gas-flow error
        q.er_gf = gasFlow_error(p,o);
        q.time = toc;
    end
    
end