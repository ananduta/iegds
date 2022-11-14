function [p,o,q] = algSOCP(p)
    


    
    for r = 1
        tic
        p.solverFlag = 2; %1 = GUROBI, 2= fmincon
        p.scpFlag = 0;
        [p,o]= solveCentralized_SOCP(p);
        disp('SOCP method done')
        % compute cost
        
        for i=1:p.n
            [q.J(i,r),q.P(i,r)] = cost_compute(o,p,i);
        end

        q.Jt(r) = sum(q.J(:,r));
        q.Pt(r) = sum(q.P(:,r));

        % compute gas-flow error
        q.er_gf(r) = gasFlow_error(p,o);
        q.time(r) = toc;
    end
    
end