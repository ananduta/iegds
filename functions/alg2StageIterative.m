function [p,out,q] = alg2StageIterative(p)
% Iterative 2-stage method [W. Ananduta & S. Grammatico, 22]
% on 33b-20n network
% centralized steps
% 20/09/2022

    p.gn.r = p.r;
    
    r_max = 10;
    gamma(1) = 0;
    gammaUpper = inf;
    gammaLower = gamma(1);
    p.Gamma_pen = 0; %for SOCP
    p.solver = 2; % 2= OSQP, 1= Gurobi
    %%
    
    for r = 1:r_max
        tic
        disp(['iteration=',num2str(r)])
        
        p.pen = gamma(r);


        %% Stage 1
        [p,o]= solveCentralized_GNEPc(p);
        disp('Stage 1 done')
        %% Stage 2
        o = solve_binary(o,p);

        [o,~] = solve_pressure3(o,p);

        disp('Stage 2 done')

        o.pen_w(r) = p.pen;

        er(r) = o.gfv_max;


        % compute cost
        for i=1:p.n
            [q.J(i,r),q.P(i,r)] = cost_compute(o,p,i);
        end

        q.Jt(r) = sum(q.J(:,r));
        q.Pt(r) = sum(q.P(:,r));

        q.gamma = gamma;
        q.Jpsi = er;
        
        
        q.gammaLower = gammaLower;
        q.gammaUpper = gammaUpper;

        % compute gas-flow error
        q.er_gf(r) = gasFlow_error(p,o);
        
        if r == 1
            out = o;
        else
            if q.er_gf(r).max < q.er_gf(r-1).max
                out = o;
            end
        end
        %% Update gamma
        [gamma(r+1),gammaLower,gammaUpper]=gammaRules(er(r),gamma(r),gammaLower,gammaUpper,r);
        q.time(r) = toc;

    end
    
        
end
