function [p,out,q] = alg2StageSOCP2(p)
% Iterative 2-stage method with SOCP model [W. Ananduta & S. Grammatico, 22]
% on 33b-20n network
% centralized steps
% 20/09/2022

    %p.gn.r = p.r;
    disp(['Proposed algorithm with SOCP model'])
    r_max = 5;
    gamma(1) = 0;
    gammaUpper = inf;
    gammaLower = gamma(1);
    p.pen = 0;
    p.Gamma_pen = 0;
    %%
    
    for r = 1:r_max
        tic
        disp(['iteration=',num2str(r)])
        
        p.Gamma_pen = gamma(r);
        %p.pen = gamma(r);

        %% Stage 1
        p.solverFlag = 1; %1 = GUROBI, 2= fmincon
        p.scpFlag = 0;
        p.fixedInt_flag = 0; 
        [p,o]= solveCentralized_SOCP(p);
        disp('Stage 1 done')
        
        %% Stage 2
        o = solve_binarySOCP(o,p);

        [o] = solve_pressureSOCP(o,p);

        disp('Stage 2 done')

        o.pen_w(r) = p.pen;

        er(r) = o.Jpsi;


        % compute cost
        for i=1:p.n
            [q.J(i,r),q.P(i,r)] = cost_compute(o,p,i);
        end

        q.Jt(r) = sum(q.J(:,r));
        q.Pt(r) = sum(q.P(:,r));

        q.gamma = gamma;
        q.Jpsi = er;
        
        if r == 1
            out = o;
        else
            if er(r) < er(r-1)
                out = o;
            end
        end
        
        q.gammaLower = gammaLower;
        q.gammaUpper = gammaUpper;

        % compute gas-flow error
        q.er_gf(r) = gasFlow_error(p,o);

        %% Update gamma
        [gamma(r+1),gammaLower,gammaUpper]=gammaRules(er(r),gamma(r),gammaLower,gammaUpper,r);
        q.time(r) = toc;

    end
    
        
end
