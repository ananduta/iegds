function [p,out,q] = algSOCP_SCP(p)


% Define algorithm parameters
p.Gamma_pen = 0; % no penalty/cost on auxiliary var. Gamma
p.scpFlag = 1;   % 1=include auxiliary constraint for SCP method
p.tauMax = 1e5;  % maximum weight of auxiliary penalty variable
p.pen = 0.1;       % initial weight of auxiliary penalty variable
p.mu = 5;        % increment for the weight of auxiliary penalty variable
p.epsS = 1e-3; p.epsC = 1e-2;     % stopping criterion
 p.solverFlag = 2; %1 = GUROBI, 2= fmincon
 
% Initialize flow
for i = 1:p.n
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        p.flow{i,j} = 0*zeros(p.h,1);
    end
end
costValNow = 1e3;
%sumPenNow = 1e3;
tic
k = 1;
iterMax = 30;
while 1
    costValOld = costValNow;
    
   
    
    [p,o]= solveCentralized_SOCP(p);
    disp('SOCP method done')
    
    % compute cost

    for i=1:p.n
        [q.J(i,k),q.P(i,k)] = cost_compute(o,p,i);
    end

    q.Jt(k) = sum(q.J(:,k));
    q.Pt(k) = sum(q.P(:,k));
    
    costValNow = o.costVal;
    
    
    % compute gas-flow error
    q.er_gf(k) = gasFlow_error(p,o);
    disp(['gas flow error:'])
    q.er_gf(k)
    
    if k == 1
        out = o;
    else
        if q.er_gf(k).max <= q.er_gf(k-1).max
            out = o;
        end
    end
    
    % update penalty tau
    p.pen = min(p.pen*p.mu,p.tauMax);
    
    % update flow decision and check slack variables
    sumPen = 0;
    for i = 1:p.n
        for jj = 1:p.gn.noN(i)
            j = p.gn.N{i}(jj);
            p.flow{i,j} = o.phi{i,j};
            sumPen = sumPen + sum(o.pen{i,j});
        end
    end
    
    disp(['sum of slack vars:',num2str(sumPen)])
    % stopping condition
    disp(['difference consecutive cost value:',num2str(norm(costValNow - costValOld))])
    if norm(costValNow - costValOld) <= p.epsC && sumPen <= p.epsS
        break
    end
    if k > iterMax
        break
    end
    k = k+1;
end
q.time = toc;
end