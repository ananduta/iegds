function [o] = solve_pressureSOCP(o,p)
% Compute pressure decisions with minimum violation of gas-flow SOCP model
% Centralized scheme
% W. Ananduta
% 27/09/2022

sumid = 0;
for i = 1:p.gn.n
    dim(i) = p.h*(1+p.gn.noN(i))+p.gn.noN(i); %[psi, [tau_i,j], [s^psi_i,j]]
    idStart(i) = sumid + 1;
    idEnd(i) = sumid + dim(i);
    sumid = idEnd(i);
end

Aall = [];
ball = [];
for i=1:p.gn.n 
        
        % bounds of psi
            A3 = zeros(1,1+p.gn.noN(i));
            A3(1) = 1;
            At3 = [kron(eye(p.h),A3) zeros(p.h,p.gn.noN(i)) ];
            bt3 = p.gn.psi_max(i)*ones(p.h,1);
            
            A4 = zeros(1,1+p.gn.noN(i));
            A4(1) = -1;
            At4 = [kron(eye(p.h),A4) zeros(p.h,p.gn.noN(i)) ];
            bt4 = -p.gn.psi_min(i)*ones(p.h,1);
            
            Ali = [At3;At4];
            bli = [bt3;bt4];
            
            Ai = zeros(size(Ali,1),idEnd(p.n));
            Ai(:,idStart(i):idEnd(i)) = Ali;
            Aall = [Aall;Ai];
            ball = [ball;bli];
            
        
        for jj = 1:p.gn.noN(i)
            j = p.gn.N{i}(jj);
    
            % local constraints
          
            % non-negative tau_ij
            A1 =  zeros(1,1+p.gn.noN(i));
            A1(1+jj) = -1;
            
            At1 = [kron(eye(p.h),A1) zeros(p.h,p.gn.noN(i)) ];
            bt1 = zeros(p.h,1);
            
            
            Alj{j} = [At1];
            blj{j} = [bt1];
                        
        end
        if p.gn.noN(i)>0
            Al{i} = cat(1,Alj{:});
            bl{i} = cat(1,blj{:});

            % non-negative s_ij^psi
            At2 = [zeros(p.gn.noN(i),p.h*(1+p.gn.noN(i))) -eye(p.gn.noN(i)) ];
            bt2 = zeros(p.gn.noN(i),1);
            Al{i} = [Al{i}; At2];
            bl{i} = [bl{i}; bt2];
            
        else
            Al{i} = [];
            bl{i} = [];
        end
        clearvars('Alj','blj');
        Ai = zeros(size(Al{i},1),idEnd(p.n));
        Ai(:,idStart(i):idEnd(i)) = Al{i};
        Aall = [Aall; Ai];
        ball = [ball; bl{i}];
end

% coupling constraints
for i=1:p.gn.n 
        id_psi = 1;
        for jj = 1:p.gn.noN(i)
            j = p.gn.N{i}(jj);
    
            id_tau = 1+jj;
            
            
            for h = 1:p.h
                % (1-2alpha)(psi_i - psi_j) - 1 s_ij^psi <= phi_ij^2/cf_ij
                E1i{h} = zeros(1,1+p.gn.noN(i));
                E1i{h}(id_psi) = 1-2*o.alpha_s{i,j}(h);
                                
                E1j{h} = zeros(1,1+p.gn.noN(j));
                E1j{h}(id_psi) = -(1-2*o.alpha_s{i,j}(h));

                b1(h,1) = o.phi{i,j}(h)^2/p.gn.cf(i,j);
                
                % -(1-2alpha)(psi_i - psi_j) - 1 s_ij^psi <= -phi_ij^2/cf_ij
                E3i{h} = zeros(1,1+p.gn.noN(i));
                E3i{h}(id_psi) = -(1-2*o.alpha_s{i,j}(h));
                                
                E3j{h} = zeros(1,1+p.gn.noN(j));
                E3j{h}(id_psi) = (1-2*o.alpha_s{i,j}(h));

                b3(h,1) = -(o.phi{i,j}(h)^2/p.gn.cf(i,j));
                
                % -(1-2alpha)(psi_i - psi_j) - tau_ij <= -phi_ij^2/cf_ij
                E2i{h} = zeros(1,1+p.gn.noN(i));
                E2i{h}(id_psi) = -(1-2*o.alpha_s{i,j}(h));
                E2i{h}(id_tau) = -1;
                
                E2j{h} = zeros(1,1+p.gn.noN(j));
                E2j{h}(id_psi) = (1-2*o.alpha_s{i,j}(h));

                b2(h,1) = -(o.phi{i,j}(h)^2/p.gn.cf(i,j));
                
            end
            % coefficients of s_ij^psi
            E1ia = zeros(1,p.gn.noN(i));
            E1ia(jj) = -1;
            
            E3ia = zeros(1,p.gn.noN(i));
            E3ia(jj) = -1;
            
            E2ia = zeros(1,p.gn.noN(i));
            
            % concatenation of the matrices
            E1it = [blkdiag(E1i{:}) kron(ones(p.h,1),E1ia) ];
            
            E3it = [blkdiag(E3i{:}) kron(ones(p.h,1),E3ia) ];
            
            E2it = [blkdiag(E2i{:}) kron(ones(p.h,1),E2ia) ];
            
            E{i,j}{1} = [E1it;E2it;E3it];
            
            E{i,j}{2} = [blkdiag(E1j{:}) zeros(p.h,p.gn.noN(j));
                         blkdiag(E2j{:}) zeros(p.h,p.gn.noN(j));
                         blkdiag(E3j{:}) zeros(p.h,p.gn.noN(j))];
            
            Ac = zeros(3*p.h,idEnd(p.n));
            Ac(:,idStart(i):idEnd(i)) = E{i,j}{1};
            Ac(:,idStart(j):idEnd(j)) = E{i,j}{2};
            
            bc = [b1;b2;b3];
            
            Aall = [Aall; Ac];
            ball = [ball; bc];
            
            clearvars('E1j','E2j','E1i','E2i');
        end
    
end

% Solve with OSQP


Aineq = Aall;
bineq = ball;

% create matrix for cost
for i = 1:p.gn.n
    fi = [0; ones(p.gn.noN(i),1)];
    fit{i} = [kron(ones(p.h,1),fi); ones(p.gn.noN(i),1)];
end

f = cat(1,fit{:});
H = zeros(length(f));

A = sparse(Aineq);
b = sparse(bineq);

% Create an OSQP object
prob = osqp;

settings = prob.default_settings();
settings.eps_abs= 1e-10;
settings.eps_rel= 1e-10;
settings.eps_prim_inf = 1e-10;
settings.eps_dual_inf = 1e-10;
settings.max_iter = 1e7;
settings.verbose = 0;
% Setup workspace and change alpha parameter
prob.setup(H, f, A, [], b, settings);
res = prob.solve();

u = res.x;
o.Jpsi = res.info.obj_val;

% [x,fval,exitflag,output] = quadprog(H, f, A, b);
% o.psi_solver = x;


% % Gurobi
% model.A = sparse(Aall);
% model.rhs = ball;
% %J = x'*Hloc*x + f'*x;
% model.Q = sparse(H);
% model.obj = f;
% 
% 
% % lower bound
% model.lb = -inf*ones(idEnd(p.n),1);
% 
% params.outputflag = 1;
%         params.BarHomogeneous = 1;
%         params.DualReductions = 0;
%         params.NumericFocus = 3;
% res = gurobi(model,params);

%u = res.x;
%o.Jpsi = o.res.objval;
for i = 1:p.gn.n
    ui{i} = u(idStart(i):idEnd(i),1);
    
    S = zeros(1,1+p.gn.noN(i));
    S(1) = 1;
    Spsi = [kron(eye(p.h),S) zeros(p.h,p.gn.noN(i))];
    o.psi{i} = Spsi*ui{i};
    
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        
        S = zeros(1,1+p.gn.noN(i));
        S(1+jj) = 1;
        Stau = [kron(eye(p.h),S) zeros(p.h,p.gn.noN(i))];
        o.tau{i,j} = Stau*ui{i};
        
    end
    o.spsi{i} = [zeros(p.gn.noN(i),p.h*(1+p.gn.noN(i))) eye(p.gn.noN(i))]*ui{i};
end


end