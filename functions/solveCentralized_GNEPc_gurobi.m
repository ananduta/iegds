function [p,o]= solveCentralized_GNEPc_gurobi(p)
% Computing v-GNE of convexified game centrally via minimizing potential function
% W. Ananduta
% 07/02/2022
    

    % Identify dimensions of decision variables (per time step, h)
    for i = 1:p.n
    p.nx(i) = 8+ p.en.noN(i);
    p.nt(i) = p.gn.noN(i);
    p.ny(i) = 2+(2+p.gn.r)*p.gn.noN(i)+p.nt(i);
    p.nz(i) = (1+3*p.gn.r)*p.gn.noN(i);
    p.nu(i) = p.nx(i) + p.ny(i) + p.nz(i);
    end
    p.id_dg = 1;
    p.id_ch = 2;
    p.id_dh = 3;
    p.id_mg = 4;
    p.id_gu = 5;
    p.id_th = 6;
    p.id_v = 7;
    p.id_eg = 8;
    
    % Generate agents' indices in the concatenated decision variable u
    p.idAgent = index_decision(p.nu,p.n,p.h);
    idStart = p.idAgent(:,1);
    idEnd = p.idAgent(:,2);
    
    % Generate matrices for  constraints       
    p = build_mat_exP2P_pen(p);  
    
    % local constraints
    A_all = [];
    b_all = [];
    for i = 1:p.n
        A = p.m.Aineq{i};
        b = p.m.bineq{i};
        Aeq = p.m.Aeq{i};
        beq = p.m.beq{i};
        
        At{i} = [A;Aeq;-Aeq];
        bt{i} = [b;beq;-beq];
        
        Ai = sparse(zeros(size(At{i},1),idEnd(p.n)));
        Ai(:,idStart(i):idEnd(i)) = At{i};
        bi = bt{i};
        
        A_all = sparse([A_all; Ai]);
        b_all = [b_all; bi];
    end
    
    % coupling constraints
    for i = 1:p.n
        for jj=1:length(p.en.N{i})
            j = p.en.N{i}(jj);
            
            %power flow constraints
            Ac = sparse(zeros(p.h,idEnd(p.n)));
            Ac(:,idStart(i):idEnd(i)) = p.m.PF{i,j}{1};
            Ac(:,idStart(j):idEnd(j)) = p.m.PF{i,j}{2};
            
            Ac = [Ac;-Ac];
            bc = zeros(2*p.h,1);
            
            A_all = sparse([A_all;Ac]);
            b_all = [b_all;bc];
        end
        
        for jj = 1:length(p.gn.N{i})
            j = p.gn.N{i}(jj);
            
            % reciprocity constraints
            Ac = sparse(zeros(p.h,idEnd(p.n)));
            Ac(:,idStart(i):idEnd(i)) = p.m.Sphi{i,j};
            Ac(:,idStart(j):idEnd(j)) = p.m.Sphi{j,i};
            
            Ac = [Ac;-Ac];
            bc = zeros(2*p.h,1);
            
            A_all = sparse([A_all;Ac]);
            b_all = [b_all;bc];
            
            % gas flow constraints
            Ac = sparse(zeros(p.h,idEnd(p.n)));
            Ac(:,idStart(i):idEnd(i)) = p.m.Hc{i,j}{1};
            Ac(:,idStart(j):idEnd(j)) = p.m.Hc{i,j}{2};
            
            Ac = [Ac;-Ac];
            bc = zeros(2*p.h,1);
            
            A_all = sparse([A_all;Ac]);
            b_all = [b_all;bc];
            
            % psi coupling constraints
            Ac = sparse(zeros(size(p.m.G{i,j}{1},1),idEnd(p.n)));
            Ac(:,idStart(i):idEnd(i)) = p.m.G{i,j}{1};
            Ac(:,idStart(j):idEnd(j)) = p.m.G{i,j}{2};
            
            bc = p.m.g{i,j};
            
            A_all = sparse([A_all;Ac]);
            b_all = [b_all;bc];
            
            
        end
        
                
    end
    % aggregative constraints
    % electrical network p_mg
    Ac = cat(2,p.m.Smg{:});
    Ac = [Ac;-Ac];
    bc = [p.en.pmg_max*ones(p.h,1);-p.en.pmg_min*ones(p.h,1)];
    A_all = sparse([A_all;Ac]);
    b_all = [b_all;bc];
    
    Ac = cat(2,p.m.Sgu{:});
    Ac = [Ac;-Ac];
    bc = [p.gn.dg_max*ones(p.h,1);-p.gn.dg_min*ones(p.h,1)];
    A_all = sparse([A_all;Ac]);
    b_all = [b_all;bc];
    
    %p.m.A_all = A_all;
    %p.m.b_all = b_all;
    %p.m.At = At;
    %p.m.bt = bt;
    
 
    
    
    % Generate matrices for cost function
    
    p = alg_param3(p); % nothing is used from the output of this function, but must be run.
    
    p = build_mat_cost_pen(p);
    
    % linear term
    f = cat(1,p.m.ch{:});
    
    % quadratic term
    D = sparse(zeros(idEnd(p.n)));
    Hloc = sparse(zeros(idEnd(p.n)));
    Qe = eye(p.h)*p.en.d_mg;
    Qg = eye(p.h)*p.gn.d_gu;
    for i = 1:p.n
        for j = 1:p.n
            if i==j
                D(idStart(i):idEnd(i),idStart(j):idEnd(j)) = 2*(p.m.Smg{i}'*Qe*p.m.Smg{j}+p.m.Sgu{i}'*Qg*p.m.Sgu{j});
            else
                D(idStart(i):idEnd(i),idStart(j):idEnd(j)) = (p.m.Smg{i}'*Qe*p.m.Smg{j}+p.m.Sgu{i}'*Qg*p.m.Sgu{j});
            end
        end
        Hloc(idStart(i):idEnd(i),idStart(i):idEnd(i)) = 2*p.m.Qph{i};
    end
    %p.m.Qe = Qe;
    %p.m.Qg = Qg;
    %p.m.D = D;
    %p.m.Hloc = Hloc;
    
    %p.m.Hall = sparse(p.m.Hloc + p.m.D); 
    
    %J = x'*Hloc*x + f'*x;
    
    %% solve with Gurobi
    model.Q = Hloc+D;
    model.obj = f;
    
    
    % lower bound
    model.lb = -inf*ones(idEnd(p.n),1);
    
    % constraints
    model.A = A_all;
    model.rhs = b_all;
    
    

    res = gurobi(model);
    o.res = res;
    u = res.x;
    o.u_all= u;
    o.res = res;
    % Assign solution
    for i = 1:p.n
        o.u{i} = u(idStart(i):idEnd(i),1);

        o.p_di{i} = p.m.Sdg{i}*o.u{i};
        o.p_ch{i} = p.m.Sch{i}*o.u{i};
        o.p_dh{i} = p.m.Sdh{i}*o.u{i};
        o.p_mg{i} = p.m.Smg{i}*o.u{i};
        o.d_gu{i} = p.m.Sgu{i}*o.u{i};
        o.th{i} = p.m.Sth{i}*o.u{i}/p.scaling;
        o.v{i} = p.m.Sv{i}*o.u{i}/p.scaling;
        for jj = 1:length(p.en.N{i})
            j = p.en.N{i}(jj);
            o.p_l{i,j} = p.m.Spl{i,j}*o.u{i};
        end

    % y
        o.psi{i} = p.m.Spsi{i}*o.u{i};

        for jj = 1:length(p.gn.N{i})
            j = p.gn.N{i}(jj);
            o.phi{i,j} = p.m.Sphi{i,j}*o.u{i};
            o.pen{i,j} = p.m.Spen{i,j}*o.u{i};
        end
    end
    
end