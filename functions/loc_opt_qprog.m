function s = loc_opt_qprog(p,s,k,i)
% Local optimization (using quadprog)
% W. Ananduta
% 02/12/2019

%for i = 1:p.n
    
    % Construct coefficient of the linear cost
    
    % linear term associated with primal cost function
    %cc = 1;
    c = zeros(p.nu(i),1);
    c(1) = p.en.c_dg(i);
    c(2) = p.en.c_st(i);
    

    for jj = 1:length(p.en.N{i})
        j = p.en.N{i}(jj);
        c(4+jj,1) = p.en.c_tr(i,j);
        
    end
    
    c1 = [];
    
    % linear term associated with pmg and dgu
    pmg_i_t = p.m.Smg{i}*s.u{i}(:,k);
    dgu_i_t = p.m.Sgu{i}*s.u{i}(:,k);
    for h = 1:p.h
        c_1 = c;
        c_1(3) = p.en.d_mg*(s.sigma_mg(h,k)-pmg_i_t(h));
        c_1(4) = p.gn.d_gu*(s.sigma_gu(h,k)-dgu_i_t(h));
        c1 = [c1;c_1];
    end 
    
    
    % linear term associated with grid constraints
    c2 = [p.m.Smg{i};-p.m.Smg{i}]'*s.lambda_mg(:,k);
    
    c3 = zeros(p.h*p.nu(i),1);
    % linear term associated with gas coupling constraint
    c3 = [p.m.Sgu{i};-p.m.Sgu{i}]'*s.lambda_gu(:,k);
    
    % linear term associated with reciprocity constraints
    c4 = zeros(p.h*p.nu(i),1);
    for jj=1:length(p.en.N{i})
        j = p.en.N{i}(jj);
        c4 = c4 + p.m.Str{i,j}'*s.mu_tr{i,j}(:,k);           
    end
    
    
    
    c5 = zeros(p.h*p.nu(i),1);
    c6 = zeros(p.h*p.nu(i),1);
    c7 = zeros(p.h*p.nu(i),1);
    for jj=1:length(p.gn.N{i})
        j = p.gn.N{i}(jj);
        
        % linear term associated with reciprocity of \phi
        c5 = c5 + p.m.Sphi{i,j}'*s.mu_phi{i,j}(:,k);     
        
        % linear term associated with gas-flow equation
        c6 = c6 + p.m.Hc{i,j}{1}'*s.mu_gf{i,j}(:,k) + p.m.Hc{j,i}{2}'*s.mu_gf{j,i}(:,k);
        
        % linear term associated with psi coupling const
        c7 = c7 + p.m.G{i,j}{1}'*s.lambda_psi{i,j}(:,k) + p.m.G{j,i}{2}'*s.lambda_psi{j,i}(:,k);
    end
    
    
    % linear term associated with proximal term
    c8 = -p.Alpha{i}*s.u{i}(:,k);
    
    p.m.c{i} = c1+  c2+  c3+  c4 + c5 + c6 + c7+  c8;
    
    
    
%end

A = p.m.Aineq{i};
b = p.m.bineq{i};
Aeq = p.m.Aeq{i};
beq = p.m.beq{i};

%quadprog
H = p.m.H{i};
f = p.m.c{i};
options = optimset('Display','off');
%options = optimset('Display','on');
[u_i,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq,[],[],[],options);

%lsqlin
%C = p.H_half{i};
%d = -p.H_half_inv{i}*p.c{i};
%options = optimset('Display','off');
%[u_i,~,~,exitflag] = lsqlin(C,d,A,b,Aeq,beq,[],[],[],options);

% %OSQP
% H = p.m.H{i};
% f = p.m.c{i};
% A = sparse([A;Aeq;-Aeq]);
% b = sparse([b;beq;-beq]);

% % Create an OSQP object
% prob = osqp;
% 
% settings = prob.default_settings();
% settings.eps_abs= 1e-7;
% settings.eps_rel= 1e-7;
% settings.max_iter = 1e5;
% settings.verbose = 0;
% % Setup workspace and change alpha parameter
% prob.setup(H, f, A, [], b, settings);
% res = prob.solve();
% u_i = res.x;

if exitflag ~= 1 && exitflag ~= 0 %quadprog
%    disp(st.problem) %in case solver cannot find the solution
   exitflag
   disp('not solved properly')
%   st.problem

% if res.info.status_val ~= 1 && res.info.status_val ~= -2 % && res.info.status_val ~= 2  
%     res.info.status_val
%     disp('not solved properly')
%     pause
%     s.u{i}(:,t+1) = s.u{i}(:,t);
%     s.pmg{i}(:,t+1) = s.pmg{i}(:,t);
%     u_i = s.u{i}(:,k);
%     s.u{i}(:,k+1) = u_i;
%         s.p_di{i}(:,k+1) = p.m.Sdg{i}*u_i;
%         s.p_st{i}(:,k+1) = p.m.Sst{i}*u_i;
%         s.p_mg{i}(:,k+1) = p.m.Smg{i}*u_i;
%         s.d_gu{i}(:,k+1) = p.m.Sgu{i}*u_i;
%         for jj = 1:length(p.en.N{i})
%             j = p.en.N{i}(jj);
%             s.p_tr{i,j}(:,k+1) = p.m.Str{i,j}*u_i;
%         end
%         
%         % y
%         s.psi{i}(:,k+1) = p.m.Spsi{i}*u_i;
%         for jj = 1:length(p.gn.N{i})
%             j = p.gn.N{i}(jj);
%             s.phi{i,j}(:,k+1) = p.m.Sphi{i,j}*u_i;
%         end
else
    %Assigning the decisions of agent i
    %s.ph_mg(:,t+1) = 0;
    %dim_prev_ag = 0;
    %for i = 1:p.n
    %    dim_i = p.h*(3+sum(p.Adj(i,:)));
     %   u_i = u_all(dim_prev_ag+1:dim_prev_ag+dim_i);
        s.u{i}(:,k+1) = u_i;
        s.p_di{i}(:,k+1) = p.m.Sdg{i}*u_i;
        s.p_st{i}(:,k+1) = p.m.Sst{i}*u_i;
        s.p_mg{i}(:,k+1) = p.m.Smg{i}*u_i;
        s.d_gu{i}(:,k+1) = p.m.Sgu{i}*u_i;
        for jj = 1:length(p.en.N{i})
            j = p.en.N{i}(jj);
            s.p_tr{i,j}(:,k+1) = p.m.Str{i,j}*u_i;
        end
        
        % y
        s.psi{i}(:,k+1) = p.m.Spsi{i}*u_i;
        for jj = 1:length(p.gn.N{i})
            j = p.gn.N{i}(jj);
            s.phi{i,j}(:,k+1) = p.m.Sphi{i,j}*u_i;
        end
        
    %    s.ph_mg(:,t+1) = s.ph_mg(:,t+1) + s.pmg{i}(:,t+1);
    %    dim_prev_ag = dim_prev_ag+dim_i;
    %end
end


end