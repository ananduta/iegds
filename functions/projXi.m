function u_out = projXi(y,np,i)
% Projection of y onto the local set
% W. Ananduta
% 15/07/2021



A = np.A_ineq{i};
b = np.b_ineq{i};
Aeq = np.Aeq{i};
beq = np.beq{i};

%quadprog
%H = np.H{i};
%f = np.c{i};
%options = optimset('Display','off');
%options = optimset('Display','on');
%[u_i,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq,[],[],[],options);

%lsqlin
%C = np.H_half{i};
%d = -np.H_half_inv{i}*np.c{i};
%options = optimset('Display','off');
%[u_i,~,~,exitflag] = lsqlin(C,d,A,b,Aeq,beq,[],[],[],options);

%OSQP
H = sparse(eye(length(y)));
f = -y;
A = sparse([A;Aeq;-Aeq]);
b = sparse([b;beq;-beq]);

% Create an OSQP object
prob = osqp;

settings = prob.default_settings();
settings.eps_abs= 1e-10;
settings.eps_rel= 1e-10;
%settings.max_iter = 1e5;
settings.verbose = 0;
% Setup workspace and change alpha parameter
prob.setup(H, f, A, [], b, settings);
res = prob.solve();
u_i = res.x;

%if exitflag ~= 1 && exitflag ~= 0 %quadprog
%     disp(st.problem): in case solver cannot find the solution
%    exitflag
%    disp('not solved properly')
%    st.problem

if res.info.status_val ~= 1 && res.info.status_val ~= -2    
    res.info.status_val
    disp('not solved properly')
    pause
%     s.u{i}(:,t+1) = s.u{i}(:,t);
%     s.pmg{i}(:,t+1) = s.pmg{i}(:,t);
else
    %Assigning the decisions of agent i
    %s.ph_mg(:,t+1) = 0;
    %dim_prev_ag = 0;
    %for i = 1:np.n
    %    dim_i = np.h*(3+sum(np.Adj(i,:)));
     %   u_i = u_all(dim_prev_ag+1:dim_prev_ag+dim_i);
        u_out = u_i;
        
    %    s.ph_mg(:,t+1) = s.ph_mg(:,t+1) + s.pmg{i}(:,t+1);
    %    dim_prev_ag = dim_prev_ag+dim_i;
    %end
end


end