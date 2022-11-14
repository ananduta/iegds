function u_i = init_prosumer(p,i)
% Initialize prosumers; ED IEGDS
% W. Ananduta
% 11/08/2021


% Local constraints
A = p.tn.A_ineq{i};
b = p.tn.b_ineq{i};


% Coupling constraints
Aco1 = p.tn.Sgu{i};
bco1 = p.tn.eps_g*ones(p.h,1);

%Sphia = [0 ones(1,p.gn.noN(i))];
%Sphiat = kron(eye(p.h),Sphia);
%bco1 = bco1 + Sphiat*s.u_no{i}(:,k) - p.gn.Gdem(i,:)';

Aco2 = -Aco1;
bco2 = p.tn.eps_g*ones(p.h,1);% - Sphiat*s.u_no{i}(:,k)+ p.gn.Gdem(i,:)';


for jj = 1:length(p.tn.N{i})
    j = p.tn.N{i}(jj);
    
    Atr1 = p.tn.Str{i,j};
    btr1 = p.tn.eps_tr*ones(p.h,1);
%    btr1 = btr1 - p.tn.Str{j,i}*s.u{j}(:,k);
    
    Atr2 = -p.tn.Str{i,j};
    btr2 = p.tn.eps_tr*ones(p.h,1);
%    btr2 = btr2 + p.tn.Str{j,i}*s.u{j}(:,k);
    
    Atr{jj} = [Atr1;Atr2];
    btr{jj} = [btr1;btr2];
end
Aco3 = cat(1,Atr{:});
bco3 = cat(1,btr{:});


H = zeros(size(p.tn.A_ineq{i},2));
f = zeros(size(p.tn.A_ineq{i},2),1);
A = sparse([A;Aco1;Aco2;Aco3]);
b = sparse([b;bco1;bco2;bco3]);

% Solve with osqp
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
    u_i = res.x;
%    Ji = res.obj_val;
%     %Assigning the decisions of agent i
%     %s.ph_mg(:,t+1) = 0;
%     %dim_prev_ag = 0;
%     %for i = 1:np.n
%     %    dim_i = np.h*(3+sum(np.Adj(i,:)));
%      %   u_i = u_all(dim_prev_ag+1:dim_prev_ag+dim_i);
%         s.u{i}(:,k+1) = u_i;
%         s.p_di{i}(:,k+1) = p.tn.Sdi{i}*u_i;
%         s.p_st{i}(:,k+1) = p.tn.Sst{i}*u_i;
%         s.p_mg{i}(:,k+1) = p.tn.Smg{i}*u_i;
%         s.d_gu{i}(:,k+1) = p.tn.Sgu{i}*u_i;
%         for jj = 1:length(p.tn.N{i})
%             j = p.tn.N{i}(jj);
%             s.p_tr{i,j}(:,k+1) = p.tn.Str{i,j}*u_i;
%         end
%     %    s.ph_mg(:,t+1) = s.ph_mg(:,t+1) + s.pmg{i}(:,t+1);
%     %    dim_prev_ag = dim_prev_ag+dim_i;
%     %end
end


end