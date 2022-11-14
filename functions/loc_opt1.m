function  s = loc_opt1(np,s,i,t)
% Local optimization
% W. Ananduta
% 24/07/2019

% Inputs:
% np: Network parameters
% s: simulation results
% i: agent

yalmip('clear');
cons = [];
Jc = 0;

Ni = sum(np.Adj(i,:));

z = sdpvar(np.h*(3+Ni),1);

x0 = np.x0(i);


% Cost function

p_mg = np.Smg{i}*z;

f_mg = np.q_mg*s.ph_mg(:,t)'*p_mg; %cost of trading with main grid

% cost associated to reciprocity consts.
f_rep = 0;
for j=1:np.n
    if np.Adj(i,j) == 1
        f_rep = f_rep+s.mu{i,j}(:,t)'*np.Str{i,j}*z;
    end
end

Ji = z'*np.Qh{i}*z + np.ch{i}'*z + f_mg;

  Jc = Ji + s.lambda{i}(:,t)'*[p_mg;-p_mg]+f_rep+0.5*(z'-s.u{i}(:,t)')*np.a_bar{i}*(z-s.u{i}(:,t));
%  Jc = Ji + s.lambda{i}(:,t)'*[p_mg;-p_mg]+f_rep;
% constraints U_i
% Load power balance
cons = [cons,  np.E1t{i}*z <=  np.Pd(i,1:np.h)'];
cons = [cons, -np.E1t{i}*z <= -np.Pd(i,1:np.h)'];

% DG unit
cons = [cons,  np.E2t{i}*z <=  np.pdg_max(i)*ones(np.h,1)];
cons = [cons, -np.E2t{i}*z <= -np.pdg_min(i)*ones(np.h,1)];

% Storage unit
cons = [cons,  np.E3t{i}*z <= np.p_dh(i)*ones(np.h,1)];
cons = [cons, -np.E3t{i}*z <= np.p_ch(i)*ones(np.h,1)];
if np.st_un(i) == 1
    cons = [cons, -np.At{i}*x0+np.Bt{i}*z <= -np.x_min(i)*ones(np.h,1)]; 
    cons = [cons,  np.At{i}*x0+np.Bt{i}*z <=  np.x_max(i)*ones(np.h,1)];
end


% Power traded
cons = [cons, np.E4t{i}*z <= np.Pt_m{i}];
cons = [cons,-np.E4t{i}*z <= np.Pt_m{i}];

% solve the optimization problem
ops = sdpsettings('verbose',0,'solver','quadprog');
st=optimize(cons,Jc,ops);

if st.problem ~= 0
%     disp(st.problem)
    disp('not solved properly') 
    s.u{i}(:,t+1) = s.u{i}(:,t);
    s.pmg{i}(:,t+1) = s.pmg{i}(:,t);
else
    u_i = value(z);
    s.u{i}(:,t+1) = u_i(:);
    s.pmg{i}(:,t+1) = np.Smg{i}*u_i;
end

end

