function  s = init_u(np,s,i)
% Initialization of u
% W. Ananduta
% 25/07/2019

% Inputs:
% np: Network parameters
% s: simulation results
% i: agent

yalmip('clear');
cons = [];
Jc = 0;

Ni = sum(np.Adj(i,:));

z = sdpvar(3+Ni,np.h);

z_c = z(:);

x = sdpvar(np.h+1);
x(1) = np.x0(i);

% constraints
for k=1:np.h
    
    % Load power balance
    cons = [cons, ones(1,3+Ni)*z(:,k) == np.Pd(i,k)]; 
    
    % DG unit
    cons = [cons, np.pdg_min(i) <= z(1,k) <= np.pdg_max(i)]; %(
   
    % Storage unit
    cons = [cons, -np.p_ch(i) <= z(2,k) <= np.p_dh(i)]; 
    cons = [cons, x(k+1) == np.a(i)*x(k) + np.b(i)*z(2,k)];
    cons = [cons, np.x_min(i) <= x(k+1) <= np.x_max(i)];
    
    % Power traded with neighbors
    cc=1;
    for j=1:np.n
        if np.Adj(i,j) == 1
            cons = [cons, -np.pt_max(i,j) <= z(3+cc,k) <= np.pt_max(i,j)];
            cc = cc+1;
        end
    end   
end

% solve the optimization problem
ops = sdpsettings('verbose',0,'solver','quadprog');
optimize(cons,Jc,ops);

u_i = value(z);
s.u{i}(:,1) = u_i(:);
s.pmg{i}(:,1) = u_i(3,:)';
end

