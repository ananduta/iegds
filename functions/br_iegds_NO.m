function [u_no,Jval] = br_iegds_NO(s,p,k)
% Best response for network operator; ED IEGDS
% W. Ananduta
% 11/08/2021


yalmip('clear')
for i = 1:p.gn.n
    u{i} = sdpvar(p.h*(1+p.gn.noN(i)),1);
    z{i} = sdpvar(p.h*(1+p.gn.r)*p.gn.noN(i),1);
    d{i} = binvar(p.h*(1+3*p.gn.r)*p.gn.noN(i),1);

end
J = 0;
cons = [];

for i = 1:p.gn.n
    
    %% Gas network
    % local cost
    q = kron(ones(1,p.h),[0;ones(p.gn.noN(i),1)]);
    Q = diag(q);
    J = J + u{i}'*Q*u{i};
     
    
    % linear constraints (3)-(5)
    cons = [cons,p.gn.Aineq{i}*u{i} <= p.gn.bineq{i} ];
    
    % local constraints (18)-(26), 
    cons = [cons, p.gn.Gu{i}*u{i} + p.gn.Gz{i}*z{i} + p.gn.Gd{i}*d{i} <= p.gn.gt{i}];
    
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        
        % simplex (15)
        cons = [cons, p.gn.Gsdi{i,j}*d{i} == p.gn.gsdi{i,j}];
        
        % coupling constraints
        % Gas-flow (14)
        cons = [cons, p.gn.GFui_eq{i,j}*u{i} + p.gn.GFzi_eq{i,j}*z{i} + p.gn.GFdi_eq{i,j}*d{i}+ p.gn.GFuj_eq{i,j}*u{j} + p.gn.GFzj_eq{i,j}*z{j}== zeros(size(p.gn.GFui_eq{i,j},1),1)];
        
        % (16)-(17)
        cons = [cons, p.gn.Hui{i,j}*u{i} + p.gn.Hzi{i,j}*z{i} + p.gn.Hdi{i,j}*d{i} + p.gn.Huj{i,j}*u{j} <= p.gn.hij{i,j}];
        
        % between neighboring nodes (reciprocity constraints)
        cons = [cons,p.gn.GRcu{i,j}*u{i} + p.gn.GRcuj{i,j}*u{j} == 0];
        
        % replacement of reciprocity const
        cons = [cons,p.gn.GRu{i,j}*u{i} + p.gn.GRd{i,j}*d{i} <= p.gn.gR{i,j}];
    end
    
    % Coupling constraints

    Aco1 = kron(eye(p.h),[0, -ones(1,p.gn.noN(i))]);
    bco1 = p.tn.eps_g*ones(p.h,1);
    bco1 = bco1 - p.gn.Gdem(i,:)' - p.tn.Sgu{i}*s.u{i}(:,k);
    
    Aco2 = kron(eye(p.h),[0, ones(1,p.gn.noN(i))]);
    bco2 = p.tn.eps_g*ones(p.h,1);
    bco2 = bco2 + p.gn.Gdem(i,:)' + p.tn.Sgu{i}*s.u{i}(:,k);
    
    cons = [cons, [Aco1;Aco2]*u{i} <= [bco1;bco2]]; 
    
    
end


% Solve with yalmip
ops = sdpsettings('verbose',1,'solver','gurobi','gurobi.OptimalityTol',1e-9,'osqp.eps_abs',1e-10);
ops = sdpsettings(ops,'solver','gurobi');
ops = sdpsettings(ops,'gurobi.OptimalityTol',1e-9);
ops = sdpsettings(ops,'gurobi.PerturbValue',1e-6);
%ops = sdpsettings(ops,'gurobi.MIPGap',1e-5);
sol = optimize(cons,J,ops);

for i = 1:p.gn.n
    u_no = value(u{i});
    
end
 
Jval = value(J);
end