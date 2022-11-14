function s = mi_pwa(p)

p = pwa_gasflow(p);


yalmip('clear')
for i = 1:p.gn.n
    u{i} = sdpvar(p.h*(2+p.gn.noN(i)),1);
    z{i} = sdpvar(p.h*(1+p.gn.r)*p.gn.noN(i),1);
    d{i} = binvar(p.h*(1+3*p.gn.r)*p.gn.noN(i),1);
    %d{i} = sdpvar((1+3*p.gn.r)*p.gn.noN(i),1);
    
    ad{i} = sdpvar(p.h,1);
end
J = 0;
cons = [];

for i = 1:p.gn.n
    
    %% Gas network
    % local cost
    J = J + p.gn.po.c(i)*kron(ones(1,p.h),[1 0 zeros(1,p.gn.noN(i))])*u{i};
    
    
    
    % linear constraints (3)-(5)
    cons = [cons,p.gn.Aineq{i}*u{i} <= p.gn.bineq{i} ];
    
    % local constraints (18)-(26), 
    cons = [cons, p.gn.Gu{i}*u{i} + p.gn.Gz{i}*z{i} + p.gn.Gd{i}*d{i} <= p.gn.gt{i}];
    
    % gas balance (2)
    % local
    cons = [cons,p.gn.Aeq{i}*u{i} - ad{i} == p.gn.beq{i} ];
      
    % active demand
    if isempty(find(p.gn.Ngu==i))
        cons = [cons, ad{i} == 0];
    end
    
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
    
    
    
    
    
    
end


%% Electrical network
for i = 1:p.en.n
    v{i} = sdpvar(p.h*2,1);
end

for i = 1:p.en.n
    % local cost
    if isempty(find(p.en.Ngu == i))
        for k = 1:p.h
            J = J + p.en.q1(i)*v{i}(2*(k-1)+1)^2 + p.en.q2(i)*v{i}(2*(k-1)+1)+ p.en.q3(i);
        end
    end
    
    % local constraints
    cons = [cons,  v{i} <= kron(ones(p.h,1),[p.en.pmax(i);p.en.theta_max])];
    cons = [cons, -v{i} <= kron(ones(p.h,1),[-p.en.pmin(i);-p.en.theta_min])];
    
    
    % DC power-flow coupling constraints
    Cc = kron(eye(p.h),[1 -sum(p.en.B(i,:))])*v{i};
    for jj = 1:p.en.noN(i)
        j = p.en.N{i}(jj);
        Cc = Cc + kron(eye(p.h),[0 p.en.B(i,j)])*v{j};
    end
    
    cons = [cons, Cc == p.en.d(:,i)];
    
    % coupling between gas and electrical network
    if ~isempty(find(p.en.Ngu == i))
        Gge = kron(eye(p.h),[1 0]);
        id_i = find(p.en.Ngu==i);
        i_g = p.gn.Ngu(id_i);
%       cons = [cons, ad{i_g}- p.en.q1(i)*(Gge*v{i})^2 - p.en.q2(i)*Gge*v{i} - p.en.q3(i) >= 0];
         cons = [cons, ad{i_g} - p.en.q2(i)*Gge*v{i} - p.en.q3(i) == 0];
    end
    
end
ops = sdpsettings('verbose',1,'solver','gurobi','gurobi.OptimalityTol',1e-9,'osqp.eps_abs',1e-10);
ops = sdpsettings(ops,'solver','gurobi');
ops = sdpsettings(ops,'gurobi.OptimalityTol',1e-9);
ops = sdpsettings(ops,'gurobi.PerturbValue',1e-6);
%ops = sdpsettings(ops,'gurobi.MIPGap',1e-5);
sol = optimize(cons,J,ops);

for i = 1:p.gn.n
    us{i} = value(u{i});
    zs{i} = value(z{i});
    ds{i} = value(d{i});
    ads{i} = value(ad{i});
end

for i = 1:p.en.n
    vs{i} = value(v{i});
end

s.u = us;
s.z = zs;
s.d = ds;
s.ad = ads;
s.v = vs;
s.J = value(J);

% Evaluate satisfaction of gas flow equation
c=1;
for i = 1:p.gn.n
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
               
        % Gas-flow (32)
        gfu = zeros(1,2+p.gn.noN(i));
        gfu(2+jj) = 1;
        gfz = zeros(1,p.gn.noN(i));
        gfz(jj) = 1;
        if us{i}(2) >= us{j}(2)
            sig = -1;
        else
            sig = 1;
        end
        gflow(c) = 1/sqrt(p.gn.pc.cf(i,j))*(us{i}(2+jj)) - sig*sqrt(abs(us{i}(2)-us{j}(2)));
        g_fl(i,j) = gflow(c);
        c = c+1;
    end
    
end
s.error_gf = mean(mean(abs(gflow)))

end