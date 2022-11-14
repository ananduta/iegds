function [o,e] = solve_pressure2(o,p)
% Compute pressure decisions with minimum violation of gas-flow PWA model
% Centralized scheme
% W. Ananduta
% 11/10/2021

eps = 1e-7;



E = [];
b = [];
c = 1;
for i=1:p.n
    
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        if j > i
            % define the parameters of the affine function at each region
            nf = @(y) y^2/p.gn.cf(i,j);
            pwaf = pwa_approx_nf(p.r,-p.gn.phi_max(i,j),p.gn.phi_max(i,j),nf);

            b_ij = zeros(p.h,1);
            for m = 1:p.gn.r
                b_ij = b_ij - o.delta_s{i,j}{m}.*(pwaf.a(m)*o.phi{i,j} + pwaf.b(m)*ones(p.h,1));
            end

            A_ij = diag(2*o.delta_psi_s{i,j}-ones(p.h,1));
            E_ij = zeros(1,p.n);
            E_ij(i) = A_ij;
            E_ij(j) = -A_ij;
            
            E = [E;E_ij];
            
            b = [b;b_ij];
        end
            
    end
    

    
end
Et = kron(E,eye(p.h));
e.Et = Et;


e.b = b;
o.psi_c = -pinv(e.Et)*e.b;
% solve with Yalmip
psi = sdpvar(p.n*p.h,1);

J = 0.5*(Et*psi + b)'*(Et*psi + b);

cons = [];
for i=1:p.n
    cons = [cons, psi(p.h*(i-1)+1:p.h*i,1) <=p.gn.psi_max(i)*ones(p.h,1)];
    cons = [cons, psi(p.h*(i-1)+1:p.h*i,1) >= p.gn.psi_min(i)*ones(p.h,1)];
end



% options = sdpsettings('verbose',1,'solver','quadprog');%,'gurobi.OptimalityTol',1e-9,'osqp.eps_abs',1e-10);
options = sdpsettings('verbose',1,'solver','osqp','osqp.eps_abs',1e-9,'osqp.eps_rel',1e-9);
sol = optimize(cons,J,options);

o.psi_solver = value(psi);

for i = 1:p.gn.n
    o.psi{i} = value(psi(p.h*(i-1)+1:p.h*i,1) );
end

% calculate violation
gfv_v=[];
for i = 1:p.gn.n
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        if j > i
            % define the parameters of the affine function at each region
            nf = @(y) y^2/p.gn.cf(i,j);
            pwaf = pwa_approx_nf(p.r,-p.gn.phi_max(i,j),p.gn.phi_max(i,j),nf);

            b_ij = zeros(p.h,1);
            for m = 1:p.gn.r
                b_ij = b_ij - o.delta_s{i,j}{m}.*(pwaf.a(m)*o.phi{i,j} + pwaf.b(m)*ones(p.h,1));
            end

            A_ij = diag(2*o.delta_psi_s{i,j}-ones(p.h,1));


            er = A_ij*(o.psi{i}-o.psi{j}) + b_ij;
            
            er_in = A_ij*(p.m.Spsi{i}*o.u{i} - p.m.Spsi{j}*o.u{j}) + b_ij;
            
            gfv_M(i,j) = sqrt(er'*er);
            gfv_v = [gfv_v;sqrt(er'*er)];
            
            gfv_Min(i,j) = sqrt(er_in'*er_in);
            gfv_vin = [gfv_v;sqrt(er_in'*er_in)];
        end
    end
end
o.gfv_M = gfv_M;
o.gfv_v = gfv_v;
o.gfv_max = max(gfv_v);


o.gfv_Min = gfv_Min;
o.gfv_vin = gfv_vin;
o.gfv_maxin = max(gfv_vin);
end