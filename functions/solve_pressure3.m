function [o,e] = solve_pressure3(o,p)
% Compute pressure decisions with minimum violation of gas-flow PWA model
% Centralized scheme
% W. Ananduta
% 02/02/2022

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
            E_ij = zeros(p.h,p.h*p.n);
            E_ij(:,p.h*(i-1)+1:p.h*i) = A_ij;
            E_ij(:,p.h*(j-1)+1:p.h*j) = -A_ij;
            
            E = [E;E_ij];
            
            b = [b;b_ij];
        end
            
    end
    

    
end
Et = E;
%Et = kron(E,eye(p.h));
e.Et = Et;


e.b = b;
o.psi_c = -pinv(e.Et)*e.b;


% Solve with OSQP

% create matrices for constraints
A1 = [Et -ones(size(Et,1),1)];
b1 = -b;
A2 = [-Et -ones(size(Et,1),1)];
b2 = b;
A3 = [eye(size(Et,2)) zeros(size(Et,2),1)];
b3 = kron(p.gn.psi_max,ones(p.h,1));
A4 = [-eye(size(Et,2)) zeros(size(Et,2),1)];
b4 = kron(-p.gn.psi_min,ones(p.h,1));

Aineq = [A1;A2;A3;A4];
bineq = [b1;b2;b3;b4];

% create matrix for cost
f = [zeros(size(Et,2),1); 1];
H = zeros(size(Et,2)+1);

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

o.psi_solver = res.x;
o.Jpsi_solver = res.info.obj_val;
% [x,fval,exitflag,output] = quadprog(H, f, A, b);
% o.psi_solver = x;
for i = 1:p.gn.n
    o.psi{i} = o.psi_solver(p.h*(i-1)+1:p.h*i,1);
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
            
            gfv_M(i,j) = max(abs(er));
            gfv_v = [gfv_v;max(abs(er))];
            
            gfv_Min(i,j) = max(abs(er));
            gfv_vin = [gfv_v;max(abs(er))];
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