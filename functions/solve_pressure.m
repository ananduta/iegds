function o = solve_pressure(o,p)
% Compute pressure decisions with minimum violation of gas-flow PWA model
% Centralized scheme
% W. Ananduta
% 11/10/2021

eps = 1e-7;
yalmip('clear')
for i=1:p.n
    psi{i} = sdpvar(p.h,1);
end

J = 0;
cons = [];

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

           J = J + 0.5*psi{i}'*(A_ij'*A_ij)*psi{i} - psi{j}'*(A_ij'*A_ij)*psi{i} ...
               + 0.5*psi{j}'*(A_ij'*A_ij)*psi{j} + b_ij'*A_ij*psi{i} - b_ij'*A_ij*psi{j};
            
%             J = J + 0.5*(A_ij*psi{i} - A_ij*psi{j} + b_ij)'*(A_ij*psi{i} - A_ij*psi{j} + b_ij);
            %cons = [cons, -psi{i} + psi{j} <= - (p.gn.psi_min(i) - p.gn.psi_max(j))*(ones(p.h,1) - o.delta_psi_s{i,j})];
            %cons = [cons, -psi{i} + psi{j} >= eps*ones(p.h,1) + (-p.gn.psi_max(i) + p.gn.psi_min(j)-eps)*o.delta_psi_s{i,j}];
        end
            
    end
    
%     cons = [cons, psi{i} <= p.gn.psi_max(i)*ones(p.h,1)];
%     cons = [cons, psi{i} >= p.gn.psi_min(i)*ones(p.h,1)];
    
end

%options = sdpsettings('verbose',1,'solver','quadprog');%,'gurobi.OptimalityTol',1e-9,'osqp.eps_abs',1e-10);
options = sdpsettings('verbose',1,'solver','osqp','osqp.eps_abs',1e-6,'osqp.eps_rel',1e-6);
sol = optimize(cons,J,options);

for i = 1:p.gn.n
    o.psi{i} = value(psi{i});
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
            gfv_M(i,j) = sqrt(er'*er);
            gfv_v = [gfv_v;sqrt(er'*er)];
        end
    end
end
o.gfv_M = gfv_M;
o.gfv_v = gfv_v;
o.gfv_max = max(gfv_v);

end