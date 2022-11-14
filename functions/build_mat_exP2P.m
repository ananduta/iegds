function p = build_mat_exP2P(p)
% Building matrices for local optimization 
% for prosumers of extended P2P (in the cost and constraints)
% W. Ananduta
% 05/08/2021

%np = p.en;

m = gen_Smat(p);
p.m = m;

for i=1:p.n
    Ni = p.en.noN(i);

    %% Cost function 
    Q = diag([p.en.q_dg(i) p.en.q_st(i) p.en.d_mg p.gn.d_gu zeros(1,p.en.noN(i))]);

    
    Qh{i} = p.m.Sx{i}'*kron(eye(p.h),Q)*p.m.Sx{i};
    p.m.Qh{i} = sparse(Qh{i});
%    p.m.H{i} = 1/p.alpha(i)*eye(size(Qh{i},1)) + 2*Qh{i};
    p.m.H{i} = p.Alpha{i} + 2*Qh{i};
    
%    p.H{i} = p.A{i} + 2*Qh{i}; % Coefficient of the quadratic term
%    p.H_half{i} = sqrt(p.H{i});
%    p.H_half_inv{i} = inv(p.H_half{i});
    
    cc = 1;
    c = [p.en.c_dg(i) p.en.c_st(i) 0 0 zeros(1,p.en.noN(i))]';
    for jj=1:p.en.noN(i)
        j = p.en.N{i}(jj);
        if p.en.Adj(i,j) == 1
            c(4+cc,1) = p.en.c_tr(i,j);
            cc = cc+1;
        end
    end
    p.m.ch{i} = p.m.Sx{i}'*kron(ones(p.h,1),c);
    
    %% Constraints U_i
    
    % Equality constraints
    
    % Electrical power balance (equality constraint) eq. (10)
    E1 = [ones(1,p.nx(i)) zeros(1,p.ny(i)+p.nz(i))];
    E1(1,4) = 0;
    Aeq1 = sparse(kron(eye(p.h),E1));
    beq1 = p.en.Pd(i,1:p.h)';
    
    
    % Gas balance eq. (11)
    E3 = zeros(1,p.nu(i));
    E3(1,4) = -1;
    E3(1,p.nx(i)+2) = 1;
    for j = 1:p.gn.noN(i)
        E3(1,p.nx(i)+2+(2+p.gn.r)*(j-1)+1) = 1;
    end
    Aeq3 = sparse(kron(eye(p.h),E3));
    beq3 = p.gn.Gdem(i,1:p.h)';
    
    
    % Coupling between gas load and power dispatchable unit (3)
    if p.en.dgu_un(i) == 1
        E2 = zeros(1,p.nu(i));
        E2(1,1) = -p.en.q2(i);
        E2(1,4) = 1;
        Aeq2 = sparse(kron(eye(p.h),E2));
        beq2 = p.en.q3(i)*ones(p.h,1);
    else
        % otherwise, set decision of d_gu = 0
        E2 = zeros(1,p.nu(i));
        E2(1,4) = 1;
        Aeq2 = sparse(kron(eye(p.h),E2));
        beq2 = zeros(p.h,1);
    end
    
    
    % g_t for not-connected nodes (17)
    if p.gn.gt(i) == 0
        E4 = zeros(1,p.nu(i));
        E4(1,p.nx(i)+2) = 1;
        Aeq4 = sparse(kron(eye(p.h),E4));
        beq4 = zeros(p.h,1);
        
        
        
    else
        Aeq4 = [];
        beq4 = [];
    end
    

    p.m.Aeq{i} = sparse([Aeq1;Aeq2;Aeq3;Aeq4]);
    p.m.beq{i} = sparse([beq1;beq2;beq3;beq4]);
    
    % Inequality constraints
    
    % DG unit eq(1)
    A1 = [ p.m.Sdg{i};
          -p.m.Sdg{i}];
    b1 = [ p.en.pdg_max(i)*ones(p.h,1);
          -p.en.pdg_min(i)*ones(p.h,1)];
    

    % Storage unit eq(??)
    A2 = [ p.m.Sst{i};
          -p.m.Sst{i}];
    b2 = [p.en.p_dh(i)*ones(p.h,1);
          p.en.p_ch(i)*ones(p.h,1)];

    if p.en.st_un(i) == 1
        % constructing the matrices for the dynamic equation (At and Bt)
        Am = p.en.a_st(i);
        Bm = [0 -1 0 0 zeros(1,p.en.noN(i)) zeros(1,p.ny(i)+p.nz(i))]; %be careful when sampling time is not 1 hour

        Btcol = zeros(p.h,size(Bm,2));
        Bt = zeros(p.h,p.h*size(Bm,2));
        for l=1:p.h
            At(l,:) =Am^l;
            Btcol(l,:) = Am^(l-1)*Bm;
        end

        for l=1:p.h
            Bt(:,(l-1)*(length(Bm))+1:l*length(Bm)) = [zeros((l-1),length(Bm)); Btcol(1:size(Btcol,1)-(l-1),:)];
        end
        A3 = [ Bt; 
              -Bt];
        b3 = [p.en.x_max(i)*ones(p.h,1)-At*p.en.x0(i);
              -p.en.x_min(i)*ones(p.h,1)+At*p.en.x0(i)];
    else
        A3 = [];
        b3 = [];
    end
    
    % Import from main grid 
    A5 = [-p.m.Smg{i}];
    %b5 = [-p.pmg_min_a*ones(p.h,1)];
    b5 = [-zeros(p.h,1)];
    
    % gas demand 
    A6 = [-p.m.Sgu{i}];
    %b5 = [-p.pmg_min_a*ones(p.h,1)];
    b6 = [-zeros(p.h,1)];
    
    % Power traded eq. (8)
    A4 = [zeros(Ni,4) eye(Ni) zeros(Ni,p.ny(i)+p.nz(i))];
    A4 = kron(eye(p.h),A4);
    
    pt_m = zeros(Ni,1);
    for jj=1:p.en.noN(i)
        j = p.en.N{i}(jj);
        pt_m(jj,1) = p.en.pt_max(i,j);
    end
    b4 = kron(ones(p.h,1),pt_m);
    
    
    
    % flow limits (15)
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        
        A7j{jj} = [ p.m.Sphi{i,j};
                   -p.m.Sphi{i,j}];
        b7j{jj} = [ p.gn.phi_max(i,j)*ones(p.h,1);
                    p.gn.phi_max(i,j)*ones(p.h,1)];
    end
    A7 = cat(1,A7j{:});
    b7 = cat(1,b7j{:});
    clearvars('A7j','b7j');
    
    % pressure limits (16)
    A8 = [ p.m.Spsi{i};
          -p.m.Spsi{i}];
    b8 = [ p.gn.psi_max(i)*ones(p.h,1);
          -p.gn.psi_min(i)*ones(p.h,1)];  
      
    p.m.Aineq{i} = sparse([A1;A2;A3;A4;A5;A6;A7;A8]);
    p.m.bineq{i} = sparse([b1;b2;b3;b4;b5;b6;b7;b8]);
    
    
    % PWA gas flow model
    pwa = pwa_mod(p,i);
    
    
    
    p.m.Aeq{i} = [p.m.Aeq{i}; pwa.Aeq];
    p.m.beq{i} = [p.m.beq{i}; pwa.beq];
    p.m.Aineq{i} = [p.m.Aineq{i}; pwa.Aineq];
    p.m.bineq{i} = [p.m.bineq{i}; pwa.bineq];
    
    
    % Coupling constraints of gas flow model
    eps = 1e-7;
    
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        
        % gas-flow equations (35)
        % define the parameters of the affine function at each region
        nf = @(y) y^2/p.gn.cf(i,j);
        pwaf = pwa_approx_nf(p.r,-p.gn.phi_max(i,j),p.gn.phi_max(i,j),nf);
        
        
        Hih = zeros(1,p.nu(i));
        
        for m =1:p.gn.r
            id_delta_m = p.nx(i)+p.ny(i)+(1+3*p.gn.r)*(jj-1) + 1+3*(m-1) + 3;
            
            id_z_m = p.nx(i) + 2 + (2+p.gn.r)*(jj-1) + 2 + m;
            
            Hih(1,id_z_m) = pwaf.a(m);
            Hih(1,id_delta_m) = pwaf.b(m);
            
            
        end
        
        id_z_psi = p.nx(i) + 2 + (2+p.gn.r)*(jj-1) + 2;
        Hih(1,id_z_psi) = -2;
        
        p.m.Hc{i,j}{1} = kron(eye(p.h),Hih) + p.m.Spsi{i};
        
        
        
        ii = find(p.gn.N{j}==i);
        Hjh = zeros(1,p.nu(j));
        
        id_z_psi_j = p.nx(j) + 2 + (2+p.gn.r)*(ii-1) + 2;
        Hjh(1,id_z_psi_j) = -2;
        
        p.m.Hc{i,j}{2} = kron(eye(p.h),Hjh) + p.m.Spsi{j};
        p.m.hc{i,j} = 0*ones(p.h,1);
        
        % (38)-(39)
        id_delta_psi = p.nx(i)+p.ny(i)+(1+3*p.gn.r)*(jj-1) + 1;
        
        Gih38 = zeros(1,p.nu(i));        
        Gih38(1,id_delta_psi)= -(p.gn.psi_min(i) - p.gn.psi_max(j));
        gih38 = -(p.gn.psi_min(i) - p.gn.psi_max(j));
        
        Gih39 = zeros(1,p.nu(i));        
        Gih39(1,id_delta_psi)= -(p.gn.psi_max(i) - p.gn.psi_min(j)) - eps;
        gih39 = -eps;
        
        p.m.G{i,j}{1} = [kron(eye(p.h),Gih38) - p.m.Spsi{i};
                         kron(eye(p.h),Gih39) + p.m.Spsi{i}];
        p.m.G{i,j}{2} = [ p.m.Spsi{j};
                         -p.m.Spsi{j}];
        p.m.g{i,j} = [gih38*ones(p.h,1);
                      gih39*ones(p.h,1)];
        
        
    end
end
% p.At_ineq = blkdiag(p.A_ineq{:});
% p.bt_ineq = cat(1,p.b_ineq{:});
% p.At_eq = blkdiag(p.Aeq{:});
% p.bt_eq = cat(1,p.beq{:});
% p.Ht = blkdiag(p.H{:});

end
