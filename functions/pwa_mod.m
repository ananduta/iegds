function pwa = pwa_mod(p,i)
% Linear constraints for PWA model of gas-flow equations
% W. Ananduta
% 04/10/2021


eps = 1e-8;


    
    
    
        
    
    %% Equality constraints


    % simplex of \delta^m (37): \sum_{m=1}^r \delta_(i,j)^m = 1
    for jj = 1:p.gn.noN(i)
        j=p.gn.N{i}(jj);
        
        
        
        E1 = zeros(1,p.nu(i));
        
        for m =1:p.r
            id_delta_m = p.nx(i)+p.ny(i)+(1+3*p.r)*(jj-1) + 1+3*(m-1) + 3;
            
            E1(1,id_delta_m) = 1;
            Aeq1j{jj} = sparse(kron(eye(p.h),E1));
            beq1j{j} = ones(p.h,1);
        end
    end
    if p.gn.noN(i) > 0
        Aeq1 = cat(1,Aeq1j{:});
        beq1 = cat(1,beq1j{:});
    else
        Aeq1 = [];
        beq1 = [];
    end
    
    pwa.Aeq = Aeq1;
    pwa.beq = beq1;

    %% Inequality constraints

    % (40)-(41),  (49)-(50)
    for jj = 1:p.gn.noN(i)
        j=p.gn.N{i}(jj);
        
        
        
        
        
        id_delta_psi = p.nx(i)+p.ny(i)+(1+3*p.gn.r)*(jj-1) + 1;
        
        % (40)-(41) : [\delta_(i,j)^\psi=1] <=> [\phi_(i,j) <= 0]
        %(40)
        A1d = zeros(1,p.nu(i));
        A1d(1,id_delta_psi)=p.gn.phi_max(i,j);
        
        A1j{jj} = kron(eye(p.h),A1d) + p.m.Sphi{i,j};
        b1j{jj} = p.gn.phi_max(i,j)*ones(p.h,1);
        
        %(41)
        A2d = zeros(1,p.nu(i));
        A2d(1,id_delta_psi)=-p.gn.phi_max(i,j)-eps;
        
        A2j{jj} = kron(eye(p.h),A2d) - p.m.Sphi{i,j};
        b2j{jj} = -eps*ones(p.h,1);
        
        
        id_z_psi = p.nx(i) + 2 + (2+p.gn.r)*(jj-1) + 2;
        
        
        % (49)-(50) : [y_(i,j)^{\psi_i} = \delta_(i,j)^\psi * \psi_i]
        %(49)
        A3h = zeros(1,p.nu(i));
        A3h(1,id_delta_psi) = p.gn.psi_min(i);
        A3h(1,id_z_psi) = -1;
        
        A3j{jj} = kron(eye(p.h),A3h);        
        b3j{jj} = 0*ones(p.h,1);
        
        A4h = zeros(1,p.nu(i));
        A4h(1,id_delta_psi) = -p.gn.psi_min(i);
        A4h(1,id_z_psi) = 1;
        
        A4j{jj} = kron(eye(p.h),A4h) - p.m.Spsi{i};        
        b4j{jj} = -p.gn.psi_min(i)*ones(p.h,1);
        
        
        %(50)
        A5h = zeros(1,p.nu(i));
        A5h(1,id_delta_psi) = -p.gn.psi_max(i);
        A5h(1,id_z_psi) = 1;
        
        A5j{jj} = kron(eye(p.h),A5h);        
        b5j{jj} = 0*ones(p.h,1);
        
        A6h = zeros(1,p.nu(i));
        A6h(1,id_delta_psi) = p.gn.psi_max(i);
        A6h(1,id_z_psi) = -1;
        
        A6j{jj} = kron(eye(p.h),A6h) + p.m.Spsi{i};        
        b6j{jj} = p.gn.psi_max(i)*ones(p.h,1);
        
        
    end
    if p.gn.noN(i) > 0
        A1 = cat(1,A1j{:});
        A2 = cat(1,A2j{:});
        A3 = cat(1,A3j{:});
        A4 = cat(1,A4j{:});
        A5 = cat(1,A5j{:});
        A6 = cat(1,A6j{:});

        b1 = cat(1,b1j{:});
        b2 = cat(1,b2j{:});
        b3 = cat(1,b3j{:});
        b4 = cat(1,b4j{:});
        b5 = cat(1,b5j{:});
        b6 = cat(1,b6j{:});
    else
        A1 = [];
        A2 = [];
        A3 = [];
        A4 = [];
        A5 = [];
        A6 = [];

        b1 = [];
        b2 = [];
        b3 = [];
        b4 = [];
        b5 = [];
        b6 = [];
    end
    
    clearvars('A1j','A2j','A3j','A4j','A5j','A6j');
    clearvars('b1j','b2j','b3j','b4j','b5j','b6j');
    
    %(42)-(48)
    for jj = 1:p.gn.noN(i)
        j=p.gn.N{i}(jj);
        
        
        % define the parameters of the affine function at each region
        nf = @(y) y^2/p.gn.cf(i,j);
        pwaf = pwa_approx_nf(p.r,-p.gn.phi_max(i,j),p.gn.phi_max(i,j),nf);
        pwa.f = pwaf;
    
        for m =1:p.r
            
            % indices
            id_delta_psi = p.nx(i)+p.ny(i)+(1+3*p.gn.r)*(jj-1) + 1;
            id_alpha_m = p.nx(i)+p.ny(i)+(1+3*p.gn.r)*(jj-1) + 1+3*(m-1) + 1;
            id_beta_m = p.nx(i)+p.ny(i) +(1+3*p.gn.r)*(jj-1) + 1+3*(m-1) + 2;
            id_delta_m = p.nx(i)+p.ny(i)+(1+3*p.gn.r)*(jj-1) + 1+3*(m-1) + 3;
            
            id_z_m = p.nx(i) + 2 + (2+p.gn.r)*(jj-1) + 2 + m;
            
                        
            %(47)-(48): [y_(i,j)^{m} = \delta_(i,j)^m * \phi_(i,j)]
            %(47)
            A7h = zeros(1,p.nu(i));
            A7h(1,id_delta_m) = -p.gn.phi_max(i,j);
            A7h(1,id_z_m) = -1;

            A7m{m} = kron(eye(p.h),A7h);        
            b7m{m} = 0*ones(p.h,1);

            A8h = zeros(1,p.nu(i));
            A8h(1,id_delta_m) = p.gn.phi_max(i,j);
            A8h(1,id_z_m) = 1;

            A8m{m} = kron(eye(p.h),A8h) - p.m.Sphi{i,j};        
            b8m{m} = p.gn.phi_max(i,j)*ones(p.h,1);


            %(48)
            A9h = zeros(1,p.nu(i));
            A9h(1,id_delta_m) = -p.gn.phi_max(i,j);
            A9h(1,id_z_m) = 1;

            A9m{m} = kron(eye(p.h),A9h);        
            b9m{m} = 0*ones(p.h,1);

            A10h = zeros(1,p.nu(i));
            A10h(1,id_delta_m) = p.gn.phi_max(i,j);
            A10h(1,id_z_m) = -1;

            A10m{m} = kron(eye(p.h),A10h) + p.m.Sphi{i,j};        
            b10m{m} = p.gn.phi_max(i)*ones(p.h,1);
            
            
            %(42)-(46): 
            %(42)
            A11h = zeros(1,p.nu(i));
            A11h(1,id_alpha_m) = p.gn.phi_max(i,j) - pwaf.M(m);
            b11h = p.gn.phi_max(i,j);
            
            A11m{m} = kron(eye(p.h),A11h) + p.m.Sphi{i,j};
            b11m{m} = b11h*ones(p.h,1);
            
            
            %(43)
            A12h = zeros(1,p.nu(i));
            A12h(1,id_alpha_m) = -p.gn.phi_max(i,j) - pwaf.M(m)-eps;
            b12h = - pwaf.M(m)-eps;
            
            A12m{m} = kron(eye(p.h),A12h) - p.m.Sphi{i,j};
            b12m{m} = b12h*ones(p.h,1);
            
            %(44)
            A13h = zeros(1,p.nu(i));
            A13h(1,id_beta_m) = p.gn.phi_max(i,j) + pwaf.m(m);
            b13h = p.gn.phi_max(i,j);
            
            A13m{m} = kron(eye(p.h),A13h) - p.m.Sphi{i,j};
            b13m{m} = b13h*ones(p.h,1);
            
            %(45)
            A14h = zeros(1,p.nu(i));
            A14h(1,id_beta_m) = -p.gn.phi_max(i,j) + pwaf.m(m)-eps;
            b14h = pwaf.m(m)-eps;
            
            A14m{m} = kron(eye(p.h),A14h) + p.m.Sphi{i,j};
            b14m{m} = b14h*ones(p.h,1);
            
            %(46)a
            A15h = zeros(1,p.nu(i));
            A15h(1,id_alpha_m) = -1;
            A15h(1,id_delta_m) = 1;
            
            A15m{m} = kron(eye(p.h),A15h);
            b15m{m} = 0*ones(p.h,1);
            
            %(46)b
            A16h = zeros(1,p.nu(i));
            A16h(1,id_beta_m) = -1;
            A16h(1,id_delta_m) = 1;
            
            A16m{m} = kron(eye(p.h),A16h);
            b16m{m} = 0*ones(p.h,1);
            
            %(46)c
            A16h = zeros(1,p.nu(i));
            A16h(1,id_alpha_m) = 1;
            A16h(1,id_beta_m) = 1;
            A16h(1,id_delta_m) = -1;
            
            A16m{m} = kron(eye(p.h),A16h);
            b16m{m} = 1*ones(p.h,1);
            
        end
        
        A716j{jj} = [cat(1,A7m{:});
                     cat(1,A8m{:});
                     cat(1,A9m{:});
                     cat(1,A10m{:});
                     cat(1,A11m{:});
                     cat(1,A12m{:});
                     cat(1,A13m{:});
                     cat(1,A14m{:});
                     cat(1,A15m{:});
                     cat(1,A16m{:})];
       b716j{jj} = [cat(1,b7m{:});
                     cat(1,b8m{:});
                     cat(1,b9m{:});
                     cat(1,b10m{:});
                     cat(1,b11m{:});
                     cat(1,b12m{:});
                     cat(1,b13m{:});
                     cat(1,b14m{:});
                     cat(1,b15m{:});
                     cat(1,b16m{:})];
       clearvars('A7m','A8m','A9m','A10m','A11m','A12m','A13m','A14m','A15m','A16m');
       clearvars('b7m','b8m','b9m','b10m','b11m','b12m','b13m','b14m','b15m','b16m')
        
    end
    if p.gn.noN(i) > 0 
        A716 = cat(1,A716j{:});
        b716 = cat(1,b716j{:});
    else
        A716 = [];
        b716 = [];
    end
    clearvars('A716j','b716j');
    
    % 0 <= z <= 1
    A17 = [ p.m.Sz{i};
           -p.m.Sz{i}];
    b17 = [ones(size(p.m.Sz{i},1),1);
           zeros(size(p.m.Sz{i},1),1)];
    
    pwa.Aineq   = [A1;
                   A2;
                   A3;
                   A4;
                   A5;
                   A6;
                   A716;
                   A17];
    
    pwa.bineq    = [b1;
                    b2;
                    b3;
                    b4;
                    b5;
                    b6;
                    b716;
                    b17];
    
%     pwa.Aineq   = [A1;
%                    A2;
%                    A3;
%                    A4;
%                    A716;
%                    A17];
%     
%     pwa.bineq    = [b1;
%                     b2;
%                     b3;
%                     b4;
%                     b716;
%                     b17];


end