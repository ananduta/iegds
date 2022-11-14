function p = build_mat_iegs_SOCPi(p)
% Building matrices for constraints with MISOCP gas flow relaxation

%np = p.en;

p.m = gen_Smat_pen(p);


% MISOCP gas flow model
misocp = misocp_gf(p);
p.m.Gc = misocp.Gc;
p.m.gc = misocp.gc;

% Other constraints
for i=1:p.n
    Ni = p.en.noN(i);

    
    
    %% Constraints U_i
    
    % Equality constraints
    
    % Electrical power balance (equality constraint) eq. (10a)
    %E1 = [ones(1,p.nx(i)) zeros(1,p.ny(i)+p.nz(i))];
    E1a =  zeros(1,p.nx(i)+p.ny(i)+p.nz(i));
    E1a(1,1:p.id_mg) = ones(1,p.id_mg);
    E1a(1,p.id_ch) = -1;
    Aeq1a = sparse(kron(eye(p.h),E1a));
    beq1a = p.en.Pd(i,1:p.h)';
    
    
    E1b = zeros(1,p.nx(i)+p.ny(i)+p.nz(i));
    E1b(1,p.id_mg) = 1;
    E1b(1,p.id_eg) = -1;
    E1b(1,p.id_eg+1:p.nx(i)) = 1*ones(1,p.en.noN(i));
    Aeq1b = sparse(kron(eye(p.h),E1b));
    beq1b = zeros(p.h,1);
    
    Aeq1 = [Aeq1a;Aeq1b];
    beq1 = [beq1a;beq1b];
    
    % p_et for not-connected nodes (17)
    if p.en.et(i) == 0
        E5 = zeros(1,p.nu(i));
        E5(1,p.id_eg) = 1;
        Aeq5 = sparse(kron(eye(p.h),E5));
        beq5 = zeros(p.h,1);            
    else
        Aeq5 = [];
        beq5 = [];
    end
    
    
    Aeq3 = [];
    beq3 = []; 
    
    % Gas balance eq. (12)
    E3 = zeros(1,p.nu(i));
    E3(1,p.id_gu) = -1;
    E3(1,p.nx(i)+2) = 1;
    for j = 1:p.gn.noN(i)
        E3(1,p.id_phi{i}(j)) = 1;
    end
    Aeq3 = sparse(kron(eye(p.h),E3));
    beq3 = p.gn.Gdem(i,1:p.h)'; 
    
    
    % Coupling between gas load and power dispatchable unit (3)
    if p.en.dgu_un(i) == 1
        E2 = zeros(1,p.nu(i));
        E2(1,p.id_dg) = -p.en.q2(i);
        E2(1,p.id_gu) = 1;
        Aeq2 = sparse(kron(eye(p.h),E2));
        beq2 = p.en.q3(i)*ones(p.h,1);
    else
        % otherwise, set decision of d_gu = 0
        E2 = zeros(1,p.nu(i));
        E2(1,p.id_gu) = 1;
        Aeq2 = sparse(kron(eye(p.h),E2));
        beq2 = zeros(p.h,1);
    end
    
    
    Aeq4 = [];
    beq4 = [];
% g_t for not-connected nodes (17)
    if p.gn.gt(i) == 0
        E4 = zeros(1,p.nu(i));
        E4(1,p.nx(i)+p.id_gs) = 1;
        Aeq4 = sparse(kron(eye(p.h),E4));
        beq4 = zeros(p.h,1);
    else
        Aeq4 = [];
        beq4 = [];
    end
    

    p.m.Aeq{i} = sparse([Aeq1;Aeq2;Aeq3;Aeq4;Aeq5]);
    p.m.beq{i} =        [beq1;beq2;beq3;beq4;beq5];
    
    % Inequality constraints
    
    % DG unit eq(1)
    A1 = [ p.m.Sdg{i};
          -p.m.Sdg{i}];
    b1 = [ p.en.pdg_max(i)*ones(p.h,1);
          -p.en.pdg_min(i)*ones(p.h,1)];
    

    % Storage unit eq(5)                                                    
    A2a = [ p.m.Sch{i};
           -p.m.Sch{i}];
    b2a = [p.en.p_ch(i)*ones(p.h,1);
          0*ones(p.h,1)];
      
    A2b = [ p.m.Sdh{i};
           -p.m.Sdh{i}];
    b2b = [p.en.p_dh(i)*ones(p.h,1);
           0*ones(p.h,1)];
    
    A2 = [A2a;A2b];
    b2 = [b2a;b2b];  
    
    if p.en.st_un(i) == 1
        % constructing the matrices for the dynamic equation (At and Bt)
        Am = p.en.a_st(i);
        Bm = zeros(1,p.nu(i));
        Bm(1,p.id_ch) = p.en.b_st(i)*p.en.eta_ch(i);
        Bm(1,p.id_dh) = -p.en.b_st(i)/p.en.eta_dh(i);
        
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

    A5 = [];
    b5 = [];
    A6 = [];
    b6 = [];
%     % Import from main grid 
%     A5 = [-p.m.Smg{i}];
%     b5 = [zeros(p.h,1)];
%     
    % gas demand 
    A6 = [-p.m.Sgu{i}];
    b6 = [zeros(p.h,1)];
    
%     % Power traded eq. (8)
%     A4 = [zeros(Ni,4) eye(Ni) zeros(Ni,p.ny(i)+p.nz(i))];
%     A4 = kron(eye(p.h),A4);
%     
%     pt_m = zeros(Ni,1);
%     for jj=1:p.en.noN(i)
%         j = p.en.N{i}(jj);
%         pt_m(jj,1) = p.en.pt_max(i,j);
%     end
%     b4 = kron(ones(p.h,1),pt_m);
    
    % Physical constraints of electrical network  
    
    % Voltage angles
    A4a = [ p.m.Sth{i};
          -p.m.Sth{i}];
    b4a = [ p.en.theta_max(i)*ones(p.h,1);
           -p.en.theta_min(i)*ones(p.h,1)];
    
    % Voltage magnitudes
    A4b = [ p.m.Sv{i};
          -p.m.Sv{i}];
    b4b = [ p.en.v_max(i)*ones(p.h,1);
          -p.en.v_min(i)*ones(p.h,1)];
      
    A4 = [A4a];
    b4 = [b4a];
    


    
    % flow limits (15)
    % and non-negativity of s ============ CHECK THIS PART
    % and non-negativity of gamma ======== CHECK THIS PART
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        
        A7j{jj} = [ p.m.Sphi{i,j};
                   -p.m.Sphi{i,j}];
        b7j{jj} = [ p.gn.phi_max(i,j)*ones(p.h,1);
                    p.gn.phi_max(i,j)*ones(p.h,1)];
        
        Apj{jj} = [-p.m.Spen{i,j}-p.m.Sphi{i,j};
                   -p.m.Spen{i,j}+p.m.Sphi{i,j}];
        bpj{jj} = [zeros(p.h,1);
                   zeros(p.h,1)];
               
%         Apj{jj} = -p.m.Spen{i,j};
%         bpj{jj} = zeros(p.h,1);
               
        A6j{jj} = -p.m.Sgam{i,j};
        b6j{jj} = zeros(p.h,1);
    end

    if p.gn.noN(i) > 0
        A7 = [];
        b7 = [];
       A7 = cat(1,A7j{:});
       b7 = cat(1,b7j{:});
        
        Ap = cat(1,Apj{:});
        bp = cat(1,bpj{:});
        
        A6 = cat(1,A6j{:});
        b6 = cat(1,b6j{:});
    else
        A7 = [];
        b7 = [];
        Ap = [];
        bp = [];
        A6 = [];
        b6 = [];
    end   
    
    clearvars('A7j','b7j','Apj','bpj','A6j','b6j');
    
    % pressure limits (16)
    A8 = [];
    b8 = [];
    A8 = [ p.m.Spsi{i};
          -p.m.Spsi{i}];
    b8 = [ p.gn.psi_max(i)*ones(p.h,1);
          -p.gn.psi_min(i)*ones(p.h,1)];  
    
    
    % SCP constraint: linearization of concave constraint
    if p.scpFlag == 1
        for jj = 1:p.gn.noN(i)
            j = p.gn.N{i}(jj);
            
            A9j{jj} = p.m.Sgam{i,j} - p.m.Spen{i,j} - 2*diag(p.flow{i,j})/p.gn.cf(i,j)*p.m.Sphi{i,j};
            b9j{jj} = -p.flow{i,j}.^2/p.gn.cf(i,j);
            
        end
        
        if p.gn.noN(i) > 0
            A9 = [];
            b9 = [];
            A9 = cat(1,A9j{:});
            b9 = cat(1,b9j{:});
        else
            A9 = [];
            b9 = [];
        end
        clearvars('A9j','b9j')
    else
        A9 = [];
        b9 = [];
    end  
      
      
         
    p.m.Aineq{i} = sparse([A1;A2;A3;A4;A5;A6;A7;A8;Ap;A9]);
    p.m.bineq{i} =       ([b1;b2;b3;b4;b5;b6;b7;b8;bp;b9]);
    
    
    % combine all local constraints
    p.m.Aineq{i} = [p.m.Aineq{i}; misocp.Al{i}];
    p.m.bineq{i} = [p.m.bineq{i}; misocp.bl{i}];
   
    
    % Coupling constraints of power flow model                             
    for jj = 1:p.en.noN(i)
        j = p.en.N{i}(jj);
        
        % power flow equation
        PFi = zeros(1,p.nu(i));
        PFi(1,p.id_th) = p.en.Bnet(i,j);
        PFi(1,p.id_v) = -p.en.Gnet(i,j);
        PFi(1,p.id_th+jj) = -1;
        
        PFj = zeros(1,p.nu(j));
        PFj(1,p.id_th) = -p.en.Bnet(i,j);
        PFj(1,p.id_v) = p.en.Gnet(i,j);
        
        p.m.PF{i,j}{1} = kron(eye(p.h),PFi);
        p.m.PF{i,j}{2} = kron(eye(p.h),PFj);
        p.m.pf{i,j} = zeros(p.h,1);
        
    end
    
    
    
end


end
