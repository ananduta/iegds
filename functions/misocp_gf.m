function misocp = misocp_gf(p)
% build matrices for linear coupling constraints of mixed-integer SOCP reformulation
% of the gas-flow equations
% W. Ananduta
% 16/02/2022
for i = 1:p.n
	for jj = 1:p.gn.noN(i)
        j=p.gn.N{i}(jj);
        
        id_psi = p.nx(i) + p.id_psi;
        id_gamma =  p.nx(i) + 2 +2*jj;
        id_alpha =  p.nx(i) + p.ny(i) + jj;
        
        
        id_psi_j = p.nx(j) + p.id_psi;
        
        %(33)
        E1i = zeros(1,p.nu(i));
        E1i(1,id_gamma) = -1;
        E1i(1,id_psi) = 1;
        E1i(1,id_alpha) = 2*(p.gn.psi_min(j) - p.gn.psi_max(i));
        
        
        E1j = zeros(1,p.nu(j));
        E1j(1,id_psi_j) = -1;
        
        f1 = 0;
    
        
        %(34)
        E2i = zeros(1,p.nu(i));
        E2i(1,id_gamma) = -1;
        E2i(1,id_psi) = -1;
        E2i(1,id_alpha) = 2*(p.gn.psi_max(j) - p.gn.psi_min(i));
        
        
        E2j = zeros(1,p.nu(j));
        E2j(1,id_psi_j) = 1;
        
        f2 = 2*(p.gn.psi_max(j) - p.gn.psi_min(i));
        
        %(35)
        E3i = zeros(1,p.nu(i));
        E3i(1,id_gamma) = 1;
        E3i(1,id_psi) = -1;
        E3i(1,id_alpha) = -2*(p.gn.psi_max(j) - p.gn.psi_min(i));
        
        
        E3j = zeros(1,p.nu(j));
        E3j(1,id_psi_j) = 1;
        
        f3 = 0;
        
        %(36)
        E4i = zeros(1,p.nu(i));
        E4i(1,id_gamma) = 1;
        E4i(1,id_psi) = 1;
        E4i(1,id_alpha) = -2*(p.gn.psi_min(j) - p.gn.psi_max(i));
        
        
        E4j = zeros(1,p.nu(j));
        E4j(1,id_psi_j) = -1;
        
        f4 = -2*(p.gn.psi_min(j) - p.gn.psi_max(i));
        
        Eci = [E1i;E2i;E3i;E4i];
        Ecj = [E1j;E2j;E3j;E4j];
        f = [f1;f2;f3;f4];
        
        Ec{i,j}{1} = kron(eye(p.h),Eci);
        Ec{i,j}{2} = kron(eye(p.h),Ecj);
        fc{i,j} = kron(ones(p.h,1),f);
        
        El1 = [];
        fl1 = [];
        El2 = [];
        fl2 = [];
        El3 = [];
        fl3 = [];
        El4 = [];
        fl4 = [];
        
        
        % (29)
        
        El1 = zeros(1,p.nu(i));
        El1(1,p.id_phi{i}(jj)) = -1;
        El1(1,id_alpha) = p.gn.phi_max(i,j);
        fl1 = p.gn.phi_max(i,j);
        
        El2 = zeros(1,p.nu(i));
        El2(1,p.id_phi{i}(jj)) = 1;
        El2(1,id_alpha) = -p.gn.phi_max(i,j);
        fl2 = 0;
               
        
        % bounds of the binary variables (0,1) (in case they become
        % continuous
        El3 = zeros(1,p.nu(i));
        El3(1,id_alpha) = 1;
        fl3 = 1;
        
        El4 = zeros(1,p.nu(i));
        El4(1,id_alpha) = -1;
        fl4 = 0;
        
    
        
        El = [El1;El2;El3;El4];
        fl = [fl1;fl2;fl3;fl4];

        Elj{jj} = kron(eye(p.h),El);
        flj{jj} = kron(ones(p.h,1),fl);
    end
    if p.gn.noN(i) > 0
        misocp.Al{i} = cat(1,Elj{:});
        misocp.bl{i} = cat(1,flj{:});
    else
        misocp.Al{i} = [];
        misocp.bl{i} = [];
    end
    
    clearvars('Elj','flj');
end

misocp.Gc = Ec;
misocp.gc = fc;

end