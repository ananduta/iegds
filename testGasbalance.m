%%
% Identify dimensions of decision variables (per time step, h)
    for i = 1:p.n
        p.nx(i) = 8+ p.en.noN(i);
        p.nt(i) = p.gn.noN(i);
        p.ny(i) = 2+(2+p.gn.r)*p.gn.noN(i)+p.nt(i);
        p.nz(i) = (1+3*p.gn.r)*p.gn.noN(i);
        p.nu(i) = p.nx(i) + p.ny(i) + p.nz(i);
        
        p.id_phi{i} = 0;
        for jj = 1:p.gn.noN(i)
            p.id_phi{i}(jj) = p.nx(i) + 2+(2)*(jj-1)+1;
        end
        
    end
    
    
    p.id_dg = 1;
    p.id_ch = 2;
    p.id_dh = 3;
    p.id_mg = 4;
    p.id_gu = 5;
    p.id_th = 6;
    p.id_v = 7;
    p.id_eg = 8;
    
    % y vars. must +p.nx(i)
    p.id_psi = 1;
    p.id_gs = 2;
    
for i = 1:p.n
    
    % Gas balance eq. (12)
    E3 = zeros(1,p.nu(i));
    E3(1,p.id_gu) = -1;
    E3(1,p.nx(i)+2) = 1;
    for j = 1:p.gn.noN(i)
        E3(1,p.nx(i)+2+(2+p.gn.r)*(j-1)+1) = 1;
    end
    Aeq3 = sparse(kron(eye(p.h),E3));
    beq3 = p.gn.Gdem(i,1:p.h)';
    
    Gb(i) = Aeq3*o.u{i} - beq3;
end
max(Gb)
