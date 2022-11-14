function m = gen_Smat_pen(p)
% Generate S^dg, S^st, S^tr, S^gu and S^mg matrices
% W. Ananduta
% 04/10/2021

% x
m.Sx = cell(p.n,1);
m.Sdg = cell(p.n,1);
m.Sch = cell(p.n,1);
m.Sdh = cell(p.n,1);
m.Smg = cell(p.n,1);
m.Sgu = cell(p.n,1);
m.Sth = cell(p.n,1);
m.Sv  = cell(p.n,1);
m.Spl = cell(p.n);

% y
m.Sy = cell(p.n,1);
m.Spsi = cell(p.n,1);
m.Sphi = cell(p.n);
m.Spen = cell(p.n);
% z
m.Sz = cell(p.n,1);

for i=1:p.n
    
    % compute Sx
    a = [eye(p.nx(i)) zeros(p.nx(i),p.ny(i)+p.nz(i))];
    m.Sx{i} = kron(eye(p.h),a);
    
    %compute a_ni,1
    a = zeros(p.nu(i),1);
    a(p.id_dg,1) = 1;
    
    %compute Sdg
    m.Sdg{i} = kron(eye(p.h),a');
    
    %compute a_ni,2
    a = zeros(p.nu(i),1);
    a(p.id_ch,1) = 1;
    
    %compute Sch
    m.Sch{i} = kron(eye(p.h),a');
    
    
    %compute a_ni,3
    a = zeros(p.nu(i),1);
    a(p.id_dh,1) = 1;
    
    %compute Sdh
    m.Sdh{i} = kron(eye(p.h),a'); 
    
    %compute a_ni,4
    a = zeros(p.nu(i),1);
    a(p.id_mg,1) = 1;
    
    %compute Smg
    m.Smg{i} = kron(eye(p.h),a');
    
    %compute a_ni,5
    a = zeros(p.nu(i),1);
    a(p.id_gu,1) = 1;
    
    %compute Sgu
    m.Sgu{i} = kron(eye(p.h),a');
    
    %compute a_ni,6
    a = zeros(p.nu(i),1);
    a(p.id_th,1) = 1;
    
    %compute Sth
    m.Sth{i} = kron(eye(p.h),a');
    
    %compute a_ni,7
    a = zeros(p.nu(i),1);
    a(p.id_v,1) = 1;
    
    %compute Sv
    m.Sv{i} = kron(eye(p.h),a');
    
    %compute a_ni,8
    a = zeros(p.nu(i),1);
    a(p.id_eg,1) = 1;
    
    %compute Sv
    m.Set{i} = kron(eye(p.h),a');
    
    c=1; %just a counter
    for j=1:p.n
        
        if p.en.Adj(i,j) ==1
            %compute a_ni,r(ij)
            a = zeros(p.nu(i),1);
            a(p.id_eg+c,1) = 1;
            
            %compute Spl
            m.Spl{i,j}=kron(eye(p.h),a');
            
%             %compute a_ni,r(ij)
%             a = zeros(p.nu(i),1);
%             a(7+p.en.noN(i)+c,1) = 1;
%             
%             %compute Sql
%             m.Sql{i,j}=kron(eye(p.h),a');
            
            c=c+1;
        end
        
    end
    
    
    % compute Sy
    a = [zeros(p.ny(i),p.nx(i)) eye(p.ny(i)) zeros(p.ny(i),p.nz(i))];
    m.Sy{i} = kron(eye(p.h),a);
    
    % compute Spsi
    a = zeros(1,p.nu(i));
    a(1,p.nx(i)+1) = 1;
    m.Spsi{i} = kron(eye(p.h),a);
        
    % compute gs
    a = zeros(1,p.nu(i));
    a(1,p.nx(i)+p.id_gs) = 1;
    m.Sgs{i} = kron(eye(p.h),a);

    % compute Sphi
    c=1;
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        a = zeros(1,p.nu(i));
        a(1,p.id_phi{i}(jj)) = 1;
        m.Sphi{i,j} = kron(eye(p.h),a);
        
        a = zeros(1,p.nu(i));
        a(1,p.id_phi{i}(jj)+1) = 1;
        m.Sgam{i,j} = kron(eye(p.h),a);
               
        a = zeros(1,p.nu(i));
        a(1,p.nx(i)+p.ny(i)-p.nt(i)+jj) = 1;
        m.Spen{i,j} = kron(eye(p.h),a);
    end
    
    % compute Sz
    a = [zeros(p.nz(i),p.nx(i)+p.ny(i)) eye(p.nz(i))];
    m.Sz{i} = kron(eye(p.h),a);
    
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        a = zeros(1,p.nu(i));
        a(1,p.nx(i)+p.ny(i)+jj) = 1;
        m.Salp{i,j} = kron(eye(p.h),a);
    end
    
end

end