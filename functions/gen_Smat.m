function m = gen_Smat(p)
% Generate S^dg, S^st, S^tr, S^gu and S^mg matrices
% W. Ananduta
% 04/10/2021

% x
m.Sx = cell(p.n,1);
m.Sdg = cell(p.n,1);
m.Sst = cell(p.n,1);
m.Smg = cell(p.n,1);
m.Sgu = cell(p.n,1);
m.Str = cell(p.n);

% y
m.Sy = cell(p.n,1);
m.Spsi = cell(p.n,1);
m.Sphi = cell(p.n);

% z
m.Sz = cell(p.n,1);

for i=1:p.n
    
    % compute Sx
    a = [eye(p.nx(i)) zeros(p.nx(i),p.ny(i)+p.nz(i))];
    m.Sx{i} = kron(eye(p.h),a);
    
    %compute a_ni,1
    a = zeros(p.nu(i),1);
    a(1,1) = 1;
    
    %compute Sdg
    m.Sdg{i} = kron(eye(p.h),a');
    
    %compute a_ni,2
    a = zeros(p.nu(i),1);
    a(2,1) = 1;
    
    %compute Sst
    m.Sst{i} = kron(eye(p.h),a');
    
    
    %compute a_ni,3
    a = zeros(p.nu(i),1);
    a(3,1) = 1;
    
    %compute Smg
    m.Smg{i} = kron(eye(p.h),a'); 
    
    %compute a_ni,4
    a = zeros(p.nu(i),1);
    a(4,1) = 1;
    
    %compute Sgu
    m.Sgu{i} = kron(eye(p.h),a');
    
    c=1; %just a counter
    for j=1:p.n
        
        if p.en.Adj(i,j) ==1
            %compute a_ni,r(ij)
            a = zeros(p.nu(i),1);
            a(4+c,1) = 1;
            
            %compute Str
            m.Str{i,j}=kron(eye(p.h),a');
            
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
    
    % compute Sphi
    c=1;
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        a = zeros(1,p.nu(i));
        a(1,p.nx(i)+2+(2+p.r)*(jj-1)+1) = 1;
        m.Sphi{i,j} = kron(eye(p.h),a);
    end
    
    % compute Sz
    a = [zeros(p.nz(i),p.nx(i)+p.ny(i)) eye(p.nz(i))];
    m.Sz{i} = kron(eye(p.h),a);
    
end

end