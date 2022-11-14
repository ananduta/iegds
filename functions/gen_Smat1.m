function [Sdg,Sst,Smg,Smgs,Str] = gen_Smat1(np)
% Generate S^tr and S^mg matrices
% W. Ananduta
% 24/07/2019

Sdg = cell(np.n,1);
Sst = cell(np.n,1);
Smg = cell(np.n,1);
Smgs = cell(np.n,1);
Str = cell(np.n);
for i=1:np.n
    Ni = sum(np.Adj(i,:));
    
    %compute a_ni,1
    a = zeros(4+Ni,1);
    a(1,1) = 1;
    
    %compute Sdg
    Sdg{i} = kron(eye(np.h),a');
    
    %compute a_ni,2
    a = zeros(4+Ni,1);
    a(2,1) = 1;
    
    %compute Sst
    Sst{i} = kron(eye(np.h),a');
    
    
    %compute a_ni,3
    a = zeros(4+Ni,1);
    a(3,1) = 1;
    
    %compute Smg
    Smg{i} = kron(eye(np.h),a'); 
    
    %compute a_ni,4
    a = zeros(4+Ni,1);
    a(4,1) = 1;
    
    %compute Smgs
    Smgs{i} = kron(eye(np.h),a'); 
    
    c=1; %just a counter
    for j=1:np.n
        
        if np.Adj(i,j) ==1
            %compute a_ni,r(ij)
            a = zeros(4+Ni,1);
            a(3+c,1) = 1;
            
            %compute Str
            Str{i,j}=kron(eye(np.h),a');
            
            c=c+1;
        end
        
    end

end

end