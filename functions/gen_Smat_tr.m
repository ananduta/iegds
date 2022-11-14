function [Sdg,Sst,Smg,Str] = gen_Smat_tr(np)
% Generate S^tr and S^mg matrices
% W. Ananduta
% 24/07/2019


Smg = cell(np.n,1);
Str = cell(np.n);
Str_ax= cell(np.n);
for i=1:np.n
    Ni = sum(np.Adj(i,:));
    
    %compute a_ni,1
    a = zeros(3+Ni*2,1);
    a(1,1) = 1;
    
    %compute Sdg
    Sdg{i} = kron(eye(np.h),a');
    
    %compute a_ni,2
    a = zeros(3+Ni*2,1);
    a(2,1) = 1;
    
    %compute Sst
    Sst{i} = kron(eye(np.h),a');
    
    
    %compute a_ni,3
    a = zeros(3+Ni*2,1);
    a(3,1) = 1;
    
    %compute Smg
    Smg{i} = kron(eye(np.h),a'); 
    
    c=1; %just a counter
    for j=1:np.n
        
        if np.Adj(i,j) ==1
            %compute a_ni,r(ij)
            a = zeros(3+Ni*2,1);
            a1 = zeros(3+Ni*2,1);
            a(3+c,1) = 1;
            a1(3+Ni+c,1) = 1;
            
            %compute Str
            Str{i,j}=kron(eye(np.h),a');
            %Str_ax{i,j}=kron(eye(np.h),a1');
            c=c+1;
        end
        
    end

end

end