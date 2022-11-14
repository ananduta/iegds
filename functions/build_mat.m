function np = build_mat(np)
% For prosumer of Extended P2P market
% W. Ananduta
% 23/09/2020

%[Sdg,Sst,np.Smg,np.Str] = gen_Smat(np);
% u_i = [p^di p^st p^mg {p^tr_j}]

for i=1:np.n
    Ni = sum(np.Adj(i,:));

    % Cost function
    Q = diag([np.q_dg(i) np.q_st(i) 0 zeros(1,Ni)]);
    Qh{i} = kron(Q,eye(np.h));
    np.Qh{i} = Qh{i};
    
    np.H{i} = np.A{i} + 2*Qh{i}; % Coefficient of the quadratic term
    
    cc = 1;
    c = [np.c_dg(i) np.c_st(i) 0 zeros(1,Ni)]';
    for j=1:np.n
        if np.Adj(i,j) == 1
            c(3+cc,1) = np.c_tr(i,j);
            cc = cc+1;
        end
    end
    np.ch{i} = kron(c,ones(np.h,1));
    
    % constraints U_i
    % Load power balance (equality constraint) (eq. (6))
    E1 = ones(1,3+Ni);
    np.Aeq{i} = kron(E1,eye(np.h));
    np.beq{i} = np.Pd(i,1:np.h)';

    % DG unit eq(3)
    A1 = [ 1 0 0 zeros(1,Ni);
          -1 0 0 zeros(1,Ni)];
    A1 = kron(A1,eye(np.h));  
    b1 = [ np.pdg_max(i)*ones(np.h,1);
          -np.pdg_min(i)*ones(np.h,1)];
    

    % Storage unit eq(5)
    A2 = [ Sst{i};
          -Sst{i}];
    b2 = [np.p_dh(i)*ones(np.h,1);
          np.p_ch(i)*ones(np.h,1)];

    if np.st_un(i) == 1
        % constructing the matrices for the dynamic equation (At and Bt)
        Am = np.a(i);
        Bm = [0 -1 0 zeros(1,Ni)]; %be careful when sampling time is not 1 hour

        Btcol = zeros(np.h,size(Bm,2));
        Bt = zeros(np.h,np.h*size(Bm,2));
        for l=1:np.h
            At(l,:) =Am^l;
            Btcol(l,:) = Am^(l-1)*Bm;
        end

        for l=1:np.h
            Bt(:,(l-1)*(length(Bm))+1:l*length(Bm)) = [zeros((l-1),length(Bm)); Btcol(1:size(Btcol,1)-(l-1),:)];
        end
        A3 = [ Bt; 
              -Bt];
        b3 = [np.x_max(i)*ones(np.h,1)-At*np.x0(i);
              -np.x_min(i)*ones(np.h,1)+At*np.x0(i)];
    else
        A3 = [];
        b3 = [];
    end


    % Power traded
    A4 = [zeros(Ni,3) eye(Ni)];
    A4 = kron(eye(np.h),A4);
    
    cc=1;
    pt_m = zeros(Ni,1);
    for j=1:np.n
        if np.Adj(i,j) ==1
            pt_m(cc,1) = np.pt_max(i,j);
            cc = cc+1;
        end
    end
    b4 = kron(ones(np.h,1),pt_m);
    
    np.A_ineq{i} = [A1;A2;A3;A4];
    np.b_ineq{i} = [b1;b2;b3;b4];
    
end
% np.At_ineq = blkdiag(np.A_ineq{:});
% np.bt_ineq = cat(1,np.b_ineq{:});
% np.At_eq = blkdiag(np.Aeq{:});
% np.bt_eq = cat(1,np.beq{:});
% np.Ht = blkdiag(np.H{:});

end
