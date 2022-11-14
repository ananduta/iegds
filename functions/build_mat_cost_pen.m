function p = build_mat_cost_pen(p)

for i = 1:p.n
%% Cost function 
    q = zeros(1,p.nx(i));
    q(p.id_dg) = p.en.q_dg(i);
    q(p.id_ch) = p.en.q_st(i);
    q(p.id_dh) = p.en.q_st(i);
    qp =q;
    q(p.id_mg) = p.en.d_mg;
    q(p.id_gu) = p.gn.d_gu;
    
    Q = diag(q);
    Qp = diag(qp);

%    Q = diag([p.en.q_dg(i) p.en.q_st(i) p.en.q_st(i) p.en.d_mg p.gn.d_gu 0 0 0 zeros(1,p.en.noN(i))]);
%    Qp= diag([p.en.q_dg(i) p.en.q_st(i) p.en.q_st(i) 0 0 0 0 0 zeros(1,p.en.noN(i))]);
    
    Qh{i} = p.m.Sx{i}'*kron(eye(p.h),Q)*p.m.Sx{i};
    Qph{i} = p.m.Sx{i}'*kron(eye(p.h),Qp)*p.m.Sx{i};
    
    p.m.Qh{i} = sparse(Qh{i});
    p.m.Qph{i} = sparse(Qph{i});
    
%    p.m.H{i} = 1/p.alpha(i)*eye(size(Qh{i},1)) + 2*Qh{i};
    p.m.H{i} = p.Alpha{i} + 2*Qh{i};
    
%    p.H{i} = p.A{i} + 2*Qh{i}; % Coefficient of the quadratic term
%    p.H_half{i} = sqrt(p.H{i});
%    p.H_half_inv{i} = inv(p.H_half{i});
    
    cc = 1;
    %c = [p.en.c_dg(i) p.en.c_st(i) p.en.c_st(i) 0 0 0 0 0 zeros(1,p.en.noN(i))]';
    c = zeros(p.nx(i),1);
    c(p.id_dg) = p.en.c_dg(i);
    c(p.id_ch) = p.en.c_st(i);
    c(p.id_dh) = p.en.c_st(i);
    c(p.id_mg) = p.en.d_mg_l;
    c(p.id_gu) = p.gn.d_gu_l;
%     for jj=1:p.en.noN(i)
%         j = p.en.N{i}(jj);
%         %if p.en.Adj(i,j) == 1
%         c(4+cc,1) = p.en.c_tr(i,j);
%         cc = cc+1;
%         %end
%     end
    p.m.ch{i} = p.m.Sx{i}'*kron(ones(p.h,1),c) ;
    
    % infinity norm penalty
    for jj = 1:p.gn.noN(i)
        j = p.gn.N{i}(jj);
        p.m.ch{i} = p.m.ch{i} + p.m.Spen{i,j}'*ones(p.h,1)*p.pen;
        p.m.ch{i} = p.m.ch{i} + p.m.Sgam{i,j}'*ones(p.h,1)*p.Gamma_pen;
    end
end
end