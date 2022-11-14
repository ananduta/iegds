 function [Ji,Pi] = cost_compute(o,p,i)

    % compute aggregative terms
    o.sigma_mg = p.en.sumPd(1:p.h);
    o.sigma_gu = p.gn.sumGd(1:p.h);
    for j=1:p.n
        o.sigma_mg = o.sigma_mg + o.p_mg{j};
        o.sigma_gu = o.sigma_gu + o.d_gu{j};
    end


    xi = p.m.Sx{i}*o.u{i};
    
    c1 = [];
%     cc = 1;
%     c = [p.en.c_dg(i) p.en.c_st(i) p.en.c_st(i) p.en.d_mg_l p.gn.d_gu_l 0 0 0 zeros(1,p.en.noN(i))]';
%     for jj=1:p.en.noN(i)
%         j = p.en.N{i}(jj);
%         c(4+cc,1) = p.en.c_tr(i,j);
%         cc = cc+1;
%     end
    c = zeros(p.nx(i),1);
    c(p.id_dg) = p.en.c_dg(i);
    c(p.id_ch) = p.en.c_st(i);
    c(p.id_dh) = p.en.c_st(i);
    c(p.id_mg) = p.en.d_mg_l;
    c(p.id_gu) = p.gn.d_gu_l;
    
    for h = 1:p.h
        c_1 = c;
        c_1(p.id_mg) = p.en.d_mg*(o.sigma_mg(h)-o.p_mg{i}(h));
        c_1(p.id_gu) = p.gn.d_gu*(o.sigma_gu(h)-o.d_gu{i}(h));
        c1 = [c1;c_1];
    end
    
    Ji = o.u{i}'*p.m.Qh{i}*o.u{i} + c1'*xi;
    
    floc = o.u{i}'*p.m.Qph{i}*o.u{i} + kron(ones(p.h,1),c)'*xi;
    
    Di = p.m.Smg{i}'*p.en.d_mg*eye(p.h)*p.m.Smg{i} + p.m.Sgu{i}'*p.gn.d_gu*eye(p.h)*p.m.Sgu{i};
    
    Pi = 0.5*(Ji + floc + o.u{i}'*Di*o.u{i});
end