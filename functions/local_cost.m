function Ji_now = local_cost(s,p,i,k)
% Compute local cost of prosumer i
% W. Ananduta
% 13/08/2021


Ji_now = s.u{i}(:,k)'*p.tn.Qh{i}*s.u{i}(:,k) + p.tn.ch{i}'*s.u{i}(:,k);

pmg_i_t = p.tn.Smg{i}*s.u{i}(:,k);
dgu_i_t = p.tn.Sgu{i}*s.u{i}(:,k);
c1 = [];
for h = 1:p.h
        
        c_1 = [0; 0; p.tn.d_mg*(s.sigma_mg(h,k)-pmg_i_t(h)); p.tn.d_gu*(s.sigma_gu(h,k)-dgu_i_t(h)); zeros(sum(p.tn.Adj(i,:)),1)]; 
        c1 = [c1;c_1];
end


Ji_now = Ji_now + c1'*s.u{i}(:,k);

end