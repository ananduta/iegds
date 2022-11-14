function Jno_now = local_cost_NO(s,p,k)
% Compute local cost of network operator
% W. Ananduta
% 13/08/2021

J = 0;
for i = 1:p.gn.n
    q = kron(ones(1,p.h),[0;ones(p.gn.noN(i),1)]);
    Q = diag(q);
    J = J + s.u_no{i}(:,k)'*Q*s.u_no{i}(:,k);
end
Jno_now = J;
end