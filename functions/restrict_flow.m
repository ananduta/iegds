function [phi_max_a,flag] = restrict_flow(gn)
% Computing additional flow constraints (restricting flows)
% To guarantee finding a feasible solution
% W. Ananduta
% 02/11/2021


psi_max = min(gn.psi_max);
psi_min = max(gn.psi_min);
dpsi = 1/4*(psi_max - psi_min);







% % uniform parameters of the gas network
% cf = gn.cf(i,j);
% r = gn.r;

% non-uniform cf while others are uniform
cf = gn.cf(:);
cf = cf(cf~=0);
cf = min(cf);

% for i = gn.n
%     for jj = gn.noN(i)
%         j = gn.N{i}(jj);
%         
%         if gn.cf(i,j) == cf
%             ag_i = i;
%             ag_j = j;
%         end
%     end
% end
% 
% i = ag_i;
% j = ag_j;

i = 1;
j = gn.N{1}(1);

% check if restriction is needed
if dpsi >= gn.phi_max(i,j)^2/cf
    phi_max_a = phi_max;
    flag = 0;
else
    
    % compute restriction
    flag = 1;
    nf = @(y) y^2/cf;
    pwaf = pwa_approx_nf(gn.r,-gn.phi_max(i,j),gn.phi_max(i,j),nf);

    for m = 1:gn.r
        if dpsi >= pwaf.a(m)*pwaf.l(m) + pwaf.b(m) && dpsi <= pwaf.a(m)*pwaf.l(m+1) + pwaf.b(m)
            phi_max_a = 1/pwaf.a(m)*(dpsi - pwaf.b(m));
        end
    end
end

end