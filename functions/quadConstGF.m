function [c,ceq] = quadConstGF(x,p)
    countQC = 1;
    idStart = p.idAgent(:,1);
    idEnd = p.idAgent(:,2);

    for i = 1:p.n
        xi = x(idStart(i):idEnd(i),1);
        %if ~isempty(p.gn.N{i})
            for jj = 1:length(p.gn.N{i})
                j = p.gn.N{i}(jj);
                

                for h = 1:p.h
                    id_gamma =  p.nx(i) + 2+(2)*jj;
%                     id_gamma_ijh = idStart(i)-1+p.nu(i)*(h-1)+p.nx(i) + 2+(2)*jj;
%                     id_phi_ijh = idStart(i)-1+p.nu(i)*(h-1)+ p.id_phi{i}(jj);
%                     
%                     Qc = (zeros(idEnd(p.n)));
%                     Qc(id_phi_ijh,id_phi_ijh) = 1/p.gn.cf(i,j);
%     
%                     q = (zeros(idEnd(p.n),1));
%                     q(id_gamma_ijh) = -1;
%                     rhs = 0;
%                     c(countQC) = x'*Qc*x + q'*x - rhs;

                    id_gamma_ijh = p.nu(i)*(h-1)+p.nx(i) + 2+(2)*jj;
                    id_phi_ijh = p.nu(i)*(h-1)+ p.id_phi{i}(jj);

                    Qc = zeros(idEnd(i)-idStart(i)+1);
                    Qc(id_phi_ijh,id_phi_ijh) = 1/p.gn.cf(i,j);

                    q = (zeros(idEnd(i)-idStart(i)+1,1));
                    q(id_gamma_ijh) = -1;

                    rhs = 0;

                    c(countQC) = xi'*Qc*xi + q'*xi - rhs;

                    countQC = countQC+1;
                end

            end
        %end
    end
    ceq = [];

end