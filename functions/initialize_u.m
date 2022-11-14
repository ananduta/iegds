function  s = initialize_u(np,s,i)
% Local optimization
% W. Ananduta
% 24/07/2019

% Inputs:
% np: Network parameters
% s: simulation results
% i: agent
% t: iteration

yalmip('clear');

%initialize cost function and constraints
cons = [];
Jc = 0;

Ni = sum(np.Adj(i,:)); % number of neighbors

z = sdpvar(np.h*(3+Ni),1); % initialize decisions

x0 = np.x0(i); % state of charge of storage


% Cost function

% constraints U_i

% % Load power balance eq(1) (E1t constructed in build_mat.m)
% % E1 = ones(1,3+Ni);
% % E1t = kron(eye(np.h),E1);
 E1t = np.E1t{i};
cons = [cons,  E1t*z ==  np.Pd(i,1:np.h)'];
%cons = [cons, -E1t*z <= -np.Pd(i,1:np.h)'];

% DG unit eq(3) (E2t constructed in build_mat.m)
% E2 = [1 0 0 zeros(1,Ni)];
% E2t = kron(eye(np.h),E2);
E2t = np.E2t{i};
cons = [cons,  E2t*z <=  np.pdg_max(i)*ones(np.h,1)];
cons = [cons, -E2t*z <= -np.pdg_min(i)*ones(np.h,1)];

% Storage unit eq(7)-(8) (E3t constructed in build_mat.m)
% E3 = [0 1 0 zeros(1,Ni)];
% E3t = kron(eye(np.h),E3);
E3t= np.E3t{i};
cons = [cons,  E3t*z <= np.p_dh(i)*ones(np.h,1)];
cons = [cons, -E3t*z <= np.p_ch(i)*ones(np.h,1)];

% Storage unit eq(5)-(6) (At and Bt constructed in build_mat.m)
 if np.st_un(i) == 1
%     % constructing the matrices for the dynamic equation (At and Bt)
%     A = np.a(i);
%     B = [0 -1 0 zeros(1,Ni)]; %be careful when sampling time is not 1 hour
% 
%     Btcol = zeros(np.h,size(B,2));
%     for l=1:np.h
%         At(l,:) =A^l;
%         Btcol(l,:) = A^(l-1)*B;
%     end
% 
%     for l=1:np.h
%         Bt(:,(l-1)*(length(B))+1:l*length(B)) = [zeros((l-1),length(B)); Btcol(1:size(Btcol,1)-(l-1),:)];
%     end
    At =np.At{i};
    Bt =np.Bt{i};
    cons = [cons, -At*x0+Bt*z <= -np.x_min(i)*ones(np.h,1)]; 
    cons = [cons,  At*x0+Bt*z <=  np.x_max(i)*ones(np.h,1)];
end


% Power traded eq(10) (E4t and Pt_m constructed in build_mat.m)
% cc=1;
% for j=1:np.n
%     if np.Adj(i,j) ==1
%         pt_m(cc,1) = np.pt_max(i,j);
%         cc = cc+1;
%     end
% end
% Pt_m = kron(ones(np.h,1),pt_m);
% E4 = [zeros(Ni,3) eye(Ni)];
% E4t = kron(eye(np.h),E4);
E4t = np.E4t{i};
Pt_m = np.Pt_m{i};
cons = [cons, E4t*z <= Pt_m];
cons = [cons,-E4t*z <= Pt_m];

% solve the optimization problem
ops = sdpsettings('verbose',0,'solver','quadprog');
st=optimize(cons,Jc,ops);

if st.problem ~= 0
%     disp(st.problem): in case solver cannot find the solution
    disp('not solved properly')
    st.problem
    pause
%     s.u{i}(:,t+1) = s.u{i}(:,t);
%     s.pmg{i}(:,t+1) = s.pmg{i}(:,t);
else
    %Assigning the decisions of agent i
    u_i = value(z);
    s.u{i}(:,1) = u_i;
    s.pmg{i}(:,1) = np.Smg{i}*u_i;
end

end

