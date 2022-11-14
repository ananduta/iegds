function s = loc_opt_qprog1(np,s,sl,t,i)
% Local optimization (using quadprog)
% W. Ananduta
% 04/12/2020

%for i = 1:np.n
    
    % Construct coefficient of the linear cost
    
    % linear term associated with primal cost function
    %cc = 1;
    c = [np.c_dg(i); np.c_st(i); 0; 0; zeros(sum(np.Adj(i,:)),1)];
    
    %mu_tr=[];
    for jj = 1:length(np.N{i})
        j = np.N{i}(jj);
        c(4+jj,1) = np.c_tr(i,j);
        
     %   mu_tr = [mu_tr;sl.mu_tr{i,j}(:,t)];
        
    end
    
    c1 = [];
    
    
    %pmg_i_t = np.Smg{i}*s.u{i}(:,t);
    for h = 1:np.h
        %c_1 = c + [0; 0; np.q_mg*(s.ph_mg(h,t)+pmg_i_t(h)); zeros(sum(np.Adj(i,:)),1)]; 
        c_1 = c + [0; 0; np.d_mg*(s.sigma_mg(h,t)); np.d_mgs*(s.sigma_mg(h,t)); zeros(sum(np.Adj(i,:)),1)]; 
        c1 = [c1;c_1];
    end
    
    % linear term associated with p_di
    c2 = np.Sdi{i}'*-sl.mu_pb{np.B_n(i)}(:,t);
    
    % linear term associated with p_st
    c3 = np.Sst{i}'*-sl.mu_pb{np.B_n(i)}(:,t);
    
    
    
    % linear term associated with grid constraints
    c4a = sl.lambda_mg(:,t)'*[np.Smg{i};-np.Smg{i}]+sl.mu_tg(:,t)'*np.Smg{i};
    c4a = c4a';
    
    % linear term associated with grid constraints
    c4b = sl.lambda_mg(:,t)'*[np.Smgs{i};-np.Smgs{i}]+sl.mu_tg(:,t)'*np.Smgs{i};
    c4b = c4b';
    
    % linear term associated with reciprocity constraints
    c5 = zeros(np.h*(4+sum(np.Adj(i,:))),1);
    for j=1:np.n
        if np.Adj(i,j) == 1
            c5 = c5 + np.Str{i,j}'*sl.mu_tr{i,j}(:,t);
        end     
    end
    %c3 = c3;
    
    % linear term associated with proximal term
    c6 = -np.A{i}*s.u{i}(:,t);
    
    np.c{i} = c1+ c2 + c3 + c4a + c4b + c5 + c6;
    
    
    
%end

A = np.A_ineq{i};
b = np.b_ineq{i};
Aeq = np.Aeq{i};
beq = np.beq{i};

%quadprog
H = np.H{i};
f = np.c{i};
options = optimset('Display','off');
%options = optimset('Display','on');
[u_i,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq,[],[],[],options);

%lsqlin
%C = np.H_half{i};
%d = -np.H_half_inv{i}*np.c{i};
%options = optimset('Display','off');
%[u_i,~,~,exitflag] = lsqlin(C,d,A,b,Aeq,beq,[],[],[],options);

if exitflag ~= 1 && exitflag ~= 0
%     disp(st.problem): in case solver cannot find the solution
    exitflag
    disp('not solved properly')
%    st.problem
    pause
%     s.u{i}(:,t+1) = s.u{i}(:,t);
%     s.pmg{i}(:,t+1) = s.pmg{i}(:,t);
else
    %Assigning the decisions of agent i
    %s.ph_mg(:,t+1) = 0;
    %dim_prev_ag = 0;
    %for i = 1:np.n
    %    dim_i = np.h*(3+sum(np.Adj(i,:)));
     %   u_i = u_all(dim_prev_ag+1:dim_prev_ag+dim_i);
        s.u{i}(:,t+1) = u_i;
        s.p_di{i}(:,t+1) = np.Sdi{i}*u_i;
        s.p_st{i}(:,t+1) = np.Sst{i}*u_i;
        s.p_mg{i}(:,t+1) = np.Smg{i}*u_i;
        s.p_mgs{i}(:,t+1) = np.Smgs{i}*u_i;
        for jj = 1:length(np.N{i})
            j = np.N{i}(jj);
            s.p_tr{i,j}(:,t+1) = np.Str{i,j}*u_i;
        end
    %    s.ph_mg(:,t+1) = s.ph_mg(:,t+1) + s.pmg{i}(:,t+1);
    %    dim_prev_ag = dim_prev_ag+dim_i;
    %end
end


end