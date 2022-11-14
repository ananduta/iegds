function s = loc_opt_c_qprog(np,s,t)
% Local optimization (using quadprog)
% W. Ananduta
% 02/12/2019

for i = 1:np.n
    
    % Construct coefficient of the linear cost
    
    % linear term associated with primal cost function
    cc = 1;
    c = [np.c_dg(i); np.c_st(i); 0; zeros(sum(np.Adj(i,:)),1)];
    for j=1:np.n
        if np.Adj(i,j) == 1
            c(3+cc,1) = np.c_tr(i,j);
            cc = cc+1;
        end
    end
    c1{i} = [];
    pmg_i_t = np.Smg{i}*s.u{i}(:,t);
    for h = 1:np.h
        c_1 = c + [0; 0; np.q_mg*(s.ph_mg(h,t)+pmg_i_t(h)); zeros(sum(np.Adj(i,:)),1)]; 
        c1{i} = [c1{i};c_1];
    end
    
    % linear term associated with grid constraints
    c2{i} = s.lambda{i}(:,t)'*[np.Smg{i};-np.Smg{i}];
    c2{i} = c2{i}';
    
    % linear term associated with reciprocity constraints
    c3{i} = zeros(np.h*(3+sum(np.Adj(i,:))),1);
    for j=1:np.n
        if np.Adj(i,j) == 1
            c3{i} = c3{i} + np.Str{i,j}'*s.mu{i,j}(:,t);
        end     
    end
    c3{i} = c3{i};
    
    % linear term associated with proximal term
    c4{i} = -2*s.u{i}(:,t);
    
    np.c{i} = c1{i} + c2{i} + c3{i} + c4{i};
    
    
    
end

A = np.At_ineq;
b = np.bt_ineq;
Aeq = np.At_eq;
beq = np.bt_eq;
H = np.Ht;
f = cat(1,np.c{:});

[u_all,fval,exitflag] = quadprog(H,f,A,b,Aeq,beq);

if exitflag ~= 1
%     disp(st.problem): in case solver cannot find the solution
    disp('not solved properly')
    st.problem
    pause
%     s.u{i}(:,t+1) = s.u{i}(:,t);
%     s.pmg{i}(:,t+1) = s.pmg{i}(:,t);
else
    %Assigning the decisions of agent i
    s.ph_mg(:,t+1) = 0;
    dim_prev_ag = 0;
    for i = 1:np.n
        dim_i = np.h*(3+sum(np.Adj(i,:)));
        u_i = u_all(dim_prev_ag+1:dim_prev_ag+dim_i);
        s.u{i}(:,t+1) = u_i;
        s.pmg{i}(:,t+1) = np.Smg{i}*u_i;
        s.ph_mg(:,t+1) = s.ph_mg(:,t+1) + s.pmg{i}(:,t+1);
        dim_prev_ag = dim_prev_ag+dim_i;
    end
end


end