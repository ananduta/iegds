function s = projUy(s,a_DSO,np,y,t)
% Projection onto local set of bus y
% W. Ananduta
% 13/11/2020

if t==0
    x_aux = zeros(np.h*(3+2*length(np.B{y})),1);
else
    x_aux = s.u_DSO{y}(:,t)- np.A_DSO{y}*a_DSO;   
end

z_aux = zeros(np.h*(3+2*length(np.B{y})),1);

% project theta (18)
z_aux(1:np.h,1) = max(min(x_aux(1:np.h,1),np.theta_max(y)),np.theta_min(y));

% project v (19)
z_aux(np.h+1:2*np.h,1) = max(min(x_aux(np.h+1:2*np.h,1),np.v_max(y)),np.v_min(y));

% project p_tg (20)
if ~isempty(find(np.B_mg==y))
    z_aux(2*np.h+1:3*np.h,1) = x_aux(2*np.h+1:3*np.h,1);
end

% project p_l and q_l (17)
for zz = 1:length(np.B{y})
    z = np.B{y}(zz);
    p_l = x_aux(np.h*(3+zz-1)+1:np.h*(3+zz),1);
    q_l = x_aux(np.h*(3+length(np.B{y})+zz-1)+1:np.h*(3+length(np.B{y})+zz),1);

    for h = 1:np.h
        c_aux = np.s_bar(y,z)/max(norm([p_l(h);q_l(h)],2),np.s_bar(y,z));
        z_aux(np.h*(3+zz-1)+h,1) = c_aux*p_l(h);
        z_aux(np.h*(3+length(np.B{y})+zz-1)+h,1) = c_aux*q_l(h);
    end
end
        


s.u_DSO{y}(:,t+1) = z_aux;
s.theta{y}(:,t+1) = s.u_DSO{y}(1:np.h,t+1);
s.v{y}(:,t+1) = s.u_DSO{y}(np.h+1:np.h*2,t+1);
s.p_tg{y}(:,t+1) = s.u_DSO{y}(np.h*2+1:np.h*3,t+1);
    for zz = 1:length(np.B{y})
        z = np.B{y}(zz);
        s.p_l{y,z}(:,t+1) =  s.u_DSO{y}(np.h*(3+zz-1)+1:np.h*(3+zz),t+1);
        s.q_l{y,z}(:,t+1) =  s.u_DSO{y}(np.h*(3+length(np.B{y})+zz-1)+1:np.h*(3+length(np.B{y})+zz),t+1);
    end
end



