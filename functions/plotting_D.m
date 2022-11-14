% plotting sim D
% satisfying capacity constraints
%%

load('sim37b_D_ov_sce_3.mat')


P_line_h = zeros(np.h,np.b,np.b);
P_line_max = zeros(np.b);
h_idx_max = zeros(np.b);
for y=1:np.b
    for zz = 1:length(np.B{y})
        z = np.B{y}(zz);
%        P_line_h(:,y,z) = o.u_DSO{y}(np.h*(3+zz-1)+1:np.h*(3+zz),1);
        P_line_h(:,y,z) = s.p_l{y,z}(:,end);
        [P_line(y,z),h_idx_max(y,z)] = max(P_line_h(21,y,z));
        
    end
end

s_bar = np.s_bar;
%%
load('sim37b_D_fine_sce_3.mat')

P_line_o = zeros(np.h,np.b,np.b);
P_line_o_max = zeros(np.b);
h_idx_o_max = zeros(np.b);
for y=1:np.b
    for zz = 1:length(np.B{y})
        z = np.B{y}(zz);
        P_line_h(:,y,z) = s.p_l{y,z}(:,end);
        [P_line_o(y,z),h_idx_o_max(y,z)] = max(P_line_h(21,y,z));
    end
end

%%
c =1;
for i=1:np.b
    for jj = 1:length(np.B{i})
        j = np.B{i}(jj);
        if j > i
            Pl(c,:) = [i, j, abs(P_line(i,j)) abs(P_line_o(i,j))];
            c=c+1;
        end
    end
        
end