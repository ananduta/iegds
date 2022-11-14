%%
dat{1} = 'sim_jext_22-Jan-2021_I_D3_1.mat';
dat{2} = 'sim_jext_22-Jan-2021_I_D3_3.mat';
dat{3} = 'sim_jext_22-Jan-2021_I_D3_5.mat';


load(dat{1})
o1 = o;
load(dat{2})
o2 = o;
load(dat{3})
o3 = o;
dat{6} = 'case_jext_21-Jan-2021_I_B.mat';
load(dat{6})
% compare dat 1 and 2
for i=1:np.n
    for j=1:np.n
        if np.Adj(i,j)==1
            er(i,j) = norm(o1.p_tr{i,j}-o2.p_tr{i,j},inf);
        else
            er(i,j) = 0;
        end
    end
    err1(i) = norm(er(i,:),inf);
end
norm(err1,inf)
figure

subplot(2,1,1)

bar(err1)
xlim([0 22])
title(['$\|p^{tr}_i - p^{tr}_i"\|_{\infty}$'], 'Interpreter','Latex')
% compare dat 1 and 2
for i=1:np.n
    for j=1:np.n
        if np.Adj(i,j)==1
            er(i,j) = norm(o1.p_tr{i,j}-o3.p_tr{i,j},inf);
            else
            er(i,j) = 0;
        end
    end
    err2(i) = norm(er(i,:),inf);
end
norm(err2,inf)
subplot(2,1,2)
bar(err2)
xlim([0 22])
