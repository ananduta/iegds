%%
dat{1} = 'sim_jext_26-Jan-2021_II_D1_1.mat';
dat{2} = 'sim_jext_26-Jan-2021_II_D1_2.mat';
dat{3} = 'sim_jext_26-Jan-2021_II_D1_3.mat';


load(dat{1})
o1 = o;
load(dat{2})
o2 = o;
load(dat{3})
o3 = o;
dat{6} = 'case_jext_26-Jan-2021_II_C.mat';
load(dat{6})
% compare dat 1 and 2
for i=1:np.n
    er1(i) = norm(o1.u{i}-o2.u{i},inf);
end
norm(er1,inf)
figure

subplot(2,1,1)

bar(er1)
xlim([0 22])
title(['$\|u_i - u_i"\|_{\infty}$'], 'Interpreter','Latex')
% compare dat 1 and 3
for i=1:np.n
    er2(i) = norm(o1.u{i}-o3.u{i},inf);
end
norm(er2,inf)
subplot(2,1,2)
bar(er2)
xlim([0 22])
