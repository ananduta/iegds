%%
% dat{1} = 'sim37b_A_1203_sce_1.mat';
% dat{2} = 'sim37b_A_1203_sce_2.mat';
% dat{3} = 'sim37b_A_1203_sce_3.mat';
%dat{1} = 'sim_jext_C_ntr_09-Feb-20211.mat';
%dat{4} = 'sim_jext_C_ntr_09-Feb-20212.mat';
dat{4} = 'sim_jext_C_10-Feb-2021 _4.mat';
dat{1} = 'sim_jext_C_09-Feb-2021 _6.mat';
% dat{5} = 'sim37b_A_1203_notrade_sce_2.mat';
% dat{6} = 'sim37b_A_1203_notrade_sce_3.mat';
% dat{7} = 'case_37b_A_base27.mat';
%%%
figure
%subplot(3,1,1)
hold on, grid on, box on
 load(dat{1})
 J1 = o.J;
 load(dat{4})
 J2 = o.J;
 DJ1 = (J2-J1);%./J2;
%  load(dat{2})
%  J3 = s.J;
%  load(dat{5})
%  J4 = s.J;
%  DJ2 = (J4-J3)./J4;
%  load(dat{3})
%  J5 = s.J;
%  load(dat{6})
%  J6 = s.J;
%  DJ3 = (J6-J5)./J6;
%bar([DJ1 DJ2 DJ3])
bar([DJ1])
%legend('no storage', 'half with storage', 'all with storage')
xlabel('Agents','Interpreter','Latex')
ylabel('Cost difference ','Interpreter','Latex')
%title('Trading - No Trading','Interpreter','Latex')
xlim([0 52])


%%
figure
subplot(3,1,1)
hold on, grid on, box on,
for c=1:3
    load(dat{c});
    Pdi_t = zeros(24,1);
    for i=1:n
           Pdi_t = Pdi_t + s.p_di{i}(:,end);
           
    end
    Pdit(:,c) = Pdi_t;
    plot(Pdi_t,'LineWidth',1.5)
end
title('With trading')
subplot(3,1,2)
hold on, grid on, box on,
for c=1:3
    load(dat{c+3});
    Pdi_t = zeros(24,1);
    for i=1:n
           Pdi_t = Pdi_t + s.p_di{i}(:,end);
           
    end
    Pdit(:,c+3) = Pdi_t;
    plot(Pdi_t,'LineWidth',1.5)
end
title('Without trading')
subplot(3,1,3)
hold on, grid on, box on,
for c=1:3
    plot(Pdit(:,c)-Pdit(:,c+3),'LineWidth',1.5)
end
title('Difference')
legend('no storage', 'half with storage', 'all with storage')
%%
figure
subplot(3,1,1)
hold on, grid on, box on,
for c=1:3
    load(dat{c});
    Pmg_t = zeros(24,1);
    for i=1:n
           Pmg_t = Pmg_t + s.p_mg{i}(:,end);
    end
    Pmgt(:,c) = Pmg_t;
    plot(Pmg_t,'LineWidth',1.5)
end
title('With trading')
subplot(3,1,2)
hold on, grid on, box on,
for c=1:3
    load(dat{c+3});
    Pmg_t = zeros(24,1);
    for i=1:n
           Pmg_t = Pmg_t + s.p_mg{i}(:,end);
    end
    Pmgt(:,c+3) = Pmg_t;
    plot(Pmg_t,'LineWidth',1.5)
end
title('Without trading')
subplot(3,1,3)
hold on, grid on, box on,
for c=1:3
    plot(Pmgt(:,c)-Pmgt(:,c+3),'LineWidth',1.5)
end
title('Difference')
legend('no storage', 'half with storage', 'all with storage')