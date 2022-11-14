%% plot effect of storage units
clear all
clc

% dat{1} = 'sim37b_A_sto1_sce_1.mat';
% dat{2} = 'sim37b_A_sto1_sce_2.mat';
% dat{3} = 'sim37b_A_sto1_sce_3.mat';
% dat{4} = 'case_37b_A_base1.mat';

dat{1} = 'sim37b_B_1112noSTnoDG_1.mat';
%dat{2} = 'sim37b_A_0812_nodg_notrade_sce_1.mat';
%dat{3} = 'sim37b_A_0412_nodg_sce_3.mat';
dat{4} = 'case_37b_B1112noSTnoDG_base.mat';
%dat{1} = 'sim37b_A_1203_notrade_sce_1.mat';
%dat{2} = 'sim37b_A_1203_notrade_sce_2.mat';
%dat{3} = 'sim37b_A_1203_notrade_sce_3.mat';
%dat{4} = 'case_37b_A_base27.mat';

%%
figure




subplot(2,1,1)
hold on, grid on, box on
%load(dat{4})
%n = np.n;
%Pdt = sum(np.Pd(1:50,1:end));
%bar(Pdt)
load(dat{4})
n = np.n;
%np.h=24;
%np.n=50;
Pdt = sum(np.Pd(1:np.n,:));
%plot(Pdt,'-o','color','k','LineWidth',1.5)
%bar(Pdt)
for c=1
    load(dat{c});
    Pdi_t = zeros(np.h,1);
    Pmg_t = zeros(np.h,1);
    for i=1:n
           Pdi_t = Pdi_t + s.p_di{i}(:,end);
           Pmg_t = Pmg_t + s.p_mg{i}(:,end);
    end
    
    plot([Pdi_t+Pmg_t],'LineWidth',1.5)
end

%plot(Pdt,'--o','color','k','LineWidth',1.5)
title('\textbf{Total power from main grid and dispatchable generation}','Interpreter','latex')
ylabel('Power [KW]','Interpreter','latex')
%legend({'total demand','no agent with storage','half of agents with storage','all agents with storage'},'Interpreter','latex'  )
legend({['cost coef.=',num2str(0.01)],['cost coef.=',num2str(np.qtrl(1))],['cost coef.=',num2str(np.qtrl(2))],['cost coef.=',num2str(np.qtrl(3))]},'Interpreter','latex'  )
%ylim([0, 800])

subplot(2,1,2)
hold on, grid on, box on,

for c=1:4
    load(dat{c});

    Ptr_t = zeros(np.h,1);
    
    for i=1:n
        for j = 1:n
            if ~isempty(s.p_tr{i,j}) 
                for k=1:np.h
                    Ptr_t(k,1) = Ptr_t(k,1) + max(0,s.p_tr{i,j}(k,end));
                end
            end
        end
    end

%    Pdt = sum(np.Pd);

    plot(Ptr_t,'LineWidth',1.5)
end    
%end
%plot(Pdt,'--o','color','k','LineWidth',1.5)
title('\textbf{P2P trading}','Interpreter','latex')
ylabel('Power [KW]','Interpreter','latex')
xlabel('time step [hour]','Interpreter','latex')
legend({'no agents with storage','half agents with storage','all agents with storage'},'Interpreter','latex'  )
%ylim([0, 800])
%%
figure
hold on, grid on, box on,
for c=1:3
    load(dat{c});
    Pdi_t = zeros(np.h,1);
    for i=1:n
           Pdi_t = Pdi_t + s.p_di{i}(:,end);
    end
    plot(Pdi_t,'LineWidth',1.5)
end

%%
figure
hold on, grid on, box on,
for c=1:3
    load(dat{c});
    Pmg_t = zeros(np.h,1);
    for i=1:n
           Pmg_t = Pmg_t + max(0,s.p_mg{i}(:,end));
    end
    plot(Pmg_t,'LineWidth',1.5)
end
