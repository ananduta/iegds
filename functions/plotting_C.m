%% plot effect of storage units

figure
%subplot(2,1,1)
hold on, grid on, box on

dat{1} = 'sim37b_C_sce_3_cst1.mat';
dat{2} = 'sim37b_C_sce_3_cst5.mat';
dat{3} = 'sim37b_C_sce_3_cst10.mat';
dat{4} = 'sim37b_C_sce_3_cst50.mat';
dat{5} = 'sim37b_C_sce_3_cst100.mat';

for c=1:5


    load(dat{c});

    Pdg_t = zeros(np.h,1);
    Pmg_t = zeros(np.h,1);
    for i=1:np.n
        Ni = sum(np.Adj(i,:));

        %compute a_ni,3
        a = zeros(3+Ni,1);
        a(1,1) = 1;
        %Sdg = kron(eye(np.h),a');
        %Pdg_t = Pdg_t + Sdg*s.u{i}(:,end);
        Pmg_t = Pmg_t + s.p_mg{i}(:,end);
        Pdg_t = Pdg_t + s.p_di{i}(:,end);
    end

    

    plot([Pmg_t+Pdg_t],'LineWidth',1.2)
    
end
Pdt = sum(np.Pd(1:25,1:end));
plot(Pdt,'--o','color','k','LineWidth',1.5)
title('\textbf{Scenario a: no P2P trading}','Interpreter','latex')
ylabel('Power [KW]','Interpreter','latex')
%ylim([0, 800])


% subplot(2,1,2)
% %%
% figure
% hold on, grid on, box on,
% 
% dat{1} = 'sim37b_C_sce_4_cst1.mat';
% dat{2} = 'sim37b_C_sce_4_cst5.mat';
% dat{3} = 'sim37b_C_sce_4_cst10.mat';
% dat{4} = 'sim37b_C_sce_4_cst50.mat';
% dat{5} = 'sim37b_C_sce_4_cst100.mat';
% 
% for c=1:5
%     load(dat{c});
%     Pdg_t = zeros(np.h,1);
%     for i=1:np.n
%         Ni = sum(np.Adj(i,:));
% 
%         %compute a_ni,3
%         a = zeros(3+Ni,1);
%         a(1,1) = 1;
%         Sdg = kron(eye(np.h),a');
%         Pdg_t = Pdg_t + Sdg*s.u{i}(:,end);
%     end
% 
% 
% Pdt = sum(np.Pd);
% plot([s.sigma_mg(:,end)+Pdg_t],'LineWidth',1.2)
% end
% plot(Pdt,'--o','color','k','LineWidth',1.5)
%title('\textbf{Scenario b: with P2P trading}','Interpreter','latex')
% ylabel('Power [KW]','Interpreter','latex')
xlabel('time step [hour]','Interpreter','latex')
legend({'cost coef.=1','cost coef.=5','cost coef.=10','cost coef.=50','cost coef.=100','total demand'},'Interpreter','latex'  )
%ylim([0, 800])

