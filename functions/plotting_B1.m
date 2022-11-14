%% plot P2P trading analysis
clear all
%close all

figure
%subplot(2,1,1)
hold on, grid on, box on

dat{1} = 'sim_jext_02-Feb-2021_D_0.mat';
dat{2} = 'sim_jext_02-Feb-2021_D_1.mat';
dat{3} = 'sim_jext_02-Feb-2021_D_2.mat';
dat{4} = 'sim_jext_02-Feb-2021_D_3.mat';
dat{5} = 'sim_jext_02-Feb-2021_D_4.mat';
%dat{3} = 'sim37b_B_0412_3.mat';
%dat{4} = 'sim37b_B_0412_4.mat';
%dat{5} = 'sim37b_B_0412_5.mat';
%dat{6} = 'sim37b_B_sce_3_ct500.mat';
dat{6} = 'case_jext_02-Feb-2021_D.mat';
load(dat{6})
for c=1:5
    


    load(dat{c});

   Ptr_t = zeros(np.h,1);
    for i=1:np.n
        for jj = 1:length(np.N{i})
            j = np.N{i}(jj);
            %if j > i
                %if c==3
                %    Ptr_t = Ptr_t + max(0,o.p_tr{i,j}(:,end));
               % else
                    Ptr_t = Ptr_t + max(0,o.p_tr{i,j}(:,end));
                %end
            %end
        end
    end

%    Pdt = sum(np.Pd);

    plot(Ptr_t,'LineWidth',1.2)
    Ptr_total(c)=sum(Ptr_t);
    clearvars('Ptr_t')
end
%plot(Pdt,'--o','color','k','LineWidth',1.5)
title('\textbf{P2P trading with different prices of trading}','Interpreter','latex')
ylabel('Power [KW]','Interpreter','latex')
%legend({['cost coef.=',num2str(o.ct(1))],['cost coef.=',num2str(o.ct(2))],['cost coef.=',num2str(o.ct(1))],['cost coef.=',num2str(o.ct(4))]},'Interpreter','latex'  )
%legend({['cost coef.=',num2str(o.ct(1))],['cost coef.=',num2str(o.ct(2))]},'Interpreter','latex'  )
%ylim([100, 200])

 