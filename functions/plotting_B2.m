hold on, grid on, box on
dat{1} = 'sim_jext_C_10-Feb-2021 _3.mat';
dat{2} = 'sim_jext_C_09-Feb-2021 _1.mat';
dat{3} = 'sim_jext_C_09-Feb-2021 _2.mat';
dat{4} = 'sim_jext_C_10-Feb-2021 _4.mat';

%dat{5} = 'sim_jext_C_09-Feb-2021 _5.mat';

dat{6} = 'case_jext_C_con_11-Feb-2021.mat';
load(dat{6})
figure
hold on, grid on, box on
for c=1:4
    


    load(dat{c});

   Ptr_t = zeros(np.h,1);
    for i=1:np.n
        for j = 1:np.n
            if ~isempty(o.p_tr{i,j})
            %if j > i
                %if c==3
                %    Ptr_t = Ptr_t + max(0,o.p_tr{i,j}(:,end));
               % else
                    Ptr_t = Ptr_t + max(0,o.p_tr{i,j}(:,end));
                %end
            end
        end
    end

%    Pdt = sum(np.Pd);

    plot(Ptr_t,'LineWidth',1.2)
    clearvars('Ptr_t')
end
%plot(Pdt,'--o','color','k','LineWidth',1.5)
title('\textbf{P2P trading with different cost coefficients}','Interpreter','latex')
ylabel('Power [KW]','Interpreter','latex')
legend({[num2str(0.02)],[num2str(0.07)],[num2str(0.4)],[num2str(1)]})%,['cost coef.=',num2str(np.ctrl(4))],['heterogeneous cost coef.']},'Interpreter','latex'  )
%ylim([100, 200])