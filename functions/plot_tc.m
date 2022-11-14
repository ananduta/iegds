%%
figure
hold on; box on; grid on;
title('\textbf{Total traded power}','Interpreter','latex')
m = ['s','*','o','x','^','d','+','d'];
ca{1} = '0';
ca{2} = '1';
ca{3} = '5';
ca{4} ='10';
ca{5} = '20';
ca{6}= '50';
ca{7} = '100';
cc_m = 7;

pt_t = zeros(cc_m,24);
for cc = 1:cc_m
    dat = ['sim5_n',num2str(10),'_tr_c',ca{cc}];
    load(dat);
    
    for i = 1:np.n
        for j= 1:np.n
            if j> i && np.Adj(i,j) ==1
                pt = np.Str{i,j}*s.u{i}(:,end);
                pt_t(cc,:) = pt_t(cc,:) + abs(pt)';
            end
        end
    end


    plot(pt_t(cc,:),'Marker',m(cc),'LineWidth',1);
    tot_pow_trad(cc) = sum(pt_t(cc,:));
    J_t(cc) = s.Jt;
end
legend({'$c^{\mathrm{tr}} = 0$','$c^{\mathrm{tr}} = 1$','$c^{\mathrm{tr}} = 10$','$c^{\mathrm{tr}} = 20$','$c^{\mathrm{tr}} = 100$'},'Interpreter','latex')
xlabel('time [hour]','Interpreter','latex')
ylabel('Power [kW]','Interpreter','latex')
xlim([1 24])
%%
figure
subplot(2,1,1)
x = categorical({'1', '2', '4', '8', '10', '20', '100'});
x = reordercats(x,{'1', '2', '4', '8', '10', '20', '100'});
bar(x,tot_pow_trad(2:8)-tot_pow_trad(1))
subplot(2,1,2)
bar(J_t(2:8)'-J_t(1))
% %%
% figure
% hold on; box on; grid on;
% title('Total power exchanged')
% m = ['s','*','o','x','^'];
% ca{1} = '0';
% ca{2} = '1';
% ca{3} = '4';
% ca{4} ='_rand';
% 
% dat = ['sim5_n',num2str(5),'_tr_c',ca{1}];
% load(dat);
% J_b = s.J;
% for cc = 1:3
%     dat = ['sim5_n',num2str(5),'_tr_c',ca{cc}];
%     load(dat);
%     
% %     dJ(cc-1,:) =s.J-J_b;
% dJ(cc,:) =s.J;
%     
% end
% bar(dJ')
% %legend('no cost','uniform costs','heteregenous costs')
% xlabel('agents')
% ylabel('cost')
% % xlim([1 24])