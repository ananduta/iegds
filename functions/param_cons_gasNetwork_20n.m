function gn = param_cons_gasNetwork_20n(gn)
% set parameters of constraints
% W. Ananduta
% 27/01/2022

% % constant on gas flow equation
% cf = zeros(gn.n);
% 
% for i = 1:gn.n
%     for jj = 1:gn.noN(i)
%         j = gn.N{i}(jj);
%         
%         if j > i
%             cf(i,j) =  1+ 5*rand; % according to 20-node Belgium network, cf ranges in [1,6]
%         end
%     end
% end
% 
% cf = cf + cf';

% cf = cf/gn.scale*sqrt(gn.scale);
% 
% gn.cf = cf.^2; %squared c^f


% % supply limits
% s_min = zeros(gn.n,1);
% s_max = zeros(gn.n,1);
% for i = 1:length(gn.gfg)
%     if gn.gfg(i) == 1
%         %s_min(i) = 500 + 500*rand;
%         
%         s_max(i) = 4000 + 2000*rand;
%     end
% 
% end
% gn.s_min = s_min;
% gn.s_max = s_max;

% supply from main grid
% total gas demand constraint
gn.dg_max = gn.n*gn.Gd_max;
gn.dg_min = 0;

% flow limits
gn.phi_max = 50*ones(gn.n);  % in the datasheet 50

%pc.phi_max = 5*ones(gn.n);
% (squared) pressure limits

gn.psi_min = zeros(gn.n,1);
gn.psi_max = zeros(gn.n,1);


tab = table2array(gn.pres);

for i = 1:length(tab(:,1))
    node = find(gn.pairList==i);
    gn.psi_min(node,1) = tab(i,2)/gn.scale;
    
    gn.psi_max(node,1) = tab(i,3)/gn.scale;
end

gn.psi_min = gn.psi_min.^2;

gn.psi_max = gn.psi_max.^2;


end
