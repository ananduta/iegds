function gn = param_cons_6n(gn)
% set parameters of constraints
% W. Ananduta
% 18/06/2021

% constant on gas flow equation
cf = zeros(gn.n);
% cf(1,2) = 50.6;
% cf(2,5) = 37.5;
% cf(2,4) = 50.1;
% cf(3,5) = 43.5;
% cf(5,6) = 45.3;

cf(1,2) = 50;
cf(2,5) = 50;
cf(2,4) = 50;
cf(3,5) = 50;
cf(5,6) = 50;

cf = cf + cf';

cf = cf/gn.scale*sqrt(gn.scale);

gn.cf = cf.^2; %squared c^f


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
gn.dg_max = gn.n*240/gn.scale;
gn.dg_min = 0;

% flow limits
gn.phi_max = 240*gn.n*ones(gn.n)/gn.scale;

%pc.phi_max = 5*ones(gn.n);
% (squared) pressure limits
gn.psi_min = (([85, 100, 105, 110, 120, 130]'/gn.scale).*(1-0.2*rand(6,1))).^2;

gn.psi_max = ((1*[140, 155, 160, 175, 175, 195]'/gn.scale).*(1+0.1*rand(6,1))).^2;

% gn.psi_min = (75/gn.scale)^2*ones(gn.n,1);
% 
% gn.psi_max = (200/gn.scale)^2*ones(gn.n,1);


end
