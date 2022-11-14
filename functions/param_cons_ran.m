function gn = param_cons_ran(gn)
% set parameters of constraints
% W. Ananduta
% 18/06/2021

% constant on gas flow equation
cf = zeros(gn.n);
for i = 1:gn.n
    for j =1:gn.n
        if j > i && ~isempty(find(gn.N{i}==j))
            cf(i,j) = 35 ;%+ 15*rand;
        end
        
    end
end


cf = cf + cf';


gn.cf = cf.^2; %squared c^f


% supply limits
s_min = zeros(gn.n,1);
s_max = zeros(gn.n,1);
for i = 1:length(gn.gfg)
    if gn.gfg(i) == 1
        %s_min(i) = 500 + 500*rand;
        
        s_max(i) = 4000 + 2000*rand;
    end

end
gn.s_min = s_min;
gn.s_max = s_max;

% supply from main grid
% total gas demand constraint
gn.dg_max = gn.n*100;
gn.dg_min = 0;

% flow limits
gn.phi_max = 30*gn.n*ones(gn.n);

%pc.phi_max = 5*ones(gn.n);
% (squared) pressure limits
gn.psi_min = 30*ones(gn.n,1);%(85*(1+0.5*rand(gn.n,1))).^2;
gn.psi_max = 170*ones(gn.n,1);%(170*(1+0.5*rand(gn.n,1))).^2;


end
