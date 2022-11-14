function [Pd,Pl,Pr] = gen_Pd_ext(n,t_lpr,r_un,b)
% Generate power disturbance
% W. Ananduta
% 20/07/2019

% Inputs: 
% n = number of agents
% t_lpr = types of load profiles of all agents
% r_un = types of renewable generation units of all agents
% Ts = sampling time (in (0,1], in hour)
% l = generation length (integer, in days)

% Outputs:
% Pd = power disturbance (n x 24*l matrix)
% Pl = load ( n x 24*l matrix)
% Pr = renewable generation ( n x 24*l matrix) 

if nargin < 5
    Ts = 1;
    l=1;
end

if nargin < 6
    l=1;
end

load('pv_data_1d.mat')


for i=1:n+b
    % assign load profile and pv production
    Pl(i,:) = gen_load(t_lpr(i));
    
    if i <=n && r_un(i)>0
        % assign renewable generation profile
        Pr(i,:) = pv_prod_1d(:,r_un(i))'/1000;
        mpv = max(Pr(i,:));
        % add randomness
        rl = 0.05*rand([1,length(Pr(i,:))]);
        Pr(i,:) = Pr(i,:).*(1-rl);
        % multiply production according to load profile
        Pr(i,:) = (1+t_lpr(i)*0.2)*Pr(i,:); 
    else
        Pr(i,:) = zeros(1,24);
    end
end
Pd = Pl-Pr;

end