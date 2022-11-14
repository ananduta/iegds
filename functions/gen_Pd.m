function [Pd,Pl,Pr] = gen_Pd(n,t_lpr,r_un,Ts,l)
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

if nargin < 4
    Ts = 1;
    l=1;
end

if nargin < 5
    l=1;
end

np.Pl = zeros(n,l);
for i=1:n
    % assign load profile
    Pl(i,:) = gen_load(t_lpr(i),Ts,l);

    % assign renewable generation profile
    Pr(i,:) = gen_rprod(r_un(i),Ts,l);
    % multiply production according to load profile
    Pr(i,:) = t_lpr(i)*Pr(i,:); 
end
Pd = Pl-Pr;

end