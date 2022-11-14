function [s,p] = TS_ED(p)
% Two-stage method for economic dispatch game of IEGDS with PWA gas-flow
% model

% W. Ananduta
% 04/10/2021



%% Stage 1
% Solve convexified GNEP via proximal-point method
[s,p] = solve_GNEPc(p);

%% Stage 2

% Obtain binary decision variables
[s,p] = solve_binary(s,p);


% Minimize violation of PWA gas-flow model 
[s,p] = solve_pressure(s,p); %centralized

end