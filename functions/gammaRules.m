function [gammaNew,gammaLower,gammaUpper]=gammaRules(er,gammaNow,gammaLower,gammaUpper,r)
% Line search rules to update penalty weight gamma
% W. Ananduta
% 12/11/2021

c = 1.5;
d = 0.6;
eps = 5e-5;

    if er > eps
        gammaLower = gammaNow;
        gammaNew = min( gammaNow+ 0.05*c^(r+1), gammaUpper*d+(1-d)*gammaLower);
    else
        gammaUpper = gammaNow;
        gammaNew = d*gammaUpper + (1-d)*gammaLower;
    end


end