function Pl = gen_load(t_lpr)
% Generating a load profile
% W. Ananduta
% 09/10/2020

load('load_profile.mat')




if t_lpr > 0
    Nlp = Load_prof(t_lpr,:);
    mload = max(Nlp);
    % Add some randomness (<=5% of maximum load)
    rl = 0.1*mload*rand([1,length(Nlp)]);
    Pl = Nlp + rl;
else
    Pl = zeros(1,24);
end

end