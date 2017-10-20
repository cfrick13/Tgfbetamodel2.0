function p = TGFparameters
%

%==========================================================================
% The best-fit parameters in Schmierer et al., 2008
%==========================================================================


kdissSB = 684;      % nM
koffSB = 100;       % /s
konSB = koffSB/kdissSB;  % /nM-s

koff = 0.016;        % /s, koff, page 2 of Supp
kDiss = 8.7;         % nM, kDiss, Fig S4
kon = koff/kDiss;    % /nM-s, kon = koff / Kdiss,

%parameters
p = zeros(16,1);
p(1) = 0.074;        %/nM-s, kTGFbeta
p(2) = konSB;        % /nM-s
p(3) = koffSB;       % /s
p(4) = 0.0056;       % /s, kex, page 1 of Supp 0.0056;  
p(5) = 0.0026;       % /s, kin, page 1 of Supp 0.0026
p(6) = 0.000404;     % /nM-s, kphosp, Fig S4kex
p(7) = kon;          % /nM-s, kon = koff / Kdiss, Kdiss is in Fig S4
p(8) = koff;         % /s
p(9) = 5.7;          % CIF, Fig S4
p(10) =  0.00657;    % /nM-s, kdephosp, Fig S4
p(11) = 0.0026;      % /s, kin, page 1 of Supp 0.0026
p(12) = 0.0026;      % /s, kin, page 1 of Supp 0.0026
p(13) = 1;           % nM, PPase, Table S1 
p(14) = 2.3;         % supplements
p(15) = 60.6+(28.5/2.3); % nM, S2total  
p(16) = 50.8+(50.8/2.3); % nM, S4total  

end