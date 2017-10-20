function y0 = TGFconcentrations(p)
%
%==========================================================================
% The best-fit parameters in Schmierer et al., 2008
%==========================================================================

% k = kex/kin;
k = p(4)/p(5);
y0 = zeros(15,1);
y0(1) = 1;                  %receptors
y0(2) = 0;                  %tgfbeta
y0(3) = 0;                  %Ractive
y0(4) = 0;                  %Rinactive
y0(5) = 0;                  %SB
y0(6) = p(15)*(k/(1+k));    %smad2-cyto
y0(7) = 0;                  %psmad2-cyto
y0(8) = p(16)/2;        	%smad4-cyto
y0(9) = 0;                  %smad2-smad4-cyto
y0(10) = 0;                 %smad2-smad2-cyto
y0(11) = p(15)/(1+k); %smad2-nuc
% y0(11) = p(14)*p(15)/(1+k); %smad2-nuc
y0(12) = 0;                 %psmad2-nuc
y0(13) = p(16)/2;     %smad4-nuc
% y0(13) = p(14)*p(16)/2;     %smad4-nuc
y0(14) = 0;                 %smad2-smad4-nuc
y0(15) = 0;                 %smad2-smad2-nuc


end




