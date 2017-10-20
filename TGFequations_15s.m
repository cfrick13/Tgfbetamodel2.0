function dy = TGFequations_15s(t,y,p,tcnoise,tcnt)
%==========================================================================
% Evaluate the dy of the dimensional ODE's describing TGF-beta signalling.
% Based on Schmierer et al., 2008
%==========================================================================

% Evaluating the left hand sides of the ODE's written in standard implicit
% for, f(t,y,y')=0


s2total = interp1(tcnt,tcnoise,t); %interpolate the S2total input

%algebraic relation for conserved species
y(6) = s2total - (y(7) + y(9) + 2*y(10)) - (1/p(14))*(2*y(15) + y(14) + y(12) + y(11));
%ncrat = p(14)

%differential equations
dy = zeros(15,1);
dy(1) = -p(1)*y(1)*y(2);                                %Receptors
dy(2) = -p(1)*y(1)*y(2);                                %Tgfbeta
dy(3) = p(1)*y(1)*y(2) - p(2)*y(3)*y(5) + p(3)*y(4);     %Ract
dy(4) = p(2)*y(3)*y(5) - p(3)*y(4);                      %Rinact
dy(5) = p(3)*y(4) - p(2)*y(3)*y(5);                      %y(5)-431542


% dy(6) = p(4)*y(11) - p(5)*y(6) - p(6)*y(6)*y(3);             %smad2-cyto
dy(6) = 0;         %smad2-cyto
dy(7) = p(4)*y(12) - p(5)*y(7) + p(6)*y(6)*y(3) - p(7)*y(7)*(y(8) + 2*y(7)) + p(8)*(y(9) + 2*y(10)); %pSmad2-cyto
dy(8) = p(12)*y(13) - p(11)*y(8) - p(7)*y(8)*y(7) + p(8)*y(9);    %smad4-cyto
dy(9) = p(7)*y(7)*y(8) - p(8)*y(9) - p(5)*p(9)*y(9);        %smad2-smad4-cyto
dy(10) = p(7)*y(7)*y(7) - p(8)*y(10) - p(5)*p(9)*y(10);   %smad2-smad2-cyto

dy(11) = p(14)*p(5)*y(6) - p(14)*p(4)*y(11) + p(10)*y(12)*p(13);        %smad2-nuc
dy(12) = p(14)*p(5)*y(7) - p(14)*p(4)*y(12) - p(10)*y(12)*p(13) - p(7)*y(12)*(y(13) + 2*y(12)) + p(8)*(y(14) + 2*y(15)); %pSmad2-nuc
dy(13) = p(14)*p(11)*y(8) - p(14)*p(12)*y(13) - p(7)*y(13)*y(12) + p(8)*y(14); %Smad4-nuc
dy(14) = p(7)*y(12)*y(13) - p(8)*y(14) + p(14)*p(5)*p(9)*y(9);   %Smad2-Smad4-nuc
dy(15) = p(7)*y(12)*y(12) - p(8)*y(15) + p(14)*p(5)*p(9)*y(10);
















