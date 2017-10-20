function dy = TGFequations_15s(t,y,p,tcnoise,tcnt)
%==========================================================================
% Evaluate the dy of the dimensional ODE's describing TGF-beta signalling.
% Based on Schmierer et al., 2008
%==========================================================================

% Evaluating the left hand sides of the ODE's written in standard implicit
% for, f(t,y,y')=0


s2total = interp1(tcnt,tcnoise,t); %interpolate the S2total input

%algebraic relation for conserved species
s2c = s2total - (ps2c + s24c + 2*s22c) - (1/ncrat)*(2*s22n + s24n + ps2n + s2n);
%ncrat = ncrat

%differential equations
dy = zeros(15,1);
dRec = -ktgfb*Rec*TgfB;                                %Receptors
dTgfB = -ktgfb*Rec*TgfB;                                %Tgfbeta
dRAct = ktgfb*Rec*TgfB - konsb*RAct*SB + koffsb*RInact;     %Ract
dRInact = konsb*RAct*SB - koffsb*RInact;                      %Rinact
dSB = koffsb*RInact - konsb*RAct*SB;                      %SB-431542

% ds2c = kex2*s2n - kin2*s2c - kphos*s2c*RAct;             %smad2-cyto
ds2c = 0;         %smad2-cyto
dps2c = kex2*ps2n - kin2*ps2c + kphos*s2c*RAct - kon*ps2c*(s4c + 2*ps2c) + koff*(s24c + 2*s22c); %pSmad2-cyto
ds4c = kex4*s4n - kin4*s4c - kon*s4c*ps2c + koff*s24c;    %smad4-cyto
ds24c = kon*ps2c*s4c - koff*s24c - kin2*CIF*s24c;        %smad2-smad4-cyto
ds22c = kon*ps2c*ps2c - koff*s22c - kin2*CIF*s22c;   %smad2-smad2-cyto

ds2n = ncrat*kin2*s2c - ncrat*kex2*s2n + kdephos*ps2n*PPase;        %smad2-nuc
dps2n = ncrat*kin2*ps2c - ncrat*kex2*ps2n - kdephos*ps2n*PPase - kon*ps2n*(s4n + 2*ps2n) + koff*(s24n + 2*s22n); %pSmad2-nuc
ds4n = ncrat*kin4*s4c - ncrat*kex4*s4n - kon*s4n*ps2n + koff*s24n; %Smad4-nuc
ds24n = kon*ps2n*s4n - koff*s24n + ncrat*kin2*CIF*s24c;   %Smad2-Smad4-nuc
ds22n = kon*ps2n*ps2n - koff*s22n + ncrat*kin2*CIF*s22c;
















