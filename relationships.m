function [parameters,species] = relationships(p,y)

parameters.ktgfb = p(1);
parameters.konsb = p(2);
parameters.koffsb = p(3);
parameters.kex2 = p(4);
parameters.kin2 = p(5);
parameters.kphos = p(6);
parameters.kon = p(7);
parameters.koff = p(8);
parameters.CIF = p(9);
parameters.kdephos = p(10);
parameters.kin4 = p(11);
parameters.kex4 = p(12);
parameters.PPase = p(13);
parameters.ncrat = p(14);
parameters.S2total = p(15);
parameters.S4total = p(16);

%basic values
species.R = y(:,1);
species.Tgfb = y(:,2);
species.Ract = y(:,3);
species.Rinact = y(:,4);
species.SB = y(:,5);
species.s2c = y(:,6);
species.ps2c = y(:,7);
species.s4c = y(:,8);
species.s24c = y(:,9);
species.s22c = y(:,10);
species.s2n = y(:,11);
species.ps2n = y(:,12);
species.s4n = y(:,13);
species.s24n = y(:,14);
species.s22n = y(:,15); 

%specific values
species.s2nucObsv = species.s2n + species.s24n + (species.s22n).*2 + species.ps2n;
species.s2cytoObsv = species.s2c + species.s24c + (species.s22c).*2 + species.ps2c;
species.s2nuccytoObsv = ((species.s2nucObsv))./(species.s2cytoObsv);
% species.s2total = (species.s2nucObsv) + (1./parameters.ncrat).*(species.s2cytoObsv);
species.s2total = (1./parameters.ncrat).*(species.s2nucObsv) + (species.s2cytoObsv);
% species.s2total = (species.s2nucObsv) + (species.s2cytoObsv);


end