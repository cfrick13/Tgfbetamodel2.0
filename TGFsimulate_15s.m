function TGFsimulate_15s

close all

mdir = mfilename('fullpath');
[~,b ] = regexp(mdir,'/');
if isempty(b)
    [~,b] = regexp(mdir,'\');
end
parentdir = mdir(1:b(end-1));


%load conditions
c = feval('TGFconditions');

%load parameters
p = feval('TGFparameters');
pp=p;

%Dynamics values
tn    = c(1);     % Time span for integration,seconds
Tgfoff = c(2);
tspan = [0;tn];


%% questions to ask

%1. test the 1:1 relationship between SB inhibited nuclear-Smad level and
%nuclear complex level.

%2. test the parameter dependence of the relationships

%3. test how complex depends on fluctuations in Smad3

%4. determine effect on pSmad3 and smad complex


%% initialize for loop values

cond1num = length(pp);
cond2num = 7;
cond3num = 4;
iterations = cond1num*cond2num*cond3num;
if ~(mod(cond2num,1) == 0)
    cond2num = ceil(iterations./cond1num);
    iterations = cond2num*cond1num;
    disp('adjust iteration number')
end
% iterations = conditionnum*cond2num;
cond1vals = linspace(1,cond1num,cond1num);
cond1mat = ones(cond2num*cond3num,1)*cond1vals;
cond1vec = sort(cond1mat(:));

cond2vals = logspace(log10(0.1),log10(10),cond2num);
cond2mat = (ones(cond1num*cond3num,1)*cond2vals)';
cond2vec = cond2mat(:);

cond3vals = logspace(log10(0.2),log10(1),cond3num+1);
cond3vals(1) = [];
cond3mat = (ones(cond1num*cond2num,1)*cond3vals)';
cond3vec = cond3mat(:);



%%
%initialize smad3 expression profile variables

Tgfon = 1;



%%
%initialize par and spe structures
y0 = TGFconcentrations(pp);
[parameters,species] = relationships(pp,y0');
par = struct(); %parameters struct
spe = struct(); %timecourse struct
rfe = struct(); %response features struct

pfnames = fieldnames(parameters);
for f = 1:length(pfnames)
    pfstr = pfnames{f};
    par(iterations).(pfstr) = parameters.(pfstr);
end

sfnames = fieldnames(species);
for f = 1:length(sfnames)
    sfstr = sfnames{f};
    spe(iterations).(sfstr) = species.(sfstr);
end

species.s2nucObsv = 1:1:100;
T = zeros(size(species.s2nucObsv));
basalLength = 3;
basal = round(length(T)./2);
CC = datastructmaker(T,species,basal,basalLength,'s2nucObsv');
cfnames = fieldnames(CC);
for f = 1:length(cfnames)
    cfstr = cfnames{f};
    if isempty(CC.(cfstr))
        rfe(iterations).(cfstr) = CC.(cfstr);
    else
        rfe(iterations).(cfstr) = CC.(cfstr);
    end
end
    

tim = struct();
tim(iterations).timeVector = 0;

pst = struct();
pst(iterations).perturbationStrength = 0;
%%


 


for i = 1:iterations
    
    pval = cond1vec(i);
    varval = cond2vec(i);
    disp(i)
    p=pp;
    
%     variations = lognrnd(0,0.1,size(p));
    variations = ones(size(p));
%     variations(pval) = varval;
    variations(pval) = varval;
    p = p.*variations;
%     p(14) = 1;
    
    %random noise in Smad3 level
    tlength = tspan(end)+tspan(end)+tspan(end);
    mixtime = 1000; %100
    Mean2 = 1.*p(15); %1
    Sig2 = 0.*p(15); %.4;
    dt = 120; %seconds
    [tcnt,tcnoise] = tcnoisefun(mixtime,Mean2,Sig2,dt,tlength);
    
    %============
    %Time course for basal state
    %============
    y0 = TGFconcentrations(p);
    y0(2) = Tgfoff;
    smad3 = interp1(tcnt,tcnoise,0);
%     smad31 = smad3;
    y0(6) = smad3 - (y0(7) + y0(9) + 2*y0(10)) - (1/p(14))*(2*y0(15) + y0(14) + y0(12) + y0(11));
%     y1 = y0;
%     smad3sum1 = (y0(6) + y0(7) + y0(9) + 2*y0(10)) + (1/p(14))*(2*y0(15) + y0(14) + y0(12) + y0(11))
    
    % Solving the ODEs
    [T1,Y1] = ode15s(@(t,y) TGFequations_15s(t,y,p,tcnoise,tcnt),tspan,y0);
    
    %============
    %Time course for stimulated state
    %============
    y0 = Y1(end,:);
    smad3 = interp1(tcnt,tcnoise,T1(end));
%     smad32 = smad3;
    y0(6) = smad3 - (y0(7) + y0(9) + 2*y0(10)) - (1/p(14))*(2*y0(15) + y0(14) + y0(12) + y0(11));
    y0(2) = cond3vec(i);
%     y2 = y0;
%     smad3sum2 = (y0(6) + y0(7) + y0(9) + 2*y0(10)) + (1/p(14))*(2*y0(15) + y0(14) + y0(12) + y0(11))
    
    % Solving the ODEs
    [T2,Y2] = ode15s(@(t,y) TGFequations_15s(t,y,p,tcnoise,tcnt),tspan+tspan(end),y0);
    
    %============
    %Time course for SB-inhibited state
    %============
    y0 = Y2(end,:);
    smad3 = interp1(tcnt,tcnoise,T2(end));
%     smad33 = smad3;
    y0(6) = smad3 - (y0(7) + y0(9) + 2*y0(10)) - (1/p(14))*(2*y0(15) + y0(14) + y0(12) + y0(11));
%     y3 = y0;
%     smad3sum3 = (y0(6) + y0(7) + y0(9) + 2*y0(10)) + (1/p(14))*(2*y0(15) + y0(14) + y0(12) + y0(11))
%     y0(2) = 1;
    % y(1) = 1;
    % y(3) = 0;
    % y(4) = 0;
    y0(5) = 10*1000;
    
    % Solving the ODEs
    [T3,Y3] = ode15s(@(t,y) TGFequations_15s(t,y,p,tcnoise,tcnt),tspan+tspan(end)+tspan(end),y0);
    
    
    %============
    %Concatentate the time courses for both states
    %============
    T = vertcat(T1,T2(2:end),T3(2:end)) ;
    Ypre = vertcat(Y1,Y2(2:end,:),Y3(2:end,:));
    Y = zeros(length(tcnt),size(Ypre,2));
    for yi = 1:size(Ypre,2)
        YY = interp1(T,Ypre(:,yi),tcnt);
        Y(:,yi) = YY;
    end
    T = (tcnt' - T1(end))./3600;
    Y(:,6) = smad3 - (Y(:,7) + Y(:,9) + 2*Y(:,10)) - (1/p(14))*(2*Y(:,15) + Y(:,14) + Y(:,12) + Y(:,11));
    
%     Y = Ypre;
%     T = (vertcat(T1(1:end-1),T2(1:end-1),T3) - T1(end))./3600 ;
    basal = find(T==0);
    
    [parameters,species] = relationships(p,Y);

    
% f98 = figure(98);
% f98.Position = [1392 281 926 649];
% fnames = fieldnames(species);
% lfn = length(fnames);
% for f=1:length(fnames)
%     fstr = fnames{f};
%     subplot(ceil(sqrt(lfn)),ceil(sqrt(lfn)),f)
%     plot(T*60,species.(fstr))
%     title(fstr)
%     xlim([-200 400])
%     ylim([0 max([1 max(species.(fstr))*2])])
% end

    pfnames = fieldnames(parameters);
    for f = 1:length(pfnames)
        pfstr = pfnames{f};
        par(i).(pfstr) = parameters.(pfstr);
    end
    
    sfnames = fieldnames(species);
    for f = 1:length(sfnames)
        sfstr = sfnames{f};
        spe(i).(sfstr) = species.(sfstr);
    end
    
    
    basalLength = 3;
    CC = datastructmaker(T,species,basal,basalLength,'s2nuccytoObsv');
    cfnames = fieldnames(CC);
    for f = 1:length(cfnames)
        cfstr = cfnames{f};
        if isempty(CC.(cfstr))
            rfe(i).(cfstr) = CC.(cfstr);
        else
            rfe(i).(cfstr) = CC.(cfstr);
        end
    end
    
    tim(i).timeVector = T;
    perturbationStrength = 10.^sum(abs(log10(variations)));
    pst(i).perturbationStrength = perturbationStrength;
    
    

end


%% save data
olddir = pwd;
cd(parentdir)
savedir = [parentdir 'data'];
if isdir(savedir)
else
    mkdir(savedir)
end
cd(savedir)

dateofsim = datestr(now,'yyyy-mm-dd');
savename = ['simdata-' dateofsim '.mat'];
save(savename,'-v6');
cd(pwd)

end






function CC = datastructmaker(T,input,basal,basalLength,protein)
basalidx = (basal-basalLength):basal;
species = input.(protein);
relspecies = species./nanmean(species(basalidx));
rate = gradient(species);% rate
relrate = gradient(species./nanmean(species(basalidx)));

[mrval,mridx] = max(rate(basal+1:end));
CC.maxrate = mrval;%max rate
[mrval,mridx] = max(relrate(basal+1:end));% relative rate
CC.maxrelrate = mrval;
CC.responsetime = T(basal+mridx)*60;

fcidx = find(relspecies > prctile(relspecies,99).*0.99,1,'first');
CC.foldchange = relspecies(fcidx);
CC.foldchangetime = T(fcidx)*60;

peakidx = find(species > prctile(species,99).*0.99,1,'first');
CC.peak = species(peakidx);
CC.peaktime = T(peakidx)*60;

specvec = species(basal+1:end) - species(basal);
t50idx = find(specvec > max(specvec./2),1,'first');
CC.t50 = T(t50idx+basal).*60;
end


function [tcnt,tcnoise] = tcnoisefun(mixtime,Mean2,Sig2,dt,tlength)
tau = (60*60*mixtime)./(-log(0.5)); %50 hours is a mixing time of ~35 hours
mu2 = log((Mean2^2)./sqrt(Sig2.^2+Mean2^2));
sigma2 = sqrt(log(Sig2.^2/(Mean2^2)+1));
Tnoise = 0:dt:tlength;
N = length(Tnoise);
Var2 = normrnd(0,1,[1,N]); %(lambda,[1,N]);
Val = zeros(1,N);
Val(1) = Var2(1);
f = exp(-dt/tau);

for ii = 2:N
    Val(ii) = Val(ii-1)*f+sqrt(1-f^2)*Var2(ii);
end

Val = exp(mu2+sigma2*Val);
%     figure, plot(Tnoise,Val)
%     mean(Val)
%     std(Val)
tcnoise = Val;
tcnt = Tnoise;

end


function cT = extractTraces(exportStruct,indices,speciesstr,finalFrame,stimulationFrame,basalLength,medfiltnum)
cT = struct();

    %define how (median, mean, total) and what to quantify (EGFP, RFP)
    xTracesString = [speciesstr];
    
    % extract the cell traces for the desired number of frames
    cellTracesFull = horzcat(exportStruct(indices).(xTracesString))';
    cellTracesUnfilt = cellTracesFull(:,1:finalFrame); %88x50 [needs to be 50x88]
    % cellTraces = medfilt1(cellTracesUnfilt,medfiltnum,[],2,'omitnan','truncate'); %use median filtering
    % cellTraces = smooth(cellTracesUnfilt,0.1,'sgolay'); %use median filtering
    cellTraces = movmean(cellTracesUnfilt,medfiltnum,2,'Endpoints','shrink');
    cellTraces = cellTracesUnfilt;
    
    % normalize by basal values
    if (stimulationFrame-basalLength)<1
        basalLength = stimulationFrame-1;
    end
    %determine cell traces abundance (normalize to basal population median)
    basalValues = cellTraces(:,stimulationFrame-basalLength:stimulationFrame);
    cellTracesAbundance = cellTraces./nanmedian(basalValues(:));
    
    %determine cell traces difference (subtract basal values of individual cells)
    basalValues = cellTracesAbundance(:,stimulationFrame-basalLength:stimulationFrame);
    basalVector = nanmedian(basalValues,2);
    invBasalVector = 1./(basalVector); %88x1 [and needs to be 88x88]
    invBasalMatrix = ones(size(cellTraces,2),1)*invBasalVector';
    basalMatrix = ones(size(cellTraces,2),1)*basalVector';
    cellTracesDifference = cellTracesAbundance - basalMatrix';
    
    %determine cell traces fold-change
    cellTracesFoldChange =cellTracesAbundance.*(invBasalMatrix');
    
    
    %determine zero to one
    basalValues = cellTraces(:,stimulationFrame-basalLength:stimulationFrame);
    vmin = nanmean(basalValues,2);

    stimulatedValues = cellTraces(:,floor((finalFrame/2))-basalLength:floor((finalFrame/2)));
    vmax = nanmean(stimulatedValues,2);
    
    scaleFactor = 1./(vmax - vmin);
    scaleMatrix = ones(size(cellTraces,2),1)*scaleFactor';
    scaleFactorMin = vmin.*scaleFactor;
    scaleMatrixMin = ones(size(cellTraces,2),1)*scaleFactorMin';
    cellTracesOne = cellTraces.*(scaleMatrix');
    cellTracesZero = cellTracesOne - scaleMatrixMin;
    
    %assign cellTraces to structure for easy access to data.
    cT.(speciesstr).none = cellTraces';
    cT.(speciesstr).abundance = cellTracesAbundance';
    cT.(speciesstr).difference = cellTracesDifference';
    cT.(speciesstr).foldchange = cellTracesFoldChange';
    cT.(speciesstr).zerotoone = cellTracesZero';
end