function EFDA_conceptDemos()

figpos0 = [0 0 0 0]; % Change if plotting not on the main monitor
SampleRate = 500;
NResamp = 200;
flag_ShowSRRF = 1; % 0 - show function space and L2 distances on functions, or SRRF space and respective L2 (= Fisher-Rao on functions)
flag_EmulateL2 = 0; % Emulating L2 on original signals.
flag_EmulateDTW = 0; % Using dtw.
DTW_Metric = 'euclidean'; % euclidean, absolute, squared, symmkl

%% Initiate original functions
% This a two-gaussian curve
fH = @(tt,Pars) Pars(1) * exp(-4.*log(2) ./ (Pars(3)).^2 .* (tt-Pars(6)).^2) +...
    Pars(2) * exp(-4.*log(2) ./ (Pars(4)).^2 .* (tt-Pars(6)-Pars(5)).^2);

% Pars contain full description: 
% A1 and A2 amplitudes, B1,B2 are full-width at half-max of each peak, 
% Dt is time between the two peaks. T1 is the time of the first peak
% Trim level is imitating a real movement, setting on/off set from the max.
% The Pars(8) is not used here, it may add noise (add script!)
%      A1  A2  B1    B2   Dt   T0  On/Off Noise
Pars1 = [4; 3; 0.16; 0.16; 0.4; 0.3; 0.00; 0];
Pars2 = [2; 5; 0.16; 0.10; 0.6; 0.2; 0.00; 0];

% Change. Making Pars 1 the same as in the main big example, and inverting
% Pars 2
Pars1 = [3; 4; 0.16; 0.16; 0.4; 0.3; 0.00; 0];
Pars2 = [5; 2; 0.10; 0.16; 0.6; 0.2; 0.00; 0];

% See inside this one to change range boundaries or default means and SDs
[tforig1,~] = fMeasureImit(fH,Pars1',SampleRate);
[tforig2,~] = fMeasureImit(fH,Pars2',SampleRate);

% Color light-blue (former red)
ColLB = [.55 .8 1];
ColLB = [0.8 0.2 0.2];

% Time-norm the two functions
time1 = tforig1(:,1); time2 = tforig2(:,1);
time = linspace(0,1,NResamp)';
f1 = interp1(time1,tforig1(:,2),linspace(time1(1),time1(end),NResamp)');
f2 = interp1(time2,tforig2(:,2),linspace(time2(1),time2(end),NResamp)');
df1 = mydiff(time,f1);
df2 = mydiff(time,f2);
q1 = f_to_srvf(f1,time);
q2 = f_to_srvf(f2,time);
FullDist_Orig = sqrt(trapz(time,(f1-f2).^2) ./ time(end));
AmpDist_Orig = sqrt(trapz(time,(q1-q2).^2) ./ time(end));
fprintf('EFDA Amplitude distance before alignment: %.2f\n',AmpDist_Orig);
ManhDist_Orig = trapz(time,abs(f1-f2)) ./ time(end);
fprintf('Manhattan distance before alignment: %.2f\n',ManhDist_Orig);
DerivativeDist_Orig = sqrt(trapz(time,(df1-df2).^2) ./ time(end));
fprintf('Derivative distance before alignment: %.2f\n',DerivativeDist_Orig);


yname = 'f(t)';
if flag_ShowSRRF, yname = 'q(t)'; end

%% plot the original fn, difference, metrics, SRRFs, difference, metrics
figpos = [100 300 300 300] + figpos0;
fig1 = Create_Reuse_Figure([],'Original',figpos);
os.SpacingsOut = [20 30 50 40];
os.SpacingsBetween = [80 40];
os.YLabels = 'None';
os.XLabels = 'None'; % May change to "First" or to "FirstEvery2nd" or "BotLeftCorn" or "None"
os.XTickMarks = 'All';
os.YTickMarks = 'First';
[sp1, ~] = TiledAxesMy(fig1,2,1,[0.0 0 1 1],os);

Y1 = f1 + flag_ShowSRRF * (q1-f1);
plot(sp1(1),time,Y1,'Color', ColLB); hold on;
Y2 = f2 + flag_ShowSRRF * (q2-f2);
plot(sp1(1),time,Y2,'b'); title(sp1(1),'Original time-series');  ylabel(sp1(1),yname); xlabel(sp1(1),'Time norm (a.u.)');
if flag_ShowSRRF, title(sp1(1),'Original SRRF'); end
plot(sp1(2),time,(Y1-Y2),'k'); ylabel('Difference');xlabel('Time norm (a.u.)'); 
yline(sp1(2),0,'Color','k','LineStyle','--');
textline = sprintf('Dist_{full} = %.2f',FullDist_Orig);
if flag_ShowSRRF, textline = sprintf('Dist_{FR} = %.2f',AmpDist_Orig); end
text(0.02, 1,textline,'Units','norm','Fontsize',10);


%% Pairwise align: blue to red
LambdaEFDA = 0.001;
gam_b2r = optimum_reparam(q1,q2,time,LambdaEFDA,'DP');

% Emulate Euclidean-DTW?
if flag_EmulateL2
    g1 = cumtrapz(time,f1.^2);
    g2 = cumtrapz(time,f2.^2);
    fprintf('Emulating Euclidean DTW!')
    gam_b2r = optimum_reparam(f1,f2,time,LambdaEFDA,'DP');
end

psi_b2r = f_to_srvf(gam_b2r,time);
f2w = warp_f_gamma(f2,gam_b2r,time).';
q2w = f_to_srvf(f2w,time);

if flag_EmulateDTW
    [dist, ix,iy] = dtw(f1,f2,DTW_Metric);
    ix = interp1((1:length(ix))',ix,linspace(ix(1),length(ix),length(time))');
    iy = interp1((1:length(iy))',iy,linspace(iy(1),length(iy),length(time))');
    [~,ir,ic] = unique(ix);
    gam_b2r = interp1(ix(ir),iy(ir),linspace(ix(1),ix(end),length(time))','makima');
    gam_b2r = gam_b2r ./ max(gam_b2r);
    psi_b2r = f_to_srvf(gam_b2r,time);
    f2w = warp_f_gamma(f2,gam_b2r,time).';
    q2w = f_to_srvf(f2w,time);
end

FullDist_Blue2Red = sqrt(trapz(time,(f1-f2w).^2) ./ time(end));
TempDist_Blue2Red = acos(trapz(time,psi_b2r));
AmpDist_Blue2Red = sqrt(trapz(time,(q1-q2w).^2) ./ time(end));
fprintf('EFDA Amplitude distance blue-to-red: %.2f\n',AmpDist_Blue2Red);

ManhDist_Blue2Red = trapz(time,abs(f1-f2w)) ./ time(end);
fprintf('Manhattan distance blue-to-red: %.2f\n',ManhDist_Blue2Red);

df2w = mydiff(time,f2w);
DerivativeDist_Blue2Red = sqrt(trapz(time,(df1-df2w).^2) ./ time(end));
fprintf('Derivative distance blue-to-red: %.2f\n',DerivativeDist_Blue2Red);

TimeShift = gam_b2r - time;
TimeSpeed = mydiff(time,gam_b2r);


figpos = [400 700 300 300] + figpos0;
fig2 = Create_Reuse_Figure([],'Warping to a template',figpos);
[sp2, ~] = TiledAxesMy(fig2,2,1,[0.0 0 1 1],os);

Y1 = f1 + flag_ShowSRRF * (q1-f1);
plot(sp2(1),time,Y1,'Color', ColLB); hold on;
Y2 = f2w + flag_ShowSRRF * (q2w-f2w);
plot(sp2(1),time,Y2,'b'); title(sp2(1),'Warping blue to red');  ylabel(sp2(1),yname); xlabel(sp2(1),'Time norm (a.u.)');
plot(sp2(2),time,(Y1-Y2),'k'); ylabel('Difference');xlabel('Time norm (a.u.)'); 
yline(sp2(2),0,'Color','k','LineStyle','--');
textline = sprintf('Dist_{full} = %.2f',FullDist_Blue2Red);
if flag_ShowSRRF, textline = sprintf('Dist_{FR} = %.2f',AmpDist_Blue2Red); end
text(0.02, 1,textline,'Units','norm','Fontsize',10);
linkaxes([sp1(2) sp2(2)],'y');

% Show warp
figpos = [700 700 180 300] + figpos0;
fig3 = Create_Reuse_Figure([],'Warping to a template - warps',figpos);
[sp3, ~] = TiledAxesMy(fig3,2,1,[0.0 0 1 1],os);
plot(sp3(1),time,time,'Color', ColLB);
plot(sp3(1),time,time,'k--','LineWidth',1);
plot(sp3(1),time,gam_b2r,'b'); hold on;
title(sp3(1),'Warping functions');
text(sp3(1),0.02, 0.9,sprintf('Dist_{temp} = %.2f',TempDist_Blue2Red),'Units','norm','Fontsize',10);
xlabel(sp3(1),'Time norm (a.u.)');
ylabel(sp3(1),'Time warped (a.u.)');


%% Pairwise align: red to blue
gam_r2b = optimum_reparam(q2,q1,time,LambdaEFDA,'DP');

% Emulate Euclidean-DTW?
if flag_EmulateL2
    g1 = cumtrapz(time,f1.^2);
    g2 = cumtrapz(time,f2.^2);
    fprintf('Emulating Euclidean DTW!')
    gam_b2r = optimum_reparam(f2,f1,time,LambdaEFDA,'DP');
end


psi_b2r = f_to_srvf(gam_r2b,time);
f1w = warp_f_gamma(f1,gam_r2b,time).';
q1w = f_to_srvf(f1w,time);

if flag_EmulateDTW
    [dist, ix,iy] = dtw(f1,f2,DTW_Metric);
    ix = interp1((1:length(ix))',ix,linspace(ix(1),length(ix),length(time))');
    iy = interp1((1:length(iy))',iy,linspace(iy(1),length(iy),length(time))');
    [~,ir,ic] = unique(ix);
    gam_r2b = interp1(ix(ir),iy(ir),linspace(ix(1),ix(end),length(time))','makima');
    gam_r2b = gam_r2b ./ max(gam_r2b);
    gam_r2b = 0.95*gam_r2b + 0.05*linspace(0,1,length(time))';
    gam_r2b = invertGamma(gam_r2b)';
    gam_r2b = 1.05 * gam_r2b - 0.05 * linspace(0,1,length(time))';

    psi_b2r = f_to_srvf(gam_r2b,time);
    f1w = warp_f_gamma(f1,gam_r2b,time).';
    q1w = f_to_srvf(f1w,time);
end


FullDist_Red2Blue = sqrt(trapz(time,(f1w-f2).^2) ./ time(end));
TempDist_Red2Blue = acos(trapz(time,psi_b2r));
AmpDist_Red2Blue = sqrt(trapz(time,(q1w-q2).^2) ./ time(end));
fprintf('EFDA Amplitude distance red-to-blue: %.2f\n',AmpDist_Red2Blue);

ManhDist_Red2Blue = trapz(time,abs(f1w-f2)) ./ time(end);
fprintf('Manhattan distance red-to-blue: %.2f\n',ManhDist_Red2Blue);

df1w = mydiff(time,f1w);
DerivativeDist_Red2Blue = sqrt(trapz(time,(df1w-df2).^2) ./ time(end));
fprintf('Derivative distance red-to-blue: %.2f\n',DerivativeDist_Red2Blue);

TimeShift = gam_r2b - time;
TimeSpeed = mydiff(time,gam_r2b);


figpos = [1000 700 300 300] + figpos0;
fig4 = Create_Reuse_Figure([],'Warping to a template 2',figpos);
[sp4, ~] = TiledAxesMy(fig4,2,1,[0.0 0 1 1],os);

Y1 = f1w + flag_ShowSRRF * (q1w-f1w);
plot(sp4(1),time,Y1,'Color', ColLB); hold on;
Y2 = f2 + flag_ShowSRRF * (q2-f2);
plot(sp4(1),time,Y2,'b'); title(sp4(1),'Warping red to blue');  ylabel(sp4(1),yname); xlabel(sp4(1),'Time norm (a.u.)');
plot(sp4(2),time,(Y1-Y2),'k'); ylabel('Difference');xlabel('Time norm (a.u.)'); 
yline(sp4(2),0,'Color','k','LineStyle','--');
textline = sprintf('Dist_{full} = %.2f',FullDist_Red2Blue);
if flag_ShowSRRF, textline = sprintf('Dist_{FR} = %.2f',AmpDist_Red2Blue); end
text(0.02, 1,textline,'Units','norm','Fontsize',10);
linkaxes([sp1(2) sp4(2)],'y');

% Show warp
figpos = [1300 700 180 300] + figpos0;
fig5 = Create_Reuse_Figure([],'Warping to a template - warps 2',figpos);
[sp5, ~] = TiledAxesMy(fig5,2,1,[0.0 0 1 1],os);
plot(sp5(1),time,time,'b');
plot(sp5(1),time,time,'k--','LineWidth',1);
plot(sp5(1),time,gam_r2b,'Color', ColLB); hold on;


title(sp5(1),'Warping functions');
text(sp5(1),0.02, 0.9,sprintf('Dist_{temp} = %.2f',TempDist_Red2Blue),'Units','norm','Fontsize',10);
xlabel(sp5(1),'Time norm (a.u.)');
ylabel(sp5(1),'Time warped (a.u.)');


%% Warp simultaneously
% [psi_both, gam_both, ~] = genwarp(0.5,1,time,'InverseFlag',1,'Shift',pi/2);
% 
% f1w = warp_f_gamma(f1,gam_both,time).';
% q1w = f_to_srvf(f1w,time);
% f2w = warp_f_gamma(f2,gam_both,time).';
% q2w = f_to_srvf(f2w,time);
% 
% FullDist_Both = sqrt(trapz(time,(f1w-f2w).^2) ./ time(end));
% TempDist_Both = acos(trapz(time,psi_both));
% AmpDist_Both = sqrt(trapz(time,(q1w-q2w).^2) ./ time(end));
% fprintf('EFDA Amplitude distance arbitrary warp: %.2f\n',AmpDist_Both);
% 
% ManhDist_Both = trapz(time,abs(f1w-f2w)) ./ time(end);
% fprintf('Manhattan distance arbitrary warp: %.2f\n',ManhDist_Both);
% 
% df1w = mydiff(time,f1w);
% df2w = mydiff(time,f2w);
% DerivativeDist_Both = sqrt(trapz(time,(df1w-df2w).^2) ./ time(end));
% fprintf('Derivative distance arbitrary warp: %.2f\n',DerivativeDist_Both);
% 
% 
% figpos = [500 250 300 300] + figpos0;
% fig6 = Create_Reuse_Figure([],'Same warping',figpos);
% [sp6, ~] = TiledAxesMy(fig6,2,1,[0.0 0 1 1],os);
% 
% Y1 = f1w + flag_ShowSRRF * (q1w-f1w);
% plot(sp6(1),time,Y1,'Color', ColLB); hold on;
% Y2 = f2w + flag_ShowSRRF * (q2w-f2w);
% plot(sp6(1),time,Y2,'b'); title(sp6(1),'Same warping');  ylabel(sp6(1),yname); xlabel(sp6(1),'Time norm (a.u.)');
% plot(sp6(2),time,(Y1-Y2),'k'); ylabel('Difference');xlabel('Time norm (a.u.)'); 
% yline(sp6(2),0,'Color','k','LineStyle','--');
% textline = sprintf('Dist_{full} = %.2f',FullDist_Both);
% if flag_ShowSRRF, textline = sprintf('Dist_{FR} = %.2f',AmpDist_Both); end
% text(0.02, 1,textline,'Units','norm','Fontsize',10);
% linkaxes([sp1(2) sp6(2)],'y');
% 
% % Show warp
% figpos = [800 250 180 300] + figpos0;
% fig7 = Create_Reuse_Figure([],'Same warping - warps',figpos);
% [sp7, ~] = TiledAxesMy(fig7,2,1,[0.0 0 1 1],os);
% plot(sp7(1),time,gam_both,'k'); hold on;
% plot(sp7(1),time,time,'k--','LineWidth',1);
% title(sp7(1),'Warping functions');
% text(sp7(1),0.05, 0.5,sprintf('Dist_{temp} = %.2f',TempDist_Both),'Units','norm','Fontsize',10);
% xlabel(sp7(1),'Time norm (a.u.)');
% ylabel(sp7(1),'Time warped (a.u.)');
% 


%% Pairwise align
y = [f1 f2];

% Emulate Euclidean-DTW?
if flag_EmulateL2
    g1 = cumtrapz(time,f1.^2);
    g2 = cumtrapz(time,f2.^2);
    fprintf('Emulating Euclidean DTW!')
    %gam_b2r = optimum_reparam(f2,f1,time,LambdaEFDA,'DP');
    y = [g1 g2];
end


TimeSeries = cell(size(y,2),1);
for iitrial = 1:size(y,2)
    TimeSeries{iitrial,1} = [time y(:,iitrial)];
end
opts.FSResampleTarget = 'orig';
opts.EFDAWarpsWithDurs = 0;
opts.EFDAParallel = 1;
opts.EstimateExecutionTime = 1;
opts.EFDALambda = 0.01;
opts.EFDAGraphics = 0;
[opts, ResultStruct] = EFDA_alignmentPublishing(TimeSeries,opts);


f1w = ResultStruct.FuncAligned(:,1);
f2w = ResultStruct.FuncAligned(:,2);
if flag_EmulateL2
    f1w = warp_f_gamma(f1,ResultStruct.Warps(:,1),time).';
    f2w = warp_f_gamma(f2,ResultStruct.Warps(:,2),time).';
end

q1w = f_to_srvf(f1w,time);
q2w = f_to_srvf(f2w,time);
gam1 = ResultStruct.Warps(:,1);
gam2 = ResultStruct.Warps(:,2);
Psis(:,1) = f_to_srvf(ResultStruct.Warps(:,1),time);
Psis(:,2) = f_to_srvf(ResultStruct.Warps(:,2),time);
psi1 = Psis(:,1);
psi2 = Psis(:,2);



FullDist_Both = sqrt(trapz(time,(f1w-f2w).^2) ./ time(end));
TempDist_1 = acos(trapz(time,psi1,1));
TempDist_2 = acos(trapz(time,psi2,1));
TempDist_Both = acos(trapz(time,psi1.*psi2,1));
AmpDist_Both = sqrt(trapz(time,(q1w-q2w).^2) ./ time(end));
fprintf('EFDA Amplitude distance simultaneous: %.2f\n',AmpDist_Both);

ManhDist_Both = trapz(time,abs(f1w-f2w)) ./ time(end);
fprintf('Manhattan distance simultaneous: %.2f\n',ManhDist_Both);

df1w = mydiff(time,f1w);
df2w = mydiff(time,f2w);
DerivativeDist_Both = sqrt(trapz(time,(df1w-df2w).^2) ./ time(end));
fprintf('Derivative distance simultaneous: %.2f\n',DerivativeDist_Both);


figpos = [1100 150 300 300] + figpos0;
fig8 = Create_Reuse_Figure([],'Alignment to mean',figpos);
[sp8, ~] = TiledAxesMy(fig8,2,1,[0.0 0 1 1],os);

Y1 = f1w + flag_ShowSRRF * (q1w-f1w);
plot(sp8(1),time,Y1,'Color', ColLB); hold on;
Y2 = f2w + flag_ShowSRRF * (q2w-f2w);
plot(sp8(1),time,Y2,'b'); title(sp8(1),'Alignment to mean');  ylabel(sp8(1),yname); xlabel(sp8(1),'Time norm (a.u.)');
plot(sp8(1),time,.5*(Y1+Y2),'k','LineWidth',1.5); 
plot(sp8(2),time,(Y1-Y2),'k'); ylabel('Difference');xlabel('Time norm (a.u.)'); 
yline(sp8(2),0,'Color','k','LineStyle','--');
textline = sprintf('Dist_{full} = %.2f',FullDist_Both);
if flag_ShowSRRF, textline = sprintf('Dist_{amp} = Dist_{FR} = %.2f',AmpDist_Both); end
text(0.02, 1,textline,'Units','norm','Fontsize',10);
linkaxes([sp1(2) sp8(2)],'y');

% Show warp
figpos = [1400 150 180 300] + figpos0;
fig9 = Create_Reuse_Figure([],'Alignment to mean - warps',figpos);
[sp9, ~] = TiledAxesMy(fig9,2,1,[0.0 0 1 1],os);
plot(sp9(1),time,gam1,'Color', ColLB); hold on;
plot(sp9(1),time,gam2,'b'); hold on;
plot(sp9(1),time,time,'k--','LineWidth',1);
title(sp9(1),'Warping functions');
text(sp9(1),0.02, 0.8,sprintf('Dist_{temp} = %.2f',TempDist_Both),'Units','norm','Fontsize',10);
text(sp9(1),0.02, 0.5,sprintf('Dist_{temp1} = %.2f',TempDist_1),'Units','norm','Fontsize',10);
text(sp9(1),0.02, 0.3,sprintf('Dist_{temp2} = %.2f',TempDist_2),'Units','norm','Fontsize',10);

xlabel(sp9(1),'Time norm (a.u.)');
ylabel(sp9(1),'Time warped (a.u.)');

% Show TimeShift
figpos = [1600 150 180 300] + figpos0;
fig10 = Create_Reuse_Figure([],'Alignment to mean - TimeShifts',figpos);
[sp10, ~] = TiledAxesMy(fig10,2,1,[0.0 0 1 1],os);
plot(sp10(1),time,gam1 - time,'Color', ColLB); hold on;
plot(sp10(1),time,gam2 - time,'b'); hold on;
yline(sp10(1),0,'Color','k','LineStyle','--','LineWidth',1);
%title(sp10(1),'Time Shift');
% text(sp10(1),0.02, 0.9,sprintf('Dist_{temp} = %.2f',TempDist_Both),'Units','norm','Fontsize',10);
% text(sp10(1),0.02, 0.6,sprintf('Dist_{temp1} = %.2f',TempDist_1),'Units','norm','Fontsize',10);
% text(sp10(1),0.02, 0.3,sprintf('Dist_{temp2} = %.2f',TempDist_2),'Units','norm','Fontsize',10);
xlabel(sp10(2),'Time norm (a.u.)');
ylabel(sp10(1),'Time Shift');
plot(sp10(2),time,mydiff(time,gam1),'Color', ColLB); hold on;
plot(sp10(2),time,mydiff(time,gam2),'b'); hold on;
yline(sp10(2),1,'Color','k','LineStyle','--','LineWidth',1);
ylabel(sp10(2),'Time Speed');



% Indicate the mean template and its infinite warps
fmean = ResultStruct.FuncAlignedMean;
qmean = f_to_srvf(ResultStruct.FuncAlignedMean,time);
if 0
    [~, gam_demo(:,1), ~] = genwarp(0.2,.5,time,'InverseFlag',1,'Shift',pi/2);
    [~, gam_demo(:,2), ~] = genwarp(0.2,.5,time,'InverseFlag',2,'Shift',pi/2);
    [~, gam_demo(:,3), ~] = genwarp(0.2,1,time,'InverseFlag',1,'Shift',0);
    [~, gam_demo(:,4), ~] = genwarp(0.2,1,time,'InverseFlag',2,'Shift',0);
    for i = 1:4
        fw_demo(:,i) = warp_f_gamma(fmean,gam_demo(:,i),time);
        qw_demo(:,i) = warp_q_gamma(qmean,gam_demo(:,i),time);
        plot(time,fw_demo(:,i),'Color',[0.6 0.6 0.6],'LineWidth',1.5);
    end
end

% Show Var band and Amp Var
VarY = std([Y1 Y2],1,2);
XPatch = [time; flipud(time)];
YPatch = [fmean+VarY; flipud(fmean-VarY)];
patch(sp8(1),'XData',XPatch,'YData',YPatch,'FaceColor',[.5 .5 .5],'FaceAlpha',.4,'EdgeColor','none');

fmean_w1 = warp_f_gamma(fmean,gam1,time).';
fmean_w2 = warp_f_gamma(fmean,gam2,time).';

% Show these also for TimeShift
TS(:,1) = gam1 - time;
TS(:,2) = gam2 - time;
TSM = mean(TS,2,'omitnan');
TSSD = std(TS,1,2,'omitnan');
YPatch = [TSM+TSSD; flipud(TSM-TSSD)];
patch(sp10(1),'XData',XPatch,'YData',YPatch,'FaceColor',[.5 .5 .5],'FaceAlpha',.4,'EdgeColor','none');

end




%%%%%%
%% AUXILIARY FUNCTIONS
%%%%%%


%%% Prepare for plot (trim at 0.05 or something of max) and numerically find parameters
function [tfout,foutparams] = fMeasureImit(fOrig,SampParamSet,SampleRate,varargin)
% Assuming the original function is on [0,1] interval.
t = (1:SampleRate)' ./ SampleRate; %linspace(0,1,SampleRate)';


% Create noise if necessary - before trimming!

indN = [];
if any(strcmpi(varargin,'Noise'))
    indN = find(any(strcmpi(varargin,'Noise')));
    NoisePars = varargin{indN+1};
    NoiseType = NoisePars{1};
    NoiseFreq = NoisePars{2}; % Frequency for Sinusoidal, Full-width at half-max for Spike
    NoiseAmp = NoisePars{3};
end


%fmin = 0.05;
fmin = SampParamSet(7);

% Skewed gaussian
% A = SampParamSet(1);
% B = SampParamSet(2);
% tC = SampParamSet(3);
% a = SampParamSet(4);
%f = fOrig(t,A,B,tC,a);

% Two gaussians
% A1 = SampParamSet(1);
% A2 = SampParamSet(2);
% B1 = SampParamSet(3);
% B2 = SampParamSet(4);
% tC = SampParamSet(5);
f = fOrig(t,SampParamSet);



[fmax, ifmax] = max(f);
% Do not trim everything < fmin, but use it only at the beginning and at the end.
%fout = f(f >= fmin * fmax);
ion = find(f >= fmin,1,'first');
ioff = find(f >= fmin,1,'last');
fout = f(ion:ioff,1);

tout = (1:length(fout))' ./ SampleRate;

% Spike noise. To provide reasonable corruption, add the spike somewhere in between left
% half-height and right half-height
if ~isempty(indN)
    fNoise = zeros(length(t),1);
    if strcmpi(NoiseType,'Spike')
        iLHM = find(fout >= .5 * SampParamSet(1),1,'first');
        iRHM = find(fout >= .5 * SampParamSet(2),1,'last');
        SpikeLoc = rand .* (tout(iRHM) - tout(iLHM)) + tout(iLHM);
        fNoise = NoiseAmp .* max(SampParamSet(1:2)) .* exp(-4 .* log(2) ./ (NoiseFreq.^2) .* (tout - SpikeLoc).^2);
        fout = fout + fNoise;
    end
    
    if strcmpi(NoiseType,'HighFreqMulti')
        fi = rand .* 2*pi;
        fNoise = NoiseAmp .* fout .* sin(2*pi*NoiseFreq.*tout + fi);
        fout = fout + fNoise;
    end
    
end


tfout = [tout fout];

f = tfout(:,2);
foutparams = TwoGaussParamsEstimate(f,SampleRate,[]);

% add arcL here too>?
tfout(:,3) = cumtrapz(f) .* tout(end) ./ length(f);

end

% Parameter estimate for two-gauss curve:
% duration, peak valueS, peak locationS, FWHM-s, SlopeHML, SlopeHMR
function foutparams = TwoGaussParamsEstimate(f,SampleRate,tfnorm)

if SampleRate == 1
    if ~(min(size(f)) == 2)
        error('To determine sample rate, f must include time vector');
    end
    t = f(:,1);
    f0 = f(:,2);
    SampleRate = length(t) ./ (t(end) - t(1));
    f = f0;
end


t = (1:length(f))' ./ SampleRate;
% Find the two peaks

%[pks,ilocs,w,p] = findpeaks(f,'NPeaks',2,'WidthReference','halfheight');
[pks,ilocs,w,p] = findpeaks(f,'WidthReference','halfheight');
[p, ipk] = maxk(p,2);
pks = pks(ipk);
ilocs = ilocs(ipk);
w = w(ipk);

% What if there's only one peak?
locs = t(ilocs);
tFWHM = w ./ SampleRate;

foutparams.Dur = t(end);
foutparams.A1 = pks(1);
foutparams.TA1 = locs(1);
foutparams.tFWHM1 = tFWHM(1);

foutparams.A2 = nan;
foutparams.TAmin = nan;
foutparams.TA2 = nan;
foutparams.Amin = nan;
foutparams.Dt = nan;
foutparams.tFWHM2 = nan;
if length(pks) > 1
    foutparams.A2 = pks(2);
    foutparams.TA2 = locs(2);
    foutparams.Dt = locs(2) - locs(1);
    foutparams.tFWHM2 = tFWHM(2);
    [~, iAmin] = min(f(ilocs(1):ilocs(2)));
    if ~isempty(iAmin)
        foutparams.Amin = f(iAmin+ilocs(1));
        foutparams.TAmin = t(iAmin+ilocs(1));
    end
end


fout_d = mydiff(t,f);
% I could just find max and min of the derivative, but what if at a ~= 0 those maxima are
% not at half-max anymore.... So straightforward, with a bit of averaging.
[~,iL] = min(abs(f(1:ilocs(1)) - .5*pks(1)));
[~,iR] = min(abs(f(ilocs(end):end) - .5*pks(end)));
iR = iR + ilocs(end);

foutparams.SlopeHML = nan;
foutparams.SlopeHMR = nan;
if iL > 2 && (iL+1) < length(f)
    foutparams.SlopeHML = sum(fout_d([-1:1]+iL)) ./ 3;% Left and right slopes at Half-max level
end
if iR > 2 && (iR+1) < length(f)
    foutparams.SlopeHMR = sum(fout_d([-1:1]+iR)) ./ 3; %Averaging to smoothen the diff noise
end

% Find regular amplitude variability
if ~isempty(tfnorm)
    timemean = tfnorm{1}(:,1);
    y = nan(length(tfnorm{1}),length(tfnorm));
    for itf = 1:length(tfnorm)
        y(:,itf) = tfnorm{itf}(:,2);
    end
    foutparams.VarAmp = trapz(var(y,[],2,'omitnan')) ./ size(y,1) ./ timemean(end);
else
    foutparams.VarAmp = nan;
end


end
%


