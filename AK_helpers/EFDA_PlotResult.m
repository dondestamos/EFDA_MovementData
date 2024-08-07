function EFDAFig = EFDA_PlotResult(EFDAResult,opts,EFDAFig)

% Load handles and settings
TimeNormLabel = EFDAFig.UserData.Data{1};
durMean = EFDAFig.UserData.Data{2};
axOrig = EFDAFig.UserData.Data{3};
axOrigNorm = EFDAFig.UserData.Data{4};
axDurBP = EFDAFig.UserData.Data{5};
NSampPlot = EFDAFig.UserData.Data{6};
trialplotlist = EFDAFig.UserData.Data{7};

% aligned and warps
posCentCol = [0.4 0.75 0.21 0.2];
axAligned = subplot('Position',posCentCol,'Units','normalized'); hold on; title('Aligned data');
axWarps = subplot('Position',posCentCol - [0 0.45 0 -0.15],'Units','normalized'); hold on; title('Warps');
ylabel(axWarps,'Time of original samples norm-d (a.u.)','fontweight','bold','FontSize',10);
xlabel(TimeNormLabel,'fontweight','bold','FontSize',10);

% meanwarped, timespeed and timeshift
posRightCol = [0.7 0.75 0.21 0.2];
axTimeSpeed = subplot('Position',posRightCol - [0 0.45 0.12 0.05],'Units','normalized'); hold on;
xlabel(TimeNormLabel,'fontweight','bold','FontSize',10);
ylabel({'Time-speed of','original samples (a.u.)'},'fontweight','bold','fontsize',10);
axTimeShift = subplot('Position',posRightCol - [0 0.25 0.12 0.05],'Units','normalized'); hold on;
xticklabels({}); xticks('auto');
ylabel({'Time-shift of','original samples (a.u.)'},'fontweight','bold','fontsize',10);
ylim(axWarps,[0 1]);
if opts.EFDAWarpsWithDurs
    ylabel(axWarps,'Time of original samples (s)','fontweight','bold');
    ylabel(axTimeShift,{'Time-shift of','original samples (s)'},'fontweight','bold','fontsize',10);
    ylabel(axTimeSpeed,{'Time-speed of','original samples (s)'},'fontweight','bold','fontsize',10);
    ylim(axWarps,'auto');
end
axMeanWarped = subplot('Position',posRightCol); hold on; title('Only temporal var');
xlabel(TimeNormLabel,'fontweight','bold','FontSize',10);

% variabilities
TxVarSpatName = uicontrol ('Style','text','Units','normalized','Position',      [0.45 0.05 0.07 0.03],'String','Spat','FontSize',10,'FontWeight','normal','ForegroundColor','r','BackgroundColor',[1 1 1],'HorizontalAlignment','left');
TxVarSpat = uicontrol ('Style','text','Units','normalized','Position',          [0.45 0.02 0.07 0.03],'String','','FontSize',10,'FontWeight','normal','ForegroundColor','r','BackgroundColor',[1 1 1],'HorizontalAlignment','left');

TxVarTempName = uicontrol ('Style','text','Units','normalized','Position',      [0.55 0.12 0.07 0.03],'String','Temp','FontSize',10,'FontWeight','normal','ForegroundColor',[0.78 0.25 0.76],'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
TxVarTemp = uicontrol ('Style','text','Units','normalized','Position',          [0.55 0.09 0.07 0.03],'String','','FontSize',10,'FontWeight','normal','ForegroundColor',[0.78 0.25 0.76],'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
TxVarTemp2SpatName = uicontrol ('Style','text','Units','normalized','Position', [0.55 0.05 0.07 0.03],'String','Temp2Spat','FontSize',10,'FontWeight','normal','ForegroundColor',[.2 .82 .78],'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
TxVarTemp2Spat = uicontrol ('Style','text','Units','normalized','Position',     [0.55 0.02 0.07 0.03],'String','','FontSize',10,'FontWeight','normal','ForegroundColor',[.2 .82 .78],'BackgroundColor',[1 1 1],'HorizontalAlignment','left');

% After EFDA
% Plot the aligned, warps, timeshift, and timespeed after complete efda
for itrial = trialplotlist
    tnorm = linspace(0,durMean,length(EFDAResult.FuncAlignedMean))';
    tplot = linspace(0,durMean,NSampPlot)';
    yplot = interp1My(tnorm,EFDAResult.FuncAligned(:,itrial),tplot);
    gammplot = interp1My(tnorm,EFDAResult.Warps(:,itrial),tplot);
    tshplot = interp1My(tnorm,EFDAResult.WarpsMinusUnity(:,itrial),tplot);
    tspplot = interp1My(tnorm,EFDAResult.WarpsDdt(:,itrial),tplot);
    ymwplot = interp1My(tnorm,EFDAResult.FuncMeanWarped(:,itrial),tplot);

    line(axAligned,tplot,yplot,'Color',[.6 .6 .6],'LineWidth',1);
    line(axWarps,tplot,gammplot,'Color',[.6 .6 .6],'LineWidth',1);
    line(axTimeShift,tplot,tshplot,'Color',[.6 .6 .6],'LineWidth',1);
    line(axTimeSpeed,tplot,tspplot,'Color',[.6 .6 .6],'LineWidth',1);
    line(axMeanWarped,tplot,ymwplot,'Color',[.6 .6 .6],'LineWidth',1);
end
yplot = interp1My(tnorm,EFDAResult.FuncAlignedMean,tplot);
ysd = interp1My(tnorm,EFDAResult.FuncAlignedStd,tplot);
gammsd = interp1My(tnorm,std(EFDAResult.Warps,[],2,'omitnan'),tplot);
gammplot = tplot;
timeshiftsd = interp1My(tnorm,std(EFDAResult.WarpsMinusUnity,[],2,'omitnan'),tplot);
timeshiftmean = interp1My(tnorm,mean(EFDAResult.WarpsMinusUnity,2,'omitnan'),tplot);
ymwplot = interp1My(tnorm,EFDAResult.FuncMeanWarpedMean,tplot);
ymwsd = interp1My(tnorm,EFDAResult.FuncMeanWarpedStd,tplot);

patch(axAligned,'XData',[tplot; flipud(tplot)],'YData',[yplot+ysd; flipud(yplot-ysd)],'FaceColor','r','FaceAlpha',0.3,'EdgeColor','none');
line(axAligned,tplot,yplot,'Color','r','LineWidth',2);
patch(axWarps,'XData',[tplot; flipud(tplot)],'YData',[gammplot+gammsd; flipud(gammplot-gammsd)],'FaceColor',[0.78 0.25 0.76],'FaceAlpha',0.3,'EdgeColor','none');
line(axWarps,tplot,gammplot,'Color','k','LineWidth',2,'LineStyle','--');
patch(axTimeShift,'XData',[tplot; flipud(tplot)],'YData',[timeshiftsd; flipud(-timeshiftsd)],'FaceColor',[0.78 0.25 0.76],'FaceAlpha',0.4,'EdgeColor','none');
line(axTimeShift,tplot,timeshiftmean,'Color','k','LineStyle','--','LineWidth',2);
yline(axTimeSpeed,1,'Color','k','LineStyle','--','LineWidth',2);
patch(axMeanWarped,'XData',[tplot; flipud(tplot)],'YData',[ymwplot+ymwsd; flipud(ymwplot-ymwsd)],'FaceColor',[.2 .82 .78],'FaceAlpha',0.4,'EdgeColor','none');
line(axMeanWarped,tplot,ymwplot,'Color',[.1 .4 .45],'LineWidth',1);

linkaxes([axOrig axAligned axOrigNorm axMeanWarped],'y');
linkaxes([axOrigNorm axAligned axWarps axMeanWarped axTimeShift axTimeSpeed],'x');

TxVarSpat.String = sprintf('%.4f',EFDAResult.VarSpat);
TxVarTemp.String = sprintf('%.4f',EFDAResult.VarTemp);
TxVarTemp2Spat.String = sprintf('%.4f',EFDAResult.VarTemp2Spat);
end