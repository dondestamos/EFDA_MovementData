function EFDAFig = EFDA_PlotInput(ynorm,tnorm,dur,opts,EFDAFig,figname,figpos)
% pass also the original&dur&normed plot data if gui enabled
if isempty(EFDAFig) %%% OR IS CLOSED??? i.e. isdeletede?
    EFDAFig = Create_Reuse_Figure([],figname,figpos);
end
EFDAFig.Color = [1 1 1];

durMean = 1;
TimeNormLabel = 'Time norm-d (a.u.)';
if opts.EFDAWarpsWithDurs
    durMean = mean(dur,'omitnan');
    TimeNormLabel = 'Time mean (s)';
end

% check if already exist
posLeftCol = [0.07 0.75 0.21 0.2];
axOrig = subplot('Position',posLeftCol,'Units','normalized'); hold on; title('Original data');
xticklabels({}); xticks('auto');
axDurBP = subplot('Position',posLeftCol - [0 0.04 0 0.17],'Units','normalized'); hold on;
xlabel('Time (s)','fontweight','bold','FontSize',10);
axOrigNorm = subplot('Position',posLeftCol - [0 0.35 0 0],'Units','normalized'); hold on;
xlabel(TimeNormLabel,'fontweight','bold','FontSize',10);
% end check

% variabilities
TxVarTitle = uicontrol ('Style','text','Units','normalized','Position',         [0.4 0.16 0.2 0.03],'String','Variabilities','FontSize',10,'FontWeight','bold','ForegroundColor','k','BackgroundColor',[1 1 1],'HorizontalAlignment','center');

TxVarOrigName = uicontrol ('Style','text','Units','normalized','Position',      [0.45 0.12 0.07 0.03],'String','Orig','FontSize',10,'FontWeight','normal','ForegroundColor',[1 .6 0],'BackgroundColor',[1 1 1],'HorizontalAlignment','left');
TxVarOrig = uicontrol ('Style','text','Units','normalized','Position',          [0.45 0.09 0.07 0.03],'String','','FontSize',10,'FontWeight','normal','ForegroundColor',[1 .6 0],'BackgroundColor',[1 1 1],'HorizontalAlignment','left');

% Alignment and display parameters
% warning('Display parameters!');


if opts.NLongOutlier > 0
    TxLongOutliers = uicontrol ('Style','text','Units','normalized','Position',          [0.05 0.29 0.25 0.02],'String','','FontSize',8,'FontWeight','normal','ForegroundColor','k','BackgroundColor',[1 1 1],'HorizontalAlignment','left');
    TxLongOutliers.String = sprintf('%d trials were removed as 3SD outliers by duration',opts.NLongOutlier);
end

TxNLinesPlotMax = uicontrol ('Style','text','Units','normalized','Position',          [0.05 0.26 0.25 0.02],'String','','FontSize',8,'FontWeight','normal','ForegroundColor','k','BackgroundColor',[1 1 1],'HorizontalAlignment','left');
if ~isempty(opts.NLinesPlotMax) && ~opts.NLinesPlotMax == 0
    TxNLinesPlotMax.String = sprintf('Only %d lines are plotted out of %d',opts.NLinesPlotMax,size(ynorm,2) - opts.NLongOutlier);
elseif size(ynorm,2) > 100
    TxNLinesPlotMax.String = sprintf('All %d lines are plotted. Vector-export might not work',size(ynorm,2) - opts.NLongOutlier);
end

% Do not plot more than NLinesPlotMax, pick a uniform subsample instead
NSampPlot = opts.NSampPlot;
if isempty(NSampPlot) || isnan(NSampPlot) || NSampPlot == 0
    NSampPlot = length(tnorm);
end
NLinesPlotMax = opts.NLinesPlotMax;
if isempty(NLinesPlotMax) || isnan(NLinesPlotMax) || NLinesPlotMax == 0
    NSampPlot = length(dur);
end
M = length(dur) - length(opts.nanTriallist);
trialplotlist = 1:M;
if NLinesPlotMax > 0 && NLinesPlotMax < length(trialplotlist)
    trialplotlist = sort(unique(randperm(M,NLinesPlotMax)));
end
trialindnotnan = 1:length(dur);
trialindnotnan(opts.nanTriallist) = [];
trialplotlist = trialindnotnan(trialplotlist);

% Plot
for itrial = trialplotlist
    tplot = linspace(0,dur(itrial),NSampPlot)';
    yplot = interp1My(tnorm .* dur(itrial),ynorm(:,itrial),tplot);
    line(axOrig,tplot,yplot,'Color',[.6 .6 .6],'LineWidth',1);
    tplot = linspace(0,durMean,NSampPlot)';
    line(axOrigNorm,tplot,yplot,'Color',[.6 .6 .6],'LineWidth',1);
end

% Plot original spatiotemporal variability
ysd = interp1My(tnorm .* durMean,std(ynorm,[],2,'omitnan'),tplot);
yplot = interp1My(tnorm .* durMean,mean(ynorm,2,'omitnan'),tplot);
patch(axOrigNorm,'XData',[tplot; flipud(tplot)],'YData',[yplot+ysd; flipud(yplot-ysd)],'FaceColor',[1 .6 0],'FaceAlpha',0.3,'EdgeColor','none');
plot(axOrigNorm,tplot,yplot,'Color',[1 .6 0]);

% scatter of durations and a boxplot on top
scatter(axDurBP,dur(trialplotlist),.98+.04*rand(length(trialplotlist),1),6,[0.6 0.6 0.6],'filled','MarkerEdgeColor','none');
durbp = boxplot(axDurBP,dur,'BoxStyle','outline','MedianStyle','line','Widths',0.07,'Orientation','horizontal');
objs = axDurBP.Children.Children;
delete(findobj(objs,'Tag','Outliers'));
otemp = findobj(objs,'Tag','Median');
otemp.Color = [0 0 0]; otemp.LineWidth = 1;
otemp = findobj(objs,'Tag','Box');
otemp.Color = [0 0 0];
ylim(axDurBP,[0.95 1.05]);
yticklabels(axDurBP,{}); yticks(axDurBP,[]);

linkaxes([axOrig axDurBP],'x');
linkaxes([axOrig axOrigNorm],'y');
VarOrigSpat = trapz(std(ynorm,[],2,'omitnan')) ./ length(tnorm) ./ tnorm(end); % [original signal units]

TxVarOrig.String = sprintf('%.4f',VarOrigSpat);

% Save handles and bits of plotting settings to the figure handle
EFDAFig.UserData.Data = {TimeNormLabel,durMean,axOrig,axOrigNorm,axDurBP,NSampPlot,trialplotlist};



end