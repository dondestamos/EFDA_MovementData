function EFDA_Whip_Average_ShowVar(data,SubjType,AverMethod)

% Processes and analyzes hand speed data from a motor neuroscience experiment of 
% hitting a target with the bullwhip (Krotov & Russo 2022 in RSOS and 
% Henrot B.S. MIT Thesis 2018). It takes experimental data provided along, 
% a subject type (Novice or Expert), and an averaging method as input. The 
% function aligns time-series data using one of four methods 
% (TimePadLeft, TimePadRight, TimeNorm, EFDA), calculates mean and standard deviation,
% compares peak hand speed estimates, and visualizes spatial and temporal variability.

% Change Nresamp to a lower of a higher value to specifically monitor the effect 
% of sampling rate to time-warping alignment

% Toggle the flag PreserveIntegral for preserving the total distance travelled 
% (the original signals are speed).

% Toggle flag NormByPeakHS for further reducing amplitude differences between the subjects.

% For data descriptuion see the current paper (Krotov, S Razavian, Sadeghi, and Sternad 2024).
% Aleksei Krotov, Northeastern University, 2024



AverMethods = {'TimePadLeft','TimePadRight','TimeNorm','EFDA'};
assert(any(strcmpi(SubjType,{'Novice1','Novice2','Expert'})),'The second argument must be "Novice1", "Novice2", or "Expert"');
assert(any(strcmpi(AverMethod,AverMethods)),'The third argument must be "TimePadLeft", "TimePadRight", "TimeNorm", or "EFDA"');
Cols = [0.7451 0.333 0.5882; 0.64 0.73 0.35; 0.2392 0.5098 0.6941; 0.9961 0.3451 0.0824];
indMethod = find(strcmpi(AverMethod,AverMethods));
ColMain = Cols(indMethod,:);
xlimTS = ismember(indMethod,[3 4]) .* [0 1] + ismember(indMethod,[1 2]) .* [0 2];

Nresamp = 200;
blockrange = 1:5;
flag_PreserveIntegral = 0;
MaxNanRatio = 0.2;
figpos0 = [0 0 0 0]; % Adjust the first argument(s) if plotting not on the main monitor.
flag_NormByPeakHS = 0;


WhipTask_TimeSeries = struct;
WhipTask_TimeSeries(1).Subj = [];
peakHS = nan(201,1);
HandOn = nan(201,1);
iOff = nan(201,1);
ResultStructs = struct;
Times = nan(201,3);

StyleNames = {'Discrete', 'Rhythmic'};
if flag_NormByPeakHS
   disp('HS from each trial are normalized by the peak HS!');
end



% Pack time series and average.
for istyle = 1:2
irow = istyle;
        % Select indices of a participant.
        T = struct2table(data);
        ind = find(strcmpi(T.Participant,SubjType) & T.StyleDR == istyle & ismember(T.Block, blockrange));
        Nind = length(ind);
        TS = cell(Nind,1);
        TSnew = cell(Nind,1);
        IndRange = nan(Nind,1);
        for iitrial = 1:Nind
            itrial = ind(iitrial);
            SampleRate = data(itrial).SampleRate;

            if any(strcmpi(SubjType,{'Novice1','Novice2'}))
                HVel = data(itrial).HandVel;
                HS = sqrt(sum(HVel.^2,2));
            else
                HS = data(itrial).HandSpeed;
            end
            TStemp = HS;

            if flag_NormByPeakHS
                TStemp = TStemp ./ max(TStemp,[],1,'omitnan');
            end

            DUR = data(itrial).Duration;
            IndRange(iitrial,1) = 1; % These could be various time landmarks for selecting an interval of interest
            IndRange(iitrial,2) = DUR; 

            TS{iitrial,1} = TStemp(IndRange(iitrial,1):IndRange(iitrial,2),:);
            TSnew{iitrial,1} = [(1:DUR)' ./ SampleRate TStemp(IndRange(iitrial,1):IndRange(iitrial,2),:)];
        end


    % Using the same timenorm first, because it conveniently removes trials with
        % too many NaNs and already converts cell to double.
        [ynorm,dur,NTrials,NTrialsUsed] = TimeNormSeries(TS,SampleRate,Nresamp,flag_PreserveIntegral,'MaxNanRatio',MaxNanRatio);
        ynormAl = nan(size(ynorm(1),size(ynorm(2))));
    if strcmpi(AverMethod,'TimeNorm')
        [TimeMean,TSmean,TSstd,varmetrics] = AverageTimeSeries_TimeNorm(ynorm,dur);
        % 
    end
    if strcmpi(AverMethod,'TimePadLeft')
        [TimeMean,TSmean,TSstd, ynorm, dur, varmetrics] = AverageTimeSeries_TimePadLeft(ynorm,dur);
    end
    if strcmpi(AverMethod,'TimePadRight')
        [TimeMean,TSmean,TSstd, ynorm, dur, varmetrics] = AverageTimeSeries_TimePadRight(ynorm,dur);
    end

    if strcmpi(AverMethod,'EFDA')
        opts.EFDAWarpsWithDurs = 0;
        opts.EFDAParallel = 1;
        opts.EstimateExecutionTime = 1;
        opts.EFDALambda = 0.01;
        opts.EFDAGraphics = 0;


        [opts, EFDAResult] = EFDA_alignmentMain(TSnew,opts);
        TimeMean = EFDAResult.TimeMean;
        TSmean = EFDAResult.FuncAlignedMean;
        TSstd = EFDAResult.FuncAlignedStd;
        ResultStruct = EFDAResult;
        ResultStructs(irow).ResultStruct = ResultStruct;
        varmetrics.VarAmp = EFDAResult.VarSpat;
        varmetrics.VarTemp = EFDAResult.VarTemp; 
        varmetrics.VarTemp2Spat = EFDAResult.VarTemp2Spat;
        ynormAl = ResultStruct.FuncAligned;

        if opts.EFDAWarpsWithDurs
            fprintf('<strong>EFDA TempVar includes varying durations and is in [s]!</strong> \n');
        end


    end
    varmetrics.VarDur = std(dur,'omitnan');
    % Save it
    WhipTask_TimeSeries(irow).Subj = 1;
    WhipTask_TimeSeries(irow).StyleDR = istyle;
    WhipTask_TimeSeries(irow).NTrials = NTrials;
    WhipTask_TimeSeries(irow).NTrialsUsed = NTrialsUsed;
    WhipTask_TimeSeries(irow).AverMethod = AverMethod;
    WhipTask_TimeSeries(irow).Durations = dur;
    WhipTask_TimeSeries(irow).TimeMean = TimeMean;
    WhipTask_TimeSeries(irow).TSall = squeeze(ynorm);
    WhipTask_TimeSeries(irow).TSAligned = ynormAl;
    WhipTask_TimeSeries(irow).TSmean = TSmean;
    WhipTask_TimeSeries(irow).TSstd = TSstd;
    WhipTask_TimeSeries(irow).VarAmp = varmetrics.VarAmp;
    WhipTask_TimeSeries(irow).VarTemp = varmetrics.VarTemp;
    WhipTask_TimeSeries(irow).VarTemp2Spat = varmetrics.VarTemp2Spat;

    fprintf('Style %d. Variabilities: Amp: %.4f, Temp: %.4f, Temp2Spat: %.4f, durations sd: %.4f \n', istyle, varmetrics.VarAmp, varmetrics.VarTemp, varmetrics.VarTemp2Spat, varmetrics.VarDur);
end
srcTimeSeries = WhipTask_TimeSeries;

%% Add plots of the "Aligned" (or not) data

% if plotting more than, say, 300 lines, pick random ones to pick; otherwise cannot export as vector;
MaxTrials2Plot = 190;
%isubj = 1;
curFigHScurrent = Create_Reuse_Figure([],sprintf('EFDA Manuscript - %s - %s',SubjType,AverMethod),figpos0 + [100 200 1200 400]);
for istyle = 1:2
    irow =  istyle;
    Durations = WhipTask_TimeSeries(irow).Durations;
    ind2plot = 1:length(Durations);
    if length(Durations) > MaxTrials2Plot
        ind2plot = randi(length(Durations),MaxTrials2Plot,1);
    end

    HS = WhipTask_TimeSeries(irow).TSall;
    if strcmpi(AverMethod,'EFDA')
        HS = ResultStructs(irow).ResultStruct.FuncAligned;
    end

    TimeMean = WhipTask_TimeSeries(irow).TimeMean;
    
    spHS(istyle) = subplot(2,9,9*(istyle-1)+(1:4)); hold on;
    title(sprintf('%s',StyleNames{istyle}));
    ylabel('HS (m/s)');
    if flag_NormByPeakHS
        ylabel('HS (a.u. wrt Peak)');
    end
    NFrames = size(HS,1);
    for iitrial = 1:length(ind2plot)
        itrial = ind2plot(iitrial);
        if strcmp(AverMethod,'TimePadLeft')
            xPlot = linspace(0,Durations(itrial),NFrames)';
            yPlot = HS(:,itrial);
            plot(xPlot,yPlot,'Color',[0.5 0.5 0.5],'LineWidth',0.5);
        elseif strcmp(AverMethod,'TimePadRight')
            xPlot = linspace(max(Durations) - Durations(itrial),max(Durations),NFrames)';
            yPlot = HS(:,itrial);
            plot(xPlot,yPlot,'Color',[0.5 0.5 0.5],'LineWidth',0.5);
        else
            xPlot = linspace(0,1,NFrames)';
            yPlot = HS(:,itrial);
            plot(xPlot,yPlot,'Color',[0.5 0.5 0.5],'LineWidth',0.5);
        end
    end
    

end
linkaxes(spHS,'xy');
ylim([0 10]); xlim(xlimTS);
if flag_NormByPeakHS
    ylim([0 1]);
end


%% Plotting Mean+SD

% Plot SD patches
for istyle = 1:2
    irow = istyle;
    X = srcTimeSeries(irow).TimeMean;
    X = [X; flip(X)];

    Ym = srcTimeSeries(irow).TSmean;
    Ys = srcTimeSeries(irow).TSstd;
    Y = [Ym + Ys; flip(Ym - Ys)];
    patch(spHS(istyle),'XData',X,'YData',Y,'FaceColor',ColMain,'FaceAlpha', 0.5,'EdgeColor','none');

end

% Plot mean lines
for istyle = 1:2
    irow = istyle;
    X = srcTimeSeries(irow).TimeMean;
    Y = srcTimeSeries(irow).TSmean;
    plot(spHS(istyle),X,Y,'Color',ColMain,'LineWidth',2);

end




%% Estimates of peak speed
if 1
    for istyle = 1:2
        irow = istyle;
        HS_mean = WhipTask_TimeSeries(irow).TSmean;
        HS_std = WhipTask_TimeSeries(irow).TSstd;
        Time_mean = WhipTask_TimeSeries(irow).TimeMean;
        Durations = WhipTask_TimeSeries(irow).Durations;
        flag_normTime = Time_mean(end) == 1;

        [peak_Mean(istyle),i_peak] = max(HS_mean,[],'omitnan');
        peak_T(istyle) = Time_mean(i_peak) .* ((Time_mean(end)~=1) + (Time_mean(end)==1) * mean(Durations,'omitnan'));
        peak_Std(istyle) = HS_std(i_peak);
        
        ind = find(strcmpi(T.Participant,SubjType) & T.StyleDR == istyle & ismember(T.Block, blockrange));
        peak_T0(istyle) = mean(T.t_HSmax(ind),'omitnan');
        peak_Mean0(istyle) = mean(T.HS_HSmax(ind),'omitnan');
        peak_Std0(istyle) = std(T.HS_HSmax(ind),[],'omitnan');
        if flag_NormByPeakHS
            peak_Mean0(istyle) = 1;
            peak_Std0(istyle) = 0;
        end

    end
    peak_Mean_d = peak_Mean - peak_Mean0;
    peak_Std_d = peak_Std - peak_Std0;
    peak_T_d = peak_T - peak_T0;

    figure(curFigHScurrent);
    spBar = subplot(1,9,5:6);
    Pars2Plot = [peak_Mean_d(1,:); peak_Std_d(1,:); peak_T_d(1,:)];
    ColMap = linspecer(length(Pars2Plot));
    bp = bar(Pars2Plot);
    ParNames = {'Mean','SD','Time'};
    StyleNames = {'Discrete', 'Rhythmic'};
    legend(StyleNames,'location','best');

    bp(1).Horizontal = 'on';
    bp(2).Horizontal = 'on';
    spBar.YTickLabel = ParNames;
    spBar.YDir = 'reverse';
    title('Peak Hand estimate errors:');
    subtitle('Estimated from the mean vs. Mean of estimates');
    %xlim([-2 2]);
    %xlim([-0.4 0.4]);
    
    bp(1).FaceColor = ColMain;
    bp(2).FaceColor = ColMain;
    bp(2).EdgeColor = 'none';
    bp(1).LineWidth = 1;

end






%% Plot variabilities

spBar(1) = subplot(1,9,7); hold on; box off;
ylabel('Var Spatial','fontweight','bold');
xticklabels({});
xticks('auto');
spBar(2) = subplot(1,9,8); hold on; box off;
% ylabel('Var Phas','fontweight','bold');
ylabel('Var Temporal','fontweight','bold');
xticklabels({});
xticks('auto');
spBar(3) = subplot(1,9,9); hold on; box off;
%ylabel('Var TShift','fontweight','bold');
ylabel('Var Temp2Spat','fontweight','bold');

VarDataAmp = nan(1,2);
VarDataTemp = nan(1,2);
VarDataTemp2Spat = nan(1,2);

for iisubj = 1:1
    for istyle = 1:2
        irow = find(ismember([srcTimeSeries.Subj],1) & [srcTimeSeries.StyleDR] == istyle);

        % Variability computed as trapz(std(Aligned))
        VarDataAmp(iisubj,istyle) = srcTimeSeries(irow).VarAmp;
        % Variability computed as trapz(std(TimeShift))
        VarDataTemp(iisubj,istyle) = srcTimeSeries(irow).VarTemp;
        % Variability computed as trapz(std(AlignedMeanWarped))
        VarDataTemp2Spat(iisubj,istyle) = srcTimeSeries(irow).VarTemp2Spat;
    end
end
BarAmp = bar(spBar(1),VarDataAmp);
BarTemp = bar(spBar(2),VarDataTemp);
BarTemp2Spat = bar(spBar(3),VarDataTemp2Spat);
BarAmp.FaceColor = ColMain;
BarTemp.FaceColor = ColMain;
BarTemp2Spat.FaceColor = ColMain;

spBar(1).XTick = 1:2;
spBar(1).XTickLabels = {'D','R'};
linkaxes(spBar,'x');
xlim([0 3]);
ylim(spBar,'auto');

end






function [ynorm,dur,Ntrials,NTrialsUsed,nanTriallist] = TimeNormSeries(src,SampleRate,Nresamp,flag_PreserveIntegral,varargin)

MaxNan = 0.25; % If more than 25% of the data are NaN, the trial is not included in mean and SD calculation.
if any(strcmpi(varargin,'MaxNanRatio'))
    MaxNan = varargin{find(strcmpi(varargin,'MaxNanRatio'))+1};
end

warning('off','MATLAB:interp1:NaNstrip');
ynorm = nan(Nresamp,size(src{1},2),length(src));
dur = nan(length(src),1);
Ntrials = length(src);
flagDiscard = zeros(Ntrials,1);
for iitrial = 1:length(src)
    y0 = src{iitrial};
    if size(y0,1) < size(y0,2)
        y0 = y0'; %Vertical vector/matrix
    end
    L = length(y0);
    dur(iitrial,1) = L ./ SampleRate;

    % Check for large gaps, separately for each dimension
    for idim = 1:size(y0,2)
        nanFramelist = isnan(y0(:,idim));
        Ngaps = nnz(nanFramelist);

        if Ngaps > MaxNan * L | length(y0) == 0
            flagDiscard(iitrial,1) = 1;
            break
        end
        % Check for gaps at the ends
        if nanFramelist(1) || nanFramelist(end)
            flagDiscard(iitrial,1) = 1;
            break
        end

        % interpolate
        irangeorig = linspace(0,1,L)';
        irangenew = linspace(0,1,Nresamp)';
        y1 = interp1My(irangeorig,y0(:,idim),irangenew,'makima','extrap');


        % scale nans locations
        inanorig = find(nanFramelist);
        inannew0 = unique(round(inanorig ./ L * Nresamp));
        % propagate to +-1
        inannew1 = sort(unique([inannew0 - 1 inannew0 inannew0 + 1]));
        inannew1(inannew1>Nresamp) = [];
        inannew1(inannew1<1) = [];

        % add these nans
        y1(inannew1,:) = nan;
        ynorm(:,idim,iitrial)  = y1;
    end
end
warning('on','MATLAB:interp1:NaNstrip');
NTrialsUsed = nnz(~flagDiscard);
nanTriallist = find(flagDiscard);

if flag_PreserveIntegral
    DurNormTemp(1,1,1:length(dur)) = dur ./ mean(dur,'omitnan');
    ynorm = ynorm .*  DurNormTemp;
end
end


function [tmean, ymean, ystd,varmetrics] = AverageTimeSeries_TimeNorm(ynorm,dur)
% squeeze ynorm if needed.
Nresamp = size(ynorm,1);
ymean = mean(ynorm,3,'omitnan');
ystd = std(ynorm,[],3,'omitnan');

% returning [0..1] time
tmean = linspace(0,1,Nresamp)';
tdur = linspace(0,mean(dur),Nresamp)';
% For special use (e.g. whip HS) only: rescaling max time to the mean of the durations
% tmean = tmean .* mean(dur,'omitnan');

% Find variability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check
%varmetrics.VarAmp = trapz(tmean,ystd.^2); % Old
varmetrics.VarAmp = trapz(ystd) ./ length(ymean);
varmetrics.VarAmpWithDur = trapz(tdur,ystd) ./ tdur(end);  % Math should be the same, but due to some numerics, results are ~1% different.
varmetrics.VarTemp = nan;
varmetrics.VarTemp2Spat = nan;
%varmetrics.VarTimeShift = nan;
%varmetrics.VarTimeSpeed = nan;
end


function [tmean, ymean, ystd, ynorm, dur, varmetrics] = AverageTimeSeries_TimePadLeft(ynorm,dur)
% squeeze ynorm if needed.
% = 500; % Original whip, but shouldn't really matter.
Nresamp = size(ynorm,1);

% Remove 3-STD differences, otherwise they deform so bad.
skipflags = abs(dur - mean(dur,'omitnan')) > 3* std(dur,[],'omitnan');
nnz(skipflags);
disp(sprintf('Skipped %d trials that were longer than Mean+3SD of duration',nnz(skipflags)));
Lmax = round(max(dur(~skipflags)) * Nresamp);
ypad = nan(Lmax,size(ynorm,3));

for itrial = 1:length(dur)
    if skipflags(itrial)
        ynorm(:,1,itrial) = nan;
        dur(itrial) = nan;
        continue,
    end
    Ltrial = round(dur(itrial) * Nresamp);
    indold = linspace(0,1,Nresamp);
    indnew = linspace(0,1,Ltrial);
    ypad(1:Ltrial,itrial) = interp1My(indold,ynorm(:,1,itrial),indnew);
end

%error('Finish here!');
ymean = mean(ypad,2,'omitnan');
ystd = std(ypad,[],2,'omitnan');
tmean = (1:Nresamp)' ./ Nresamp .* max(dur);

% downscale (changing sampling rate) for consistent plotting
ymean = interp1My(linspace(0,1,Lmax)',ymean,linspace(0,1,Nresamp)');
ystd = interp1My(linspace(0,1,Lmax)',ystd,linspace(0,1,Nresamp)');

% Find variability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check
%varmetrics.VarAmp = trapz(tmean,ystd.^2); % Old
varmetrics.VarAmp = trapz(tmean,ystd) ./ tmean(end);% ./ mean(dur);
%varmetrics.VarAmp = trapz(ystd) ./ length(ymean);% ./ mean(dur);
varmetrics.VarTemp = nan;
varmetrics.VarTemp2Spat = nan;
%varmetrics.VarTimeShift = nan;
%varmetrics.VarTimeSpeed = nan;
end


function [tmean, ymean, ystd, ynorm, dur, varmetrics] = AverageTimeSeries_TimePadRight(ynorm,dur)
% squeeze ynorm if needed.
% = 500; % Original whip, but shouldn't really matter.
Nresamp = size(ynorm,1);

% Remove 3-STD differences, otherwise they deform so bad.
skipflags = abs(dur - mean(dur,'omitnan')) > 3* std(dur,[],'omitnan');
nnz(skipflags);
disp(sprintf('Skipped %d trials that were longer than Mean+3SD of duration',nnz(skipflags)));
Lmax = round(max(dur(~skipflags)) * Nresamp);
ypad = nan(Lmax,size(ynorm,3));

for itrial = 1:length(dur)
    if skipflags(itrial)
        ynorm(:,1,itrial) = nan;
        dur(itrial) = nan;
        continue,
    end
    Ltrial = round(dur(itrial) * Nresamp);
    indold = linspace(0,1,Nresamp);
    indnew = linspace(0,1,Ltrial);
    ypad((Lmax-Ltrial)+1:end,itrial) = interp1My(indold,ynorm(:,1,itrial),indnew);
end

%error('Finish here!');
ymean = mean(ypad,2,'omitnan');
ystd = std(ypad,[],2,'omitnan');
tmean = (1:Nresamp)' ./ Nresamp .* max(dur);

% downscale (changing sampling rate) for consistent plotting
ymean = interp1My(linspace(0,1,Lmax)',ymean,linspace(0,1,Nresamp)');
ystd = interp1My(linspace(0,1,Lmax)',ystd,linspace(0,1,Nresamp)');

% Find variability
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check
%varmetrics.VarAmp = trapz(tmean,ystd.^2); % Old
varmetrics.VarAmp = trapz(tmean,ystd) ./ tmean(end);% ./ mean(dur);
varmetrics.VarTemp = nan;
varmetrics.VarTemp2Spat = nan;
%varmetrics.VarTimeShift = nan;
%varmetrics.VarTimeSpeed = nan;
end





