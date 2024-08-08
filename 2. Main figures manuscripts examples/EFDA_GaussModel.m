function EFDA_GaussModel()
% A comprehensive tool for generating and analyzing synthetic data to evaluate 
% the performance of different signal processing methods, specifically alignment 
% and resampling techniques to extract meaningful information from noisy data. 
% A triple-Gaussian signal is used as an example, with the peak parameters specified 
% in the second sub-function, and optional noise values (also a Gaussian noise, 
% added to one or few of the parameters) are specified in the third sub-function. 
% This function was used for the submitted manuscript by Krotov, S. Razavian, Sadeghi, and Sternad (2024).

% 1) N (default 500) noisy realizations of the original triple-Gaussian signal 
% are generated; Each signal is trimmed with on/offset at 2% of the signal’s maximum 
% imitating an experimental measurement. Peaks are found and characterized (peak value, 
% peak location, peak width at half-maximum). 

% 2) Three approaches are then used on the whole ensemble to extract time-series 
% of mean and standard deviation. Time-padding approach aligns all the signals at 
% their start (trimmed onset) and pads all but the longest signal with NaNs. 
% Time-normalizing approach resamples all signals to the same number of samples, 
% effectively rescaling time to a relative scale. Time-warping approach performs 
% similar time-normalization, and then performs simultaneous time-warping alignment 
% via fdasrvf.time_warp method by Tucker (2014). 

% 3) After each of the three approaches, mean and SD are extracted. The results 
% are displayed. Peaks are identified and characterized. Their properties are compared 
% against “ground truth” obtained from the distributions from (1). These distributions 
% and the estimate after alignment are displayed as histograms. The resulting 
% differences (errors of estimation from the mean) are additionally displayed as 
% horizontal and vertical barcharts. 

% 4) For the time-warping approach, estimated warping functions are plotted, 
% against with “TimeShift” where identity was subtracted, and “TimeSpeed”, a derivative 
% of warping function. Additionally, the extracted mean warped with every warping 
% function is plotted as another ensemble demonstrating contribution of the 
% temporal (warp-exposed) variability to the “spatial”, signal domain. 

% 5) Variabilities are displayed as RMSE across aligned ensembles from (2), 
% with additional temporal and temporal-to-spatial components from the time-warping approach.


% Explore the NoiseOption parameter and combinations of noise imposed on the noisy signal ensemble.

% Explore the rangeSet parameter to scale the amount of corresponding noise imposed and 
% relate that noise to post-alignment parameter estimation error and variabilities.


% Aleksei Krotov
% Northeastern University, 2024


figpos0 = [0 0 0 0]; % Change the first argument(s) if not plotting on the main monitor
%figpos0 = [-1920 0 0 0];


SampleRate = 300;
N = 500; % noisy versions

% See the bottom of the script file for my exploration of possible functions.

% Specifying T1 and T2
% fH = @(tt,Pars) Pars(1) * exp(-4.*log(2) ./ (Pars(3)).^2 .* (tt-Pars(5)).^2) +...
%     Pars(2) * exp(-4.*log(2) ./ (Pars(4)).^2 .* (tt-Pars(6)).^2);

% Adding a small third peak, related to the second peak, for demonstration of small-feature-resolving.
fH = @(tt,Pars) Pars(1) * exp(-4.*log(2) ./ (Pars(3)).^2 .* (tt-Pars(5)).^2) +...
    Pars(2) * exp(-4.*log(2) ./ (Pars(4)).^2 .* (tt-Pars(6)).^2) + ...
    Pars(2)/3 * exp(-4.*log(2) ./ (Pars(4)/5).^2 .* (tt-Pars(6)-Pars(4)/1.5).^2); % 1/8 of A2, 1/4 of W2, delayed wrt T2 by W2/4


% See inside this one to change range boundaries or default means and SDs
Pars = setParams1();
[tforig,fparamsorig] = fMeasureImit(fH,Pars',SampleRate);
% See inside this one to change range boundaries SDs
[T,randpars,rangeSet,NoiseFlags] = setParams2(fH,Pars,fparamsorig);
% rangeSet from setParams2 has linspace of SD values.
% reset it here to a single set for a single-run.
rangeSet = {1};



RndSeed = 1; % Using Twister rng with this seed. Set to "shuffle" to use current time.
NoiseOption = 11;
% 1-6 single param: Amps (1 and 2), Widths (1 and 2), Peak Times (1 and 2)
% 7 Amps only
% 8 Widths only
% 9 Peak times only
% 10 Widths and Peak times
% 11 Everything
% 12 1st peak params only
% 13 everything+spike
% 14 everything+highfreq
% 15 amps+spike
% 16 amps+hifreq


% SPECIAL NOISE
NoiseParams = {'None',0,0};
% {'HighFreqMulti',50,0.2} Sine wave, freq = 50 Hz, uniformly random phaseshift, amplitude multiplicative F * 0.2
% {'Spike',0.05,0.5} A gaussian spike, random location in between onset and
% ofset (added after trimming 5%), width at half-max is 0.05 s, amplitude is 0.5 of the peak
% *** Used examples: 
% NoiseParams = {'Spike',0.05,Pars(8)};
% NoiseParams = {'HighFreqMulti',50,Pars(8)};
if NoiseOption > 12
    disp('If adding spike or highfreq noise, parameter estimation is inaccurate!');
end



% Save the profiles as a movie sequence, framerate 2 fps
%NoiseOptionNames = {'A1','A2','B1','B2','DT','A1 and A2','B1 and B2','B1,B2,DT','Everything','Special 1','Everything med + spike','Everything med + highfreq','Amps + spike','Amps + highfreq'};
NoiseOptionNames = {'A1','A2','B1','B2','T1','T2','A1 and A2','B1 and B2','T1 and T2','B1,B2,T1,T2','Everything','A1,B1,T1','Everything med + spike','Everything med + highfreq','Amps + spike','Amps + highfreq'};

figname = 'err(sigma) Everything';
figname = strcat('err(sigma)',NoiseOptionNames{NoiseOption});
%figname = 'test';
if rangeSet{1} > 1
    if ~strcmpi(figname,'test')
        disp(sprintf('<strong>Animation is recorded, is it really necessary?</strong>'));
    end
    FrameRate = 2;
    v = VideoWriter(sprintf('EFDA analyses/Method manuscript/Synthetic2Gauss_Parameters/%s',figname),'MPEG-4');
    v.FrameRate = FrameRate;
    open(v);
end

for isi = 1:rangeSet{1}
    disp(sprintf('Trying a range of induced SDs, stage %d of %d',isi,rangeSet{1}));
    % Check sizes... or wait for bugs and for a better structuring idea
    % 1. Assign values from range (or use constant ones)
    for iparrange = 1:length(rangeSet)-1
        if any(rangeSet{iparrange+1} ~= 0) % If [0 0] range specified, leave default values
            range = linspace(rangeSet{iparrange+1}(1),rangeSet{iparrange+1}(2),rangeSet{1});
            if any(randpars(iparrange,:) ~= 0) % assign range to noise parameter=
                randpars(iparrange,1) = range(isi);
                % 2. Adjust 3sigma distribution truncation
                randpars(iparrange,4:5) = [-1 1] .* 3 .* range(isi);
            else % assign range to Par, i.e. to the fix value not noise parameter.
                Pars(iparrange,1) = range(isi);
            end
        end
    end
    % 3. Multiply with NoiseFlags
    if NoiseOption == 0
        randpars(1:6,1) = zeros(5,1);
    else
        randpars(1:6,1) = randpars(1:6,1) .* NoiseFlags(NoiseOption,:)';
    end
    % 4. Specially for noise
    NoiseParams{3} = Pars(8);
    
    %%%%%%%%%%%%%%%%%% TRYING TO MAKE ITERATIVE AMOUNT OF NOISE, NEED TO CHANGE THIS LOOP
    %%%%%%%%%%%%%%%%%% LINES ABOVE, AND CHANGE THE EFDA_GAUSS ACTION TO USE RANDPARS(end)
    %%%%%%%%%%%%%%%%%% AS NOISE AMPLITUDE.
    
    flag_plot = 1; % Keep it on!
    [T,tfmeanT,tfmeanP,tfmeanE,curFig] = EFDA_GaussAction(randpars,Pars,NoiseParams,fH,N,SampleRate,RndSeed,tforig,T,rangeSet{1},isi,figpos0,flag_plot);
    
    % When many parameters, save the figures as frames, produce a slow video, i.e. 1 s per
    % frame?
    
    % Record video of profiles at changing noise
    if rangeSet{1} > 1
        drawnow limitrate
        writeVideo(v, getframe(curFig));
    end
end



% Now visualize that T
if rangeSet{1} > 1
    Pars2Plot = 10:28;
    SigmasColForPlot = {2,5}; % If {2,3}, then a column of Sigma. If 5, then T{column 5}
    [curFigRes,curFigResVar] = PlotMetricsVsSigma(T,Pars2Plot,SigmasColForPlot,figpos0);
    SaveFig_FigPng('Handle',curFigRes,'FigName',figname,'Path','EFDA analyses/Method manuscript/Synthetic2Gauss_Parameters');
    SaveFig_FigPng('Handle',curFigResVar,'FigName',strcat(figname,' Variabilities'),'Path','EFDA analyses/Method manuscript/Synthetic2Gauss_Parameters');
    close(v); % and by the way, close the video object
end


% Visualize a trajectory?
if 0
    % Let it be a line
    curFigTrajVis = Create_Reuse_Figure([],'TestTraj',[1950,200,600,400]);
    Lorig = tforig(end,3);
    spOrig = subplot(1,4,1); hold on; title('Original'); xlabel('X (m)'); ylabel('Y (m)');
    spT = subplot(1,4,2); hold on; title('TimeNorm'); xlabel('X (m)');
    spP = subplot(1,4,3); hold on; title('PathNorm'); xlabel('X (m)');
    spE = subplot(1,4,4); hold on; title('TimeNormEFDA'); xlabel('X (m)');
    linkaxes([spOrig spT spP spE],'xy');
    
    % Find original max path to define limits and a line.
    line(spOrig,[0.5 0.5],[0 Lorig],'Color','k');
    line(spT,[0.5 0.5],[0 Lorig],'Color','k');
    line(spP,[0.5 0.5],[0 Lorig],'Color','k');
    line(spE,[0.5 0.5],[0 Lorig],'Color','k');
    sc(1) = scatter(spOrig,0.5,0,50,'k','filled','MarkerEdgeColor','none');
    sc(2) = scatter(spT,0.5,0,50,'k','filled','MarkerEdgeColor','none');
    sc(3) = scatter(spP,0.5,0,50,'k','filled','MarkerEdgeColor','none');
    sc(4) = scatter(spE,0.5,0,50,'k','filled','MarkerEdgeColor','none');
    
    
    TDur = tforig(end,1);
    % Assuming the same duration!!! Which is quite true for all the methods
    PlaySpeed = 0.2;
    FPS = 25;
    NFrames = round(FPS .* TDur ./ PlaySpeed);
    
    tfPlot(:,:,1) = interp1(linspace(0,1,length(tforig))',tforig,linspace(0,1,NFrames));
    tfPlot(:,:,2) = interp1(linspace(0,1,length(tfmeanT))',tfmeanT,linspace(0,1,NFrames));
    tfPlot(:,:,3) = interp1(linspace(0,1,length(tfmeanP))',tfmeanP,linspace(0,1,NFrames));
    tfPlot(:,:,4) = interp1(linspace(0,1,length(tfmeanE))',tfmeanE,linspace(0,1,NFrames));
    
    
    % Use mean path to obtain average speed???
    for iframe = 1:NFrames
        % Oooops, should I follow the path profile (which is almost ideal for all the cases)
        % or the speed profile (which is not ideal)
    end
end

% let it be a quarter-circle arc

end



function Pars = setParams1()
%   MEAN PARAMS
Pars = [3;... %A1
    4;... %A2
    0.16;... % B1 Full-width at half-maximum
    0.16;... % B2
    0.4;    % Delta-T time between the two peaks !!!!!!!!!!!
    0.3;... % Time of the first peak
    0.02;... % Trim level - Imitating a real measurement, set onset and offset at this % of max. In ppt currently 0.02
    0.5]; % Noise (or Spike) level - amplitude of sine noise as fraction of signal (multiplicative)
% or amplitude of a Gaussian spike in percentage of Amax


% Changing on Feb 2024, to avoid incorrect estimation of Widths
Pars = [3;... %A1
    4;... %A2
    0.14;... % B1 Full-width at half-maximum
    0.14;... % B2
    0.3;    % Time of the first peak
    0.7;... % Time of the second peak
    0.02;... % Trim level - Imitating a real measurement, set onset and offset at this % of max. In ppt currently 0.02
    0.5]; % Noise (or Spike) level - amplitude of sine noise as fraction of signal (multiplicative)
% or amplitude of a Gaussian spike in percentage of Amax


% FROM Efda concept_demo
% Pars1 = [4; 3; 0.16; 0.16; 0.4; 0.3; 0.00; 0];
% Pars2 = [2; 5; 0.16; 0.10; 0.6; 0.2; 0.00; 0];

% From here (used to be before SFN)
% [4 3 .16 .16 .4 .3 .02 0];

end



function [T,randpars,rangeSet,NoiseFlags] = setParams2(fH,Pars,fparamsorig)

% Create a table of stats
Stage = {'Init'};
Sigmas = [0 0 0 0 0 0];
SigmaScale = 0;
NoiseType = {'None'};
NoiseTemp = 0;
NoiseAmp = 0;
OnOffThres = Pars(7);
Values = {'Orig'};
Dur = fparamsorig.Dur;
A1 = fparamsorig.A1;
Amin = fparamsorig.Amin;
A2 = fparamsorig.A2;
TA1 = fparamsorig.TA1;
TAmin = fparamsorig.TAmin;
TA2 = fparamsorig.TA2;
Dt = fparamsorig.Dt;
tFWHM1 = fparamsorig.tFWHM1;
tFWHM2 = fparamsorig.tFWHM2;
SlopeHML = fparamsorig.SlopeHML;
SlopeHMR = fparamsorig.SlopeHMR;
OrigRMSE = nan;
OrigIntegralRMSE = nan;
SD_A1 = nan;
SD_Amin = nan;
SD_A2 = nan;
VarSpat = nan;
VarTemp = nan;
VarTemp2Spat = nan;
VarsOrigAll = nan(1,4); % the special cases only for the expected manuscript to compare these variability expressions.
VarsAlignedAll = nan(1,4);
VarsTShiftAll = nan(1,4);
VarsAlignedMeanWarpedAll = nan(1,4);
T = table(Stage,SigmaScale,Sigmas,NoiseType,NoiseTemp,NoiseAmp,OnOffThres,Values,...
    Dur,A1,Amin,A2,TA1,TAmin,TA2,Dt,tFWHM1,tFWHM2,SlopeHML,SlopeHMR,SD_A1,SD_Amin,SD_A2,...
    OrigRMSE,OrigIntegralRMSE,VarSpat,VarTemp,VarTemp2Spat,VarsOrigAll,VarsAlignedAll,VarsTShiftAll,VarsAlignedMeanWarpedAll);


% Distribution parameters for noisy signal variables
% First 3 cols are reserved for distribution params (only 1 typically means that it is sigma)
% The two last are min/max truncate levels of a distribution to avoid aggressive outliers
% Zeros mean that no noise is imposed, only the absolute value might be sweeped
randpars = zeros(length(Pars),5);
randpars = randpars +...
    [.3 0 0 -.9 .9;... % A1
    .3 0 0 -.9 .9;... % A2
    .04 0 0 -.12 .12;... % B1
    .04 0 0 -.12 .12;... % B2
    .03 0 0 -.09 .09;... % T of the first peak
    .03 0 0 -.09 .09;... % T of the second peak
    0 0 0 0 0;... % Trim level % of max
    0 0 0 0 0]; % Noise amplitude


%%% AT SFN WORKING FOR THE POSTER, RESET THESE PAREAMETERS
% randpars = randpars +...    
%     [.6 0 0 -1.5 1.5;... % A1
%     .6 0 0 -1.5 1.5;... % A2
%     .05 0 0 -.15 .15;... % B1
%     .05 0 0 -.15 .15;... % B2
%     .05 0 0 -.1 .1;... % dT % 0.15 for 'Strong DT'
%     .04 0 0 -.09 .09;... % t0   % USED TO BE 0 0 0 0 0;... % t0
%     0 0 0 0 0;... % Trim level % of max
%     0 0 0 0 0]; % Noise amplitude




% Sweep noise magnitude or the mean of each parameter (if NoiseFlags allow)
% 1st cell is N of linspaced values. The rest are vectors for min and max.
% If more than 1 non-zero parameter for a distribution, then corresponding vector
% would be not 2x1 but 4x1 or 6x1...

% Default ranges
rangeSet = {25; % number of linspaced values
    [eps .5]; % A1
    [eps .5]; % A2
    [eps .07]; % B1
    [eps .07]; % B2
    [eps .12]; % T1
    [eps .12]; % T2
    [0 0]; % Trim level % of max
    [0 0]}; % Noise amplitude (First try 0 to 10)
% I expect to play with time of the first peak and with trim level, but without noise

% A special set for varying Everything + Noise
if 0
    rangeSet = {25; % number of linspaced values
        [.3 .3]; % A1
        [.3 .3]; % A2
        [.05 .05]; % B1
        [.05 .05]; % B2
        [.03 .03]; % T1
        [.03 .03]; % T2
        [0 0]; % Trim level % of max
        [eps 0.7]}; % Noise amplitude (First try 0 to 10)
end

% For now only 5 parameters can be noisy (add more if needed). There're 9 options for
% simultaneous noise right now. For corresponding parameter, randpars value will be used,
% but can be overwritten with sigmaSet or sigmaOnce
NoiseFlags = [
    1 0 0 0 0 0;...
    0 1 0 0 0 0;...
    0 0 1 0 0 0;...
    0 0 0 1 0 0;...
    0 0 0 0 1 0;...
    0 0 0 0 0 1;...
    1 1 0 0 0 0;... % 7 - amplitudes only
    0 0 1 1 0 0;... % widths only
    0 0 0 0 1 1;... % peak times only
    0 0 1 1 1 1;... % temporal params only
    1 1 1 1 1 1;... % 11 - all parameters noisy
    1 0 1 0 1 0;... % 12 - Special. First peak params only
    1 1 1 1 1 1;... % Everything with spike noise
    1 1 1 1 1 1;... % Everything with highfreq noise
    1 1 0 0 0 0;... % Amps with spike noise
    1 1 0 0 0 0];   % Amps with highfreq noise

end







% Main action of the script, to use it in a loop
function [T,tfmeanT,tfmeanP,tfmeanE,curFig] = EFDA_GaussAction(randpars,MeanParams,NoiseParams,fH,N,SampleRate,RndSeed,tforig,T,rangeSig,iSig,figpos0,flag_plot)

tfmeanP = [];
tfmeanT = [];
tfmeanE = [];


[tfnoisy, fparamsnoisy,tfExamp] = GenerateNormalNoisyData(MeanParams,randpars,fH,N,SampleRate,RndSeed,'Noise',NoiseParams);

% Function T = Add data to table (fparamsnoisy)
StageName = 'Sample estimates';
T = NoisySignalsEstParams(tforig,tfnoisy,fparamsnoisy,T,StageName,MeanParams,randpars,NoiseParams,rangeSig,iSig);



% Ensure constant (or almost-constant) arclength of these noisy profiles?
% I.e. after any manipulation, scale the speed profile vertically.


% For plotting
N2Plot = 50; % limit number of lines to plot
NDownS = 200; % downsample to this value for plot performance

% OBSOLETE
% Time Norm and Path averaging
% [tfnorm, tfmean] = PathTimeNormAver(tfnoisy,SampleRate); % when path-averaging, speed is computed from path
% figname = 'Time-Norm Path-based averaging';
% PlotStuff(tforig,tfnoisy,tfnorm,tfmean,figname,N2Plot,NDownS,N);

% Time-padding
% Time Norm and Speed averaging
[tfnorm, tfmean, tfsd] = SpeedTimePadAver(tfnoisy,SampleRate); % Path is recomputed from speed
figname = 'Time-Pad Speed-based averaging';
% Prepare for multiple average columns in each figure
fignameAll = {'Signal registration and averaging demo','Time-padded'};
tfnormAll(:,1) = tfnorm;
tfmeanAll{:,1} = tfmean;
tfsdAll{:,1} = tfsd;

%PlotStuff1Row2Methods(tforig,tfnoisy,tfnorm,tfmean,figname,N2Plot,NDownS,N);
foutparams = TwoGaussParamsEstimate(tfmean(:,1:2),1,tfnorm);
StageName = 'TimePad Speed Aver';
T = TwoGaussCompareMeanOrig(tfmean,tfsd,tforig,foutparams,T,StageName,MeanParams,randpars,NoiseParams,rangeSig,iSig); % Produce a table???


% Time Norm and Speed averaging
[tfnorm, tfmean, tfsd] = SpeedTimeNormAver(tfnoisy,SampleRate); % Path is recomputed from speed
figname = 'Time-Norm Speed-based averaging';
% Prepare for multiple average columns in each figure
%fignameAll = {'Signal registration and averaging demo','Time-normalized'};
fignameAll{end+1} = 'Time-normalized';
tfnormAll(:,2) = tfnorm;
tfmeanAll{:,2} = tfmean;
tfsdAll{:,2} = tfsd;

%PlotStuff1Row2Methods(tforig,tfnoisy,tfnorm,tfmean,figname,N2Plot,NDownS,N);
foutparams = TwoGaussParamsEstimate(tfmean(:,1:2),1,tfnorm);
StageName = 'TimeNorm Speed Aver';
T = TwoGaussCompareMeanOrig(tfmean,tfsd,tforig,foutparams,T,StageName,MeanParams,randpars,NoiseParams,rangeSig,iSig); % Produce a table???
%!!! When speed is time-normalized without amplitude adjustment, the
% recomputed path changes! // Not considering the hypothetical path anymore.


% OBSOLETE
% Path Norm and Path averaging
% [tfnorm, tfmean] = PathPathNormAver(tfnoisy,SampleRate); % when path-averaging, speed is computed from path
% figname = 'Path-Norm Path-based averaging';
% PlotStuff(tforig,tfnoisy,tftimenorm,tfmean,figname,N2Plot,NDownS,N);
% This one probably doesn't work properly.
% Think carefully once again what you are averaging.

% Path Norm and Speed averaging - not showing for now
% [tfnorm, tfmean] = SpeedPathNormAver2(tfnoisy,SampleRate);
% tfmeanP = tfmean;
% figname = 'Path-Norm Speed-based averaging';
% PlotStuff3Rows(tforig,tfnoisy,tfnorm,tfmean,figname,N2Plot,NDownS,N);
% foutparams = TwoGaussParamsEstimate(tfmean(:,1:2),1);
% StageName = 'PathNorm Speed Aver';
% T = TwoGaussCompareMeanOrig(tfmean,tforig,foutparams,T,StageName,randpars);

% EFDA time-norm Speed-averaging
[tfnorm, tfmean, tfsd,ElasticParams,EFDARes] = SpeedTimeNormEFDAAver(tfnoisy,SampleRate);
figname = 'Time-Norm Speed-based EFDA averaging';

fignameAll(end+1) = {'EFDA'};
tfnormAll(:,3) = tfnorm;
tfmeanAll{:,3} = tfmean;
tfsdAll{:,3} = tfsd;


foutparams = TwoGaussParamsEstimate(tfmean(:,1:2),1,tfnorm);
StageName = 'EFDA Speed Aver';
T = TwoGaussCompareMeanOrig(tfmean,tfsd,tforig,foutparams,T,StageName,MeanParams,randpars,NoiseParams,rangeSig,iSig);
T.VarSpat(height(T)) = ElasticParams.VarSpat;
T.VarTemp(height(T)) = ElasticParams.VarTemp;
T.VarTemp2Spat(height(T)) = ElasticParams.VarTemp2Spat;

% T.VarsOrigAll(height(T),:) = ElasticParams.VarsAll.VarsOrig;
% T.VarsAlignedAll(height(T),:) = ElasticParams.VarsAll.VarsAligned;
% T.VarsTShiftAll(height(T),:) = ElasticParams.VarsAll.VarsTShift;
% T.VarsAlignedMeanWarpedAll(height(T),:) = ElasticParams.VarsAll.VarsAlignedMeanWarped;
% T.EFDA_VarTSpeed(height(T)) = ElasticParams.VarTSpeed;


%curFig = PlotStuff2Rows3MethodsHistograms(tforig,tfnoisy,tfExamp,tfnormAll,tfmeanAll,tfsdAll,fignameAll,N2Plot,NDownS,N,T,fparamsnoisy);

% if flag_plot
    curFig = PlotStuff1Row3Methods(tforig,tfnoisy,tfExamp,tfnormAll,tfmeanAll,tfsdAll,fignameAll,EFDARes,N2Plot,NDownS,N,T,fparamsnoisy,figpos0);
    % Plot warps? % And also aligned-mean-warped with the SD band
    curFigW = PlotWarps3Plots(EFDARes,figpos0);
    % Plot warps with full duration info?
% end

%disp('Not plotting histograms this time');
%return


% Histograms with param estimates
curFigHistParams = Create_Reuse_Figure([],'Parameter estimates',figpos0 + [80 650 1400 200]);
ParamsList = {'A1','A2','tFWHM1','tFWHM2','TA1','TA2'};
ParamsNames = {'A1','A2','W1','W2','T1','T2'};
MethodCols = [0.7451 0.333 0.5882; 0.2392 0.5098 0.6941; 0.9961 0.3451 0.0824];
for ipar = 1:6
    subplot(1,6,ipar); hold on
    title(ParamsNames{ipar});
    TT = struct2table(fparamsnoisy);
    histogram(TT.(ParamsList{ipar}), 'FaceColor',[.6 .6 .6],'DisplayName','Sample')
    M = mean(TT.(ParamsList{ipar}));
    SD = std(TT.(ParamsList{ipar}));
    xline(M,'Color','k','LineWidth',3,'DisplayName','Sample Mean')
    xline(M + SD,'Color','k','LineStyle','--','LineWidth',1,'DisplayName','Sample SD');
    xline(M - SD,'Color','k','LineStyle','--','LineWidth',1,'HandleVisibility','off');
    T1 = M + T.(ParamsList{ipar}); 
    xline(T1(3),'Color',MethodCols(1,:),'LineWidth',2,'DisplayName','TimePad Mean');
    xline(T1(4),'Color',MethodCols(2,:),'LineWidth',2,'DisplayName','TimeNorm Mean');
    xline(T1(5),'Color',MethodCols(3,:),'LineWidth',2,'DisplayName','TimeWarp Mean');
    if ipar == 1
        TSd = T.SD_A1(3);
        xline(T1(3) + TSd * [-1 1],'Color',MethodCols(1,:),'LineStyle','--','LineWidth',1,'HandleVisibility','off')
        TSd = T.SD_A1(4);
        xline(T1(4) + TSd * [-1 1],'Color',MethodCols(2,:),'LineStyle','--','LineWidth',1,'HandleVisibility','off')
        TSd = T.SD_A1(5);
        xline(T1(5) + TSd * [-1 1],'Color',MethodCols(3,:),'LineStyle','--','LineWidth',1,'HandleVisibility','off')
    end
    if ipar == 2
        TSd = T.SD_A2(3);
        xline(T1(3) + TSd * [-1 1],'Color',MethodCols(1,:),'LineStyle','--','LineWidth',1)
        TSd = T.SD_A2(4);
        xline(T1(4) + TSd * [-1 1],'Color',MethodCols(2,:),'LineStyle','--','LineWidth',1)
        TSd = T.SD_A2(5);
        xline(T1(5) + TSd * [-1 1],'Color',MethodCols(3,:),'LineStyle','--','LineWidth',1)
    end
    if ipar == 1
        Leg = legend('location','bestoutside');
        Leg.FontSize = 9;
        Leg.Position(1) = 0.005;
        
    end
end

%%%%%% Plot amplitude and temporal elastic distances to mean?
% curFigDist = Create_Reuse_Figure([],'Elastic distances',figpos0 + [80 600 800 300]);
% %curFigDist = Create_Reuse_Figure([],'Elastic distances',[-1000 600 800 300]); % Home
% subplot(2,3,1);
% % Fisher-Rao distance between aligned SRVFs and mean SRVF, [m/s], if original function is [m/s]
% histogram(EFDARes.ElasticSpatDist); title ('SRVF Amplitude Dist, (m/s)');
% subplot(2,3,4);
% % Fisher-Rao distance between gamma SRVFs and identity, [s]
% histogram(EFDARes.ElasticPhDist); title ('SRVF Phase Dist, (s)');
% 
% for i = 1:N
%     % L2-distance between aligned fn and mean, [m]
%     AlignedDist(i) = trapz(EFDARes.TimeMean,sqrt((EFDARes.FuncAligned(:,i) - EFDARes.FuncAlignedMean).^2));
%     % L2-distance between aligned mean and the same mean warped respectively, [m]
%     AlignedPhDist(i) = trapz(EFDARes.TimeMean,sqrt((EFDARes.FuncMeanWarped(:,i) - EFDARes.FuncAlignedMean).^2));
%     % L2-distance between TimeShift and zero, [s.^2]
%     TimeShiftDist(i) = trapz(EFDARes.TimeMean,sqrt((EFDARes.WarpsMinusUnity(:,i)).^2));
%     % L2-distance between TimeSpeed and one, [s]
%     TimeSpeedDist(i) = trapz(EFDARes.TimeMean,sqrt((EFDARes.WarpsDdt(:,i) - ones(length(EFDARes.TimeMean),1)).^2));
% end
% 
% subplot(2,3,2);
% histogram(AlignedDist); title ('L2-dist Amp, (m)');
% subplot(2,3,5);
% histogram(AlignedPhDist); title ('L2-dist Phase, (m)');
% subplot(2,3,3);
% histogram(TimeShiftDist); title ('TimeShift Dist, (s^2)');
% subplot(2,3,6);
% histogram(TimeSpeedDist); title('TimeSpeedDist, (s)');


% EFDA time-norm Path-averaging
% [tfnorm, tfmean] = SpeedPathNormEFDAAver(tfnoisy,SampleRate);
% tfmeanEP = tfmean;
% figname = 'Path-Norm Speed-based EFDA averaging';
% PlotStuff(tforig,tfnoisy,tfnorm,tfmean,figname,N2Plot,NDownS,N);
% foutparams = TwoGaussParamsEstimate(tfmean(:,1:2),1);
% StageName = 'PathNorm Speed EFDA';
% T = TwoGaussCompareMeanOrig(tfmean,tforig,foutparams,T,StageName,randpars);

end









%%% Generate random functions
% Add height adjustment to maintain the same trajectory length???
function [tfnoisy, fparamsnoisy,tfExamp] = GenerateNormalNoisyData(MeanParams,randpars,fH,N,SampleRate,RndSeed,varargin)
SampParams = zeros(N,8);
tfnoisy = cell(N,1);
fparamsnoisy = struct;


rng(RndSeed,"twister");
fprintf('Using Twister rand generator with seed %d \n', RndSeed)
% Fix from noise flag to randpars(i, 1) which is usually sigma.
for ipar = 1:6
    if randpars(ipar,1) % nonzero
        flagW = 1;
        ind = 1:size(SampParams,1);
        % Sample from a distribution. For values outside of the bounds, keep resampling.
        while flagW
            SampParams(ind,ipar) = randn(length(ind),1) .* randpars(ipar,1);
            ind = find(SampParams(:,ipar) > randpars(ipar,5) | SampParams(:,ipar) < randpars(ipar,4));
            if isempty(ind), flagW = 0; end
        end
    end
end


if any(strcmpi(varargin,'Noise'))
    indN = find(any(strcmpi(varargin,'Noise')));
    NoisePars = varargin{indN+1};
    NoiseType = NoisePars{1};
    NoiseFreq = NoisePars{2}; % And if spike, what does this mean?
    NoiseAmp = NoisePars{3};
end


NoisyParams = SampParams + ones(N,8) .* MeanParams'; %DIMENSION
for iN = 1:N
    SampParamSet = NoisyParams(iN,:);
    [tfnoisy{iN,1},fparamsnoisy0] = fMeasureImit(fH,SampParamSet,SampleRate,'Noise',{NoiseType,NoiseFreq,NoiseAmp});
    fldn = fieldnames(fparamsnoisy0);
    for ifld = 1:length(fldn)
        fparamsnoisy(iN).(fldn{ifld}) = fparamsnoisy0.(fldn{ifld});
    end
end
% Find the two most different signals from the sample to use as examples of variability
% Metric - either cross-correlation or RMS?...
% To increase speed, first resample all of them.
NSamp = 100;
tfnew = nan(NSamp,N);
for iN = 1:N
    tf = tfnoisy{iN,1};
    tfnew(:,iN) = interp1(linspace(0,1,length(tf))',tf(:,2),linspace(0,1,NSamp)');
end

Corr = nan(N,N);
RMSE = nan(N,N);
for iN = 1:N
    tf1 = tfnew(:,iN);
    Corr(iN,:) = corr(tf1,tfnew);
    RMSE(iN,:) = sqrt(sum((tf1 - tfnew).^2));
end
[CorrMin,iCorrMin] = min(Corr,[],'all','linear');
[RMSEmax,iRMSEMax] = max(RMSE,[],'all','linear');

%[r,c] = ind2sub([N N],iCorrMin);
[r,c] = ind2sub([N N],iRMSEMax);
tfExamp = tfnoisy([r c],1);

% Debug?
% figure,pc = pcolor(RMSE); pc.LineStyle = 'none';
% figure,pc = pcolor(Corr); pc.LineStyle = 'none';
% figure,scatter(RMSE(:),Corr(:));
end




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




%%%%% AVERAGING METHODS %%%%%%%%%%%

% Align beginning, pad the rest with nans, and extract parameters
function [tfnorm, tfmean, tfsd] = SpeedTimePadAver(tfnoisy,SampleRate)
tdur = nan(length(tfnoisy),1);
for itf = 1:length(tfnoisy)
    tdur(itf) = tfnoisy{itf}(end,1);
    L(itf) = length(tfnoisy{itf});
end
yres = nan(max(L),length(tfnoisy));
tmean = (1:max(L))' ./ SampleRate;
for itf = 1:length(tfnoisy)
    yres(1:length(tfnoisy{itf}),itf) = tfnoisy{itf}(:,2);
    tfnorm{itf,1} = nan(max(L),2);
    tfnorm{itf,1}(:,1) = tmean;
    tfnorm{itf,1}(1:length(tfnoisy{itf}),2) = tfnoisy{itf}(:,2);
end

tfmean = [tmean mean(yres,2,'omitnan') cumtrapz(mean(yres,2,'omitnan')) ./ SampleRate];
tfsd = [tmean std(yres,[],2,'omitnan')];
% ALTERNATIVELY, IF AVER PATH IS OBTAINED FROM NORMED PATHS
% tfmean(:,3) = mean(yresInt,2);

f = tfmean(:,2);
%foutparams = TwoGaussParamsEstimate(f,SampleRate,tfnorm);

end


% Time-normalize and average *the speed* and extract parameters
function [tfnorm, tfmean,tfsd] = SpeedTimeNormAver(tfnoisy,SampleRate)
tdur = nan(length(tfnoisy),1);
for itf = 1:length(tfnoisy)
    tdur(itf) = tfnoisy{itf}(end,1);
end
durmean = mean(tdur);
Nresamp = ceil(durmean .* SampleRate);
%tmean =  %%% Make sure times are consistent in plots.....
tmean = (1:Nresamp)' ./ SampleRate;
yres = nan(Nresamp,length(tfnoisy));
for itf = 1:length(tfnoisy)
    y = tfnoisy{itf}(:,2);
    t = tfnoisy{itf}(:,1);
    tnew = linspace(t(1),t(end),Nresamp)';
    yres(:,itf) = interp1(t,y,tnew);
    % ALTERNATIVELY Adjust amplitude to maintain path?
    % yres(:,itf) = yres(:,itf) ./ tmean .* t(end);
    
    yresInt(:,itf) = cumtrapz(yres(:,itf)) ./ SampleRate;
    tfnorm{itf,1} = [tmean yres(:,itf) yresInt(:,itf)];
end
tfmean = [tmean mean(yres,2) cumtrapz(mean(yres,2)) ./ SampleRate];
tfsd = [tmean std(yres,[],2,'omitnan')];

% ALTERNATIVELY, IF AVER PATH IS OBTAINED FROM NORMED PATHS
% tfmean(:,3) = mean(yresInt,2);

f = tfmean(:,2);
%foutparams = TwoGaussParamsEstimate(f,SampleRate,tfnorm);

end





%function [yaver,ystd,ResultStruct] = AverEFDA(t,y,tmax,Amax,tHW,yorig,N,TDur,SampleRate)
% Time-normalize, EFDA-align speeds, average speed
function [tfnorm, tfmean, tfsd, ElasticParams,ResultStruct] = SpeedTimeNormEFDAAver(tfnoisy,SampleRate)

for itf = 1:length(tfnoisy)
    y{itf} = tfnoisy{itf}(:,[1 2]);
    dur(itf) = tfnoisy{itf}(end,1);
end
NSampEFDA = 100;



% Newer version of EFDA alignment
opts.FSResampleTarget = 50;
opts.EFDAParallel = 1;
opts.EFDAWarpsWithDurs = 0;
opts.EFDALambda = 0.01;
opts.EstimateExecutionTime = 1;
opts.EFDAGraphics = 0;

[~, ResultStruct] = EFDA_alignmentPublishing(y',opts);
NSampEFDA = length(ResultStruct.TimeMean);




durmean = mean(dur);
f = nan(NSampEFDA,length(tfnoisy));
tfnorm = cell(length(tfnoisy),1);
for itf = 1:length(tfnoisy)
    gamma = ResultStruct.Warps(:,itf);
    t = ResultStruct.TimeMean .* durmean; % For the display here, scale to mean time.
    ind = linspace(0,1,length(tfnoisy{itf}));
    indNew = linspace(0,1,NSampEFDA);
    f(:,itf) = interp1(ind,tfnoisy{itf}(:,3),indNew);
    tfnorm{itf} = [t ResultStruct.FuncAligned(:,itf) warp_f_gamma(f(:,itf),gamma,t)'];
end
vmean = ResultStruct.FuncAlignedMean;
vsd = ResultStruct.FuncAlignedStd;
t = ResultStruct.TimeMean .* durmean; % For the display here, scale to mean time.
pathmean = mean(f,2);

tfmean = [t vmean pathmean];
tfsd = [t vsd];

ElasticParams.ElasticSpat = ResultStruct.ElasticSpatDist;
ElasticParams.ElasticPhas = ResultStruct.ElasticPhDist;

% The three proposed ones
ElasticParams.VarSpat = ResultStruct.VarSpat;
ElasticParams.VarTemp = ResultStruct.VarTemp;
ElasticParams.VarTemp2Spat = ResultStruct.VarTemp2Spat;

% Many others
% ElasticParams.VarsAll = ResultStruct.VarsAll;

% ElasticParams.VarTShift = ResultStruct.VariabilityTimeShift;
% ElasticParams.VarTSpeed = ResultStruct.VariabilityTimeSpeed;

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
% [ilocs,indsort] = sort(ilocs);
% pks = pks(indsort);
% w = w(indsort);
% p = p(indsort);

% [p, ipk] = maxk(p,2);
% pks = pks(ipk);
% ilocs = ilocs(ipk);
% w = w(ipk);

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
    
    %foutparams.VarSpat = trapz(var(y,[],2,'omitnan')) ./ size(y,1) ./ timemean(end);
    
    % Fix from Jan 2023. Proposed var is the area of the SD(t), normalized by average time
    % duration

    foutparams.VarSpat = trapz(std(y,[],2,'omitnan')) ./ size(y,1) ./ mean(timemean,'omitnan');
    % But for TimePad, additionally divide over total time.
    % add srrf var amp
    


    % warning('Add calculation of SRRFs and compute variability over them');
    
else
    
    foutparams.VarSpat = nan;
end


end
%


function T = NoisySignalsEstParams(tforig,tfnoisy,foutparams,T,StageName,MeanParams,randpars,NoiseParams,rangeSig,iSig)
warning('off','MATLAB:table:RowsAddedExistingVars');
indR = height(T) + 1;
T.Stage(indR) = {StageName};
T.SigmaScale(indR) = (iSig-1) ./ (rangeSig-1);
T.Sigmas(indR,:) = randpars(1:6,1)';
T.Values(indR) = {'Mean of Est-s'};
T.OnOffThres(indR) = MeanParams(7);
if NoiseParams{3} == 0
    T.NoiseType(indR) = {'None'};
    T.NoiseTemp(indR) = 0;
    T.NoiseAmp(indR) = 0;
else
    T.NoiseType(indR) = NoiseParams(1);
    T.NoiseTemp(indR) = NoiseParams{2};
    T.NoiseAmp(indR) = NoiseParams{3};
end
NameList = {'Dur','A1','Amin','A2','TA1','TAmin','TA2','Dt','tFWHM1','tFWHM2','SlopeHML','SlopeHMR'};
for ipar = 1:length(NameList)
    T.(NameList{ipar})(indR) = mean([foutparams.(NameList{ipar})],'omitnan');
end
T.SD_A1(indR) = std([foutparams.A1],[],'omitnan');
T.SD_A2(indR) = std([foutparams.A2],[],'omitnan');
T.SD_Amin(indR) = std([foutparams.Amin],[],'omitnan');

T.OrigRMSE(indR) = nan;
T.OrigIntegralRMSE(indR) = nan;
T.VarSpat(indR) = nan;
T.EFDA_ElVarPhas(indR) = nan;
T.EFDA_VarTShift(indR) = nan;
T.EFDA_VarTSpeed(indR) = nan;


warning('on','MATLAB:table:RowsAddedExistingVars');
end



% Compare metrics to original, add to the table
function T = TwoGaussCompareMeanOrig(tfmean,tfsd,tforig,foutparams,T,StageName,MeanParams,randpars,NoiseParams,rangeSig,iSig)
warning('off','MATLAB:table:RowsAddedExistingVars');
indR = height(T) + 1;
T.Stage(indR) = {StageName};
T.SigmaScale(indR) = (iSig-1) ./ (rangeSig-1);
T.Sigmas(indR,:) = randpars(1:6,1)';
T.Values(indR) = {'Est Err'};
T.Dur(indR) = foutparams.Dur;
T.OnOffThres(indR) = MeanParams(7);
if NoiseParams{3} == 0
    T.NoiseType(indR) = {'None'};
    T.NoiseTemp(indR) = 0;
    T.NoiseAmp(indR) = 0;
else
    T.NoiseType(indR) = NoiseParams(1);
    T.NoiseTemp(indR) = NoiseParams{2};
    T.NoiseAmp(indR) = NoiseParams{3};
end
NameList = {'Dur','A1','Amin','A2','TA1','TAmin','TA2','Dt','tFWHM1','tFWHM2','SlopeHML','SlopeHMR'};
for ipar = 1:length(NameList)
    T.(NameList{ipar})(indR) = foutparams.(NameList{ipar});
end

iA1 = round(foutparams.TA1 ./ foutparams.Dur .* length(tfmean));
iA2 = round(foutparams.TA2 ./ foutparams.Dur .* length(tfmean));
iAmin = round(foutparams.TAmin ./ foutparams.Dur .* length(tfmean));

T.SD_A1(indR) = tfsd(iA1,2);
T.SD_A2(indR) = tfsd(iA2,2);
T.SD_Amin(indR) = tfsd(iAmin,2);


paramsRange = 9:23;
% Make sure, when sweeping through a range of variabilities, that estimate Errors are with respect to a correct Mean of Estimates
indBase = find(strcmpi(T.Stage,'Sample estimates'),1,'last');
T{indR,paramsRange} = T{indR,paramsRange} - T{indBase,paramsRange}; % Absolute error


% T{indR,paramsRange} = T{indR,paramsRange} ./ T{indBase,paramsRange} - 1; % Relative Error
% warning('Relative estimation error computed and plotted');


% Resample orig to match mean to find RMSE
NOrig = length(tforig);
NNew = length(tfmean);

NOrig = length(tforig);
NMean = length(tfmean);

indOrig = linspace(0,1,NOrig)';
indNew = linspace(0,1,NNew)';
Vnew = interp1(indNew,tfmean(:,2),indOrig);
%Pathnew = interp1(indOrig,tforig(:,3),indNew);

T.OrigRMSE(indR) = sqrt(sum((Vnew - tforig(:,2)).^2) ./ length(Vnew));
%T.OrigIntegralRMSE(indR) = sum((Pathnew - tfmean(:,3)).^2) ./ length(Vnew);

T.VarSpat(indR) = foutparams.VarSpat;
T.EFDA_ElVarPhas(indR) = nan;
T.EFDA_VarTShift(indR) = nan;
T.EFDA_VarTSpeed(indR) = nan;


warning('on','MATLAB:table:RowsAddedExistingVars');
end









function PlotStuff3Rows(tforig,tfnoisy,tftimenorm,tfmean,figname,N2Plot,NDownS,N)

curFig = Create_Reuse_Figure([],figname,figpos0 + [30,100,1200,800]);
%curFig = Create_Reuse_Figure([],figname,[-1850,100,1200,800]); % LAPTOP & HOME

% Plot middle row, speed vs T
spOrig = subplot(3,3,4); hold on; xlabel('Time (s)','fontweight','bold'); ylabel('Speed (m/s)','fontweight','bold');
line(tforig(:,1),tforig(:,2),'Color','k');

spNoisy = subplot(3,3,5); hold on; xlabel('Time (s)','fontweight','bold');
for iN = 1:min(N2Plot,N)
    tf = tfnoisy{iN};
    [~,indun] = unique(tf(:,1));
    tnew = linspace(tf(1,1),tf(end,1),NDownS)';
    tfPlot = interp1(tf(indun,1),tf(indun,2),tnew);
    line(tnew,tfPlot,'Color',[.6 .6 .6],'LineWidth',1);
end
line(tforig(:,1),tforig(:,2),'Color','k');

spTimeNorm = subplot(3,3,6); hold on; xlabel('Time (s)','fontweight','bold');
for iN = 1:min(N2Plot,N)
    tf = tftimenorm{iN};
    [~,indun] = unique(tf(:,1));
    tnew = linspace(tf(1,1),tf(end,1),NDownS)';
    tfPlot = interp1(tf(indun,1),tf(indun,2),tnew);
    line(tnew,tfPlot,'Color',[.6 .6 .6],'LineWidth',1);
end
line(tforig(:,1),tforig(:,2),'Color','k');
line(tfmean(:,1),tfmean(:,2),'Color','r');


% plot top row, arclength vs time
spOrigAL = subplot(3,3,1);hold on; title('Original'); ylabel('Path (m)','fontweight','bold');
line(tforig(:,1),tforig(:,3),'Color','k');

spNoisyAL = subplot(3,3,2);hold on; title('Noisy');
for iN = 1:min(N2Plot,N)
    tf = tfnoisy{iN};
    [~,indun] = unique(tf(:,1));
    tnew = linspace(tf(1,1),tf(end,1),NDownS)';
    tfPlot = interp1(tf(indun,1),tf(indun,3),tnew);
    line(tnew,tfPlot,'Color',[.6 .6 .6],'LineWidth',1);
end
line(tforig(:,1),tforig(:,3),'Color','k');

spTimeNormAL = subplot(3,3,3); hold on; title('Normalized');
for iN = 1:min(N2Plot,N)
    tf = tftimenorm{iN};
    [~,indun] = unique(tf(:,1));
    tnew = linspace(tf(1,1),tf(end,1),NDownS)';
    tfPlot = interp1(tf(indun,1),tf(indun,3),tnew);
    line(tnew,tfPlot,'Color',[.6 .6 .6],'LineWidth',1);
end
line(tforig(:,1),tforig(:,3),'Color','k');
line(tfmean(:,1),tfmean(:,3),'Color','r');

% plot bottom row, speed vs arclength
spOrigTAL = subplot(3,3,7);hold on;  xlabel('Path (m)','fontweight','bold'); ylabel('Speed (m/s)','fontweight','bold');
line(tforig(:,3),tforig(:,2),'Color','k');

spNoisyTAL = subplot(3,3,8);hold on;  xlabel('Path (m)','fontweight','bold');
for iN = 1:min(N2Plot,N)
    tf = tfnoisy{iN};
    [~,indun] = unique(tf(:,1));
    tnew = linspace(tf(1,1),tf(end,1),NDownS)';
    tfPlot = interp1(tf(indun,1),tf(indun,2),tnew);
    tfALPlot = interp1(tf(indun,1),tf(indun,3),tnew);
    line(tfALPlot,tfPlot,'Color',[.6 .6 .6],'LineWidth',1);
end
line(tforig(:,3),tforig(:,2),'Color','k');

spTimeNormTAL = subplot(3,3,9);hold on;  xlabel('Path (m)','fontweight','bold');
for iN = 1:min(N2Plot,N)
    tf = tftimenorm{iN};
    [~,indun] = unique(tf(:,1));
    tnew = linspace(tf(1,1),tf(end,1),NDownS)';
    tfPlot = interp1(tf(indun,1),tf(indun,2),tnew);
    tfALPlot = interp1(tf(indun,1),tf(indun,3),tnew);
    line(tfALPlot,tfPlot,'Color',[.6 .6 .6],'LineWidth',1);
end
line(tforig(:,3),tforig(:,2),'Color','k');
line(tfmean(:,3),tfmean(:,2),'Color','r');


linkaxes([spOrig spNoisy spTimeNorm spOrigAL spNoisyAL spTimeNormAL],'x');
%xlim(spOrig,[0 .6]);
linkaxes([spOrigTAL spNoisyTAL spTimeNormTAL],'x');
linkaxes([spOrig spNoisy spTimeNorm],'y');
linkaxes([spOrigAL spNoisyAL spTimeNormAL],'y');
linkaxes([spOrigTAL spNoisyTAL spTimeNormTAL],'y');

end




% Only a row, without histograms. Working backup.
function curFig = PlotStuff1Row3Methods(tforig,tfnoisy,tfExamp,tftimenorm,tfmean,tfsd,figname,EFDARes,N2Plot,NDownS,N,T,fparamsnoisy,figpos0)

% If figname is a cell, the first one is the actual figname, the others are
% names of averaging columns.
% Check if tftimenorm and tfmean are cell arrays. From dimensions,
% determine subplot layout
title_flag = 1; % When preparing stacked multi-row figure
flag_ManyCols = 0;
NCols = 3;
NAver2Plot = 1;
if iscell(figname)
    flag_ManyCols = 1;
    fignameTitle = figname{1};
    NAver2Plot = length(figname) - 1;
    NCols = 3 + length(figname) - 1;
else
    fignameTitle = figname;
end

Cols = [0.7451 0.333 0.5882; 0.2392 0.5098 0.6941; 0.9961 0.3451 0.0824]; %TimePad TimeNorm and TimeWarp

curFig = Create_Reuse_Figure([],fignameTitle,figpos0 + [30 100 1000 150]);
%curFig = Create_Reuse_Figure([],fignameTitle,[-1850,100,900,200]); % LAPTOP & HOME


os.SpacingsOut = [20 30 50 40];
os.SpacingsBetween = [0 40];
os.YLabels = 'None';
os.XLabels = 'None'; % May change to "First" or to "FirstEvery2nd" or "BotLeftCorn" or "None"
os.XTickMarks = 'First';
os.YTickMarks = 'First';
[sp, osTop] = TiledAxesMy(curFig,1,NCols,[0.0 0 1 1],os);



% Plot speed vs T
% Two first columns are the same
%spOrig = subplot(1,NCols,1); hold on;
subplot(sp(1));
xlabel('Time (s)','fontweight','bold'); ylabel('f(t) (m)','fontweight','bold');
if title_flag, title('Original'); end
line(tforig(:,1),tforig(:,2),'Color','k');

subplot(sp(2));
if title_flag,  title('Noisy example'); end
line(tforig(:,1),tforig(:,2),'Color','k');
line(tfExamp{1}(:,1),tfExamp{1}(:,2),'Color',[.2 .6 .2],'LineWidth',1.5);
line(tfExamp{2}(:,1),tfExamp{2}(:,2),'Color',[.2 .2 .6],'LineWidth',1.5);

%spNoisy = subplot(1,NCols,2); hold on;
subplot(sp(3));
if title_flag, title('Noisy all'); end %xlabel('Time (s)','fontweight','bold');
for iN = 1:min(N2Plot,N)
    tf = tfnoisy{iN};
    [~,indun] = unique(tf(:,1));
    tnew = linspace(tf(1,1),tf(end,1),NDownS)';
    tfPlot = interp1(tf(indun,1),tf(indun,2),tnew);
    line(tnew,tfPlot,'Color',[.6 .6 .6],'LineWidth',.5);
end
line(tforig(:,1),tforig(:,2),'Color','k');



for icol = 1:NAver2Plot
    %spNormAver(icol) = subplot(1,NCols,2+icol); hold on; %xlabel('Time (s)','fontweight','bold');
    subplot(sp(3+icol));
    if flag_ManyCols
        if title_flag, title(figname(1+icol)); end
    end
    for iN = 1:min(N2Plot,N)
        tf = tftimenorm{iN,icol};
        [~,indun] = unique(tf(:,1));
        tnew = linspace(tf(1,1),tf(end,1),NDownS)';
        tfPlot = interp1(tf(indun,1),tf(indun,2),tnew);
        line(tnew,tfPlot,'Color',[.6 .6 .6],'LineWidth',1);
    end
    line(tforig(:,1),tforig(:,2),'Color','k');
    XPatch = [tfmean{icol}(:,1); flip(tfmean{icol}(:,1))];
    M = tfmean{icol}(:,2);
    SD = tfsd{icol}(:,2);
    [~,locs,~,~] = findpeaks(M);
    SDA(icol,1) = SD(locs(1));
    YPatch = [M+SD; flip(M-SD)];
    patch('XData',XPatch,'YData',YPatch,'FaceColor',Cols(icol,:),'FaceAlpha',0.4,'EdgeColor','none','DisplayName','SD');
    line(tfmean{icol}(:,1),M,'Color',Cols(icol,:));
    
    %line(tfmean{icol}(:,1),tfmean{icol}(:,2),'Color','r');
    
    if any(strcmpi(figname{1+icol},{'Time-normalized','Time-norm','Norm Time', 'norm time'}))
        xlabel('Time mean (s)','fontweight','bold');
    end
end
linkaxes(sp,'xy');
ylim(sp,[0 8]); % To have uniform figures.

%N=6;
% X = linspace(0,pi*3,1000); 
% Y = bsxfun(@(x,n)sin(x+2*n*pi/N), X.', 1:N); 
% C = linspecer(N);
% axes('NextPlot','replacechildren', 'ColorOrder',C);
% plot(X,Y,'linewidth',5)
% ylim([-1.1 1.1]);

% Plot barplots with estimates
%%% OLD
if 0
curFigBar = Create_Reuse_Figure([],strcat(fignameTitle,'Errors'),figpos0 + [100 300 300 400]);
spBar = subplot(1,1,1);
Pars2Plot = [10 11 12 21 22 23 13 14 15 16 17 18 24];
ParsNames = {'A1','Amin','A2','SD1','SDmin','SD2','T1','Tmin','T2','DT','HW1','HW2','RMSEOrig'};
ColMap = linspecer(length(Pars2Plot));
Pars = T{3:5,Pars2Plot};
colororder(ColMap);
bp = bar(Pars);
legend(ParsNames,'location','eastoutside');
spBar.XTickLabel = {'NanPad','TimeNorm','EFDA'};
ylabel('Error','fontweight','bold');
end


AverMethods = {'TimePad','TimeNorm','EFDA'};

if 1
    curFigBar = Create_Reuse_Figure([],strcat(fignameTitle,'Param estimation errors'),figpos0 + [100 300 500 400]);
    for imethod = 1:3
        spBar(imethod) = subplot(1,3,imethod); % POSITION!!
        hold on;

        Pars2Plot = [10 11 12 21 22 23 13 14 15 16 17 18 24];
        ParNames = {'A1','Amin','A2','SD1','SDmin','SD2','T1','Tmin','T2','DT','HW1','HW2','RMSEOrig'};
        
        ColMain = Cols(imethod,:);
        Pars = T{2+imethod,Pars2Plot};
        bp(imethod) = bar(Pars);
        title(AverMethods{imethod});

        bp(imethod).Horizontal = 'on';
        spBar(imethod).YTick = 1:length(Pars2Plot);
        spBar(imethod).YTickLabel = ParNames;
        if imethod > 1
            spBar(imethod).YTickLabel = {};
            yticks(spBar(imethod),'auto');
        end

        spBar(imethod).YDir = 'reverse';
        xlim([-1.5 1.5]);
        % xlim([-0.4 0.4]);

        bp(imethod).FaceColor = ColMain;
        bp(imethod).EdgeColor = 'none';
    end
end

% Make a horizontal version as well, for the manuscript.
if 1
    curFigBarHoriz = Create_Reuse_Figure([],strcat(fignameTitle,'Param estimation errors Horiz'),figpos0 + [0 300 300 600]);
    for imethod = 1:3
        spBarH(imethod) = subplot(3,1,imethod); % POSITION!!
        hold on;

        Pars2Plot = [10 11 12 21 22 23 13 14 15 16 17 18 24];
        ParNames = {'A1','Amin','A2','SD1','SDmin','SD2','T1','Tmin','T2','DT','HW1','HW2','RMSEOrig'};
        
        ColMain = Cols(imethod,:);
        Pars = T{2+imethod,Pars2Plot};
        bp(imethod) = bar(Pars);
        title(AverMethods{imethod});

        %bp(imethod).Horizontal = 'on';
%         if imethod < 3
%             spBar(imethod).XTickLabel = {};
%             xticks(spBar(imethod),'auto');
%         end
        spBarH(imethod).XTick = 1:length(Pars2Plot);
        spBarH(imethod).XTickLabel = ParNames;

        %spBar(imethod).YDir = 'reverse';
        ylim([-1.5 1.5]);
        % xlim([-0.4 0.4]);
        linkaxes(spBarH,'y');

        bp(imethod).FaceColor = ColMain;
        bp(imethod).EdgeColor = 'none';
    end
end




% Plot variances
if 0
    curFigVarBar = Create_Reuse_Figure([],strcat(fignameTitle,' Variances'),figpos0 + [400 300 300 300]);
    spVar = subplot('Position',[.25 .2 .6 .65]); hold on
    
    
    %%%%%% CHANGE THIS ONE, USE ONLY 3 VARIABILITIES FROM EFDA
    % EFDARes.Var....
    
    
    
    ind3 = find(strcmpi(T.Stage,'EFDA Speed Aver'));
    % These T.Vars used to be computed differently. Now standardized, use VarSpat, VarTemp,
    % and
    VarsAll(:,:,1) = T.VarsOrigAll;
    VarsAll(:,:,2) = T.VarsAlignedAll;
    VarsAll(:,:,3) = T.VarsTShiftAll;
    VarsAll(:,:,4) = T.VarsAlignedMeanWarpedAll;
    VarsAll = VarsAll(ind3,:,:);
    
    % Plot for different methods for the variability in each plot of the first row
    VarMethods = {'Area of Var(t)','Area of SD(t)','MSE of L2-distances to the mean','RMSE of L2-distances to the mean'};
    VarTypes = {'Orig spat','Aligned spat','TimeShift temp','Temp2Spat'};
    
    % Which method to plot? Pick from the VarMethods list above
    Method2Plot = 2;
    Vars = [T.VarSpat(3); squeeze(VarsAll(:,Method2Plot,:))];
    bpVar = bar(Vars); 
    bpVar.FaceColor = [0 0.5 0.5];
    bpVar.XData = [1 3 5 6 7];
    title(VarMethods{Method2Plot});
    
    % yyaxis right
    % bpVarT = bar(Vars([3]));
    % bpVarT.XData = 5;
    spVar.XTick = [1 3 5 6 7];
    spVar.XTickLabels = {'TimePad Sp', 'TimeNorm Sp','TimeWarp Sp','TimeWarp T2Sp','TimeWarp Temp'};
   
end



if 1
    curFigVarBar = Create_Reuse_Figure([],strcat(fignameTitle,' Variances'),figpos0 + [400 300 300 300]);
    spVar = subplot('Position',[.25 .2 .6 .65]); hold on

    VarMethods = {'Spatial TimePad','Spatial TimeNorm','Spatial TimeWarp','Temp2Spat', 'Temporal TimeWarp'};
    Vars = [T.VarSpat(3:5)' T.VarTemp2Spat(5) T.VarTemp(5)];
    bpVar = bar(Vars); 
    bpVar.FaceColor = [0 0.5 0.5];
    bpVar.XData = [1 3 5 6 7];
    title('Variabilities');
    spVar.XTick = [1 3 5 6 7];
    spVar.XTickLabels = VarMethods;

end

end




% With histograms
function curFig = PlotStuff2Rows3MethodsHistograms(tforig,tfnoisy,tfExamp,tfnorm,tfmean,tfsd,figname,N2Plot,NDownS,N,T,fparamsnoisy)


% If figname is a cell, the first one is the actual figname, the others are
% names of averaging columns.
% Check if tftimenorm and tfmean are cell arrays. From dimensions,
% determine subplot layout
flag_ManyCols = 0;
NCols = 3;
NAver2Plot = 1;
if iscell(figname)
    flag_ManyCols = 1;
    fignameTitle = figname{1};
    NAver2Plot = length(figname) - 1;
    NCols = 3 + length(figname) - 1;
else
    fignameTitle = figname;
end

curFig = Create_Reuse_Figure([],fignameTitle,figpos0 + [30,100,1200,400]);
%curFig = Create_Reuse_Figure([],fignameTitle,[-1850,100,900,200]); % LAPTOP & HOME

%
os.SpacingsOut = [20 30 50 0];
os.SpacingsBetween = [0 40];
os.YLabels = 'None';
os.XLabels = 'BotLeftCorn'; % May change to "First" or to "FirstEvery2nd" or "BotLeftCorn"
os.XTickMarks = 'All';
os.YTickMarks = 'All';
[sp, osTop] = TiledAxesMy(curFig,1,NCols,[0.0 0.6 1 0.4],os);

% Plot speed vs T
% Three first columns are the same
%spOrig = subplot(1,NCols,1); hold on;
subplot(sp(1));
xlabel('Time (s)','fontweight','bold'); ylabel('f(t) (m)','fontweight','bold'); title('Original');
line(tforig(:,1),tforig(:,2),'Color','k');

subplot(sp(2));
title('Noisy example');
line(tforig(:,1),tforig(:,2),'Color','k');
line(tfExamp{1}(:,1),tfExamp{1}(:,2),'Color',[.2 .6 .2],'LineWidth',1.5);
line(tfExamp{2}(:,1),tfExamp{2}(:,2),'Color',[.2 .2 .6],'LineWidth',1.5);

%spNoisy = subplot(1,NCols,2); hold on;
subplot(sp(3));
title('Noisy all'); %xlabel('Time (s)','fontweight','bold');
for iN = 1:min(N2Plot,N)
    tf = tfnoisy{iN};
    [~,indun] = unique(tf(:,1));
    tnew = linspace(tf(1,1),tf(end,1),NDownS)';
    tfPlot = interp1(tf(indun,1),tf(indun,2),tnew);
    line(tnew,tfPlot,'Color',[.6 .6 .6],'LineWidth',1);
end
line(tforig(:,1),tforig(:,2),'Color','k');



for icol = 1:NAver2Plot
    subplot(sp(3+icol));
    %spNormAver(icol) = subplot(1,NCols,2+icol); hold on; %xlabel('Time (s)','fontweight','bold');
    if flag_ManyCols
        title(figname(1+icol));
    end
    for iN = 1:min(N2Plot,N)
        tf = tfnorm{iN,icol};
        [~,indun] = unique(tf(:,1));
        tnew = linspace(tf(1,1),tf(end,1),NDownS)';
        tfPlot = interp1(tf(indun,1),tf(indun,2),tnew);
        line(tnew,tfPlot,'Color',[.6 .6 .6],'LineWidth',1);
    end
    line(tforig(:,1),tforig(:,2),'Color','k');
    XPatch = [tfmean{icol}(:,1); flip(tfmean{icol}(:,1))];
    M = tfmean{icol}(:,2);
    SD = tfsd{icol}(:,2);
    [~,locs,~,~] = findpeaks(M);
    SDA(icol,1) = SD(locs(1));
    YPatch = [M+SD; flip(M-SD)];
    patch('XData',XPatch,'YData',YPatch,'FaceColor',[0.8 0 0],'FaceAlpha',0.4,'EdgeColor','none','DisplayName','SD');
    line(tfmean{icol}(:,1),M,'Color','r');
    
end
linkaxes([sp],'xy');
ylim(sp,[0 6]); % To have uniform figures.
spLines = sp;

% Now using position data from osTop, position histograms accordingly

os.SpacingsOut = [0 25 20 100];
os.SpacingsBetween = [20 50];
os.YLabels = 'None';
os.XLabels = 'First'; % May change to "First" or to "FirstEvery2nd" or "BotLeftCorn"
os.XTickMarks = 'All';
os.YTickMarks = 'None';
TNoisy = struct2table(fparamsnoisy);

Params2Plot = {'A1','TA1','Amin'};
ParamUnits = {' (m)',' (s)',' (m)'};
flag_legend = 0;
%ParamLimits
for icol = 1:NAver2Plot
    histpos = [osTop.AxPositions(1,icol+2) 0 osTop.AxesWidth 0.55];
    [sp, osOut] = TiledAxesMy(curFig,1,3,histpos,os);
    spHistAll(icol,:) = sp;
    BaseColor = [0.5 0.5 0.5];
    % Now rotate axes and plot parameters.
    for isp = 1:length(sp)
        sp(isp).View = [90,90];
        sp(isp).XDir ='reverse';
        ParamName = Params2Plot{isp};
        xlabel(sp(isp),sprintf('%s%s',ParamName,ParamUnits{isp}),'fontweight','bold');
        sp(isp).XLabel.FontSize = 10;
        histogram(sp(isp),TNoisy.(ParamName),6,'FaceColor',BaseColor,'FaceAlpha',1,'EdgeColor','none','DisplayName','Original data');
        xline(sp(isp),mean(TNoisy.(ParamName),'omitnan'),'Color','k','LineWidth',2,'DisplayName','Mean of estimates');
        xline(sp(isp),T.(ParamName)(1) + T.(ParamName)(1+icol),'Color','r','LineWidth',2,'DisplayName','Estimate of mean');
        %xlim(sp(isp),ParamLims(isp,:));
        % Add numbers
        if isp == 1
            text(sp(1),-0.7, -0.2,sprintf('SD(A)_{orig} = %.3f',T.Sigmas(1+icol,1)),'Units','normalized','FontSize',8,'fontweight','bold');
            
            text(sp(1),-0.7, -0.4,sprintf('SD(A)_{aver} = %.3f',SDA(icol,1)),'Units','normalized','FontSize',8,'fontweight','bold');
            if isnan(T.EFDA_ElVarPhas(icol+1))
                text(sp(1),3.4, -0.17,sprintf('<Var> = %.3g',T.VarSpat(icol+1)),'Units','normalized','FontSize',8,'fontweight','bold');
            else
                text(sp(1),3.4, -0.2,sprintf('<Var_{Amp}> = %.3g',T.VarSpat(icol+1)),'Units','normalized','FontSize',8,'fontweight','bold');
                text(sp(1),3.4, -0.4,sprintf('<Var_{Temp}> = %.3g',T.EFDA_ElVarPhas(icol+1)),'Units','normalized','FontSize',8,'fontweight','bold');
                text(sp(1),3.4, -0.6,sprintf('<Var_{TShift}> = %.3g',T.EFDA_VarTShift(icol+1)),'Units','normalized','FontSize',8,'fontweight','bold');
                text(sp(1),3.4, -0.8,sprintf('<Var_{TSpeed}> = %.3g',T.EFDA_VarTSpeed(icol+1)),'Units','normalized','FontSize',8,'fontweight','bold');
            end
            
            if flag_legend
                lg = legend(sp(isp),'Orientation','horizontal','location','southoutside');
                lg.NumColumns = 2;
                lg.Position(1:2) = [0.34 0.56];
                lg.FontSize = 8;
            end
        end
        
        % PRINT A TABLE ROW IN THE COMMAND LINE WITH MY PARAMETERS
        %         dPar = ParamsOrig.(ParamName) - ParamsAver.(ParamName);
        %         [~,p,ci,stats] = ttest(dPar);
        %         dParMPerc(isp) = mean(dPar) ./ mean(ParamsOrig.(ParamName)) .* 100;
        %         dParMPVal(isp) = p;
    end
    
    
    
end
for isp = 1:length(sp)
    linkaxes(spHistAll(:,isp),'x');
end

end




function curFigW = PlotWarps3Plots(EFDARes,figpos0)
curFigW = Create_Reuse_Figure([],'Warps',figpos0 + [1000,100,180,400]);
%curFigW = Create_Reuse_Figure([],'Warps',[-1850,100,180,400]); % LAPTOP & HOME
os.SpacingsOut = [10 20 60 50];
os.SpacingsBetween = [30 0];
os.YLabels = 'None';
os.XLabels = 'None'; % May change to "First" or to "FirstEvery2nd" or "BotLeftCorn" or "None"
os.XTickMarks = 'First';
os.YTickMarks = 'First';
[sp, osTop] = TiledAxesMy(curFigW,3,1,[0 0 1 1],os);
subplot(sp(1));
ylabel('Warps','fontweight','bold');
plot(EFDARes.TimeMean,EFDARes.Warps,'LineWidth',1,'Color',[.6 .6 .6]);
subplot(sp(2));
ylabel('TimeShift','fontweight','bold');
plot(EFDARes.TimeMean,EFDARes.WarpsMinusUnity,'LineWidth',1,'Color',[.6 .6 .6]);
subplot(sp(3));
ylabel('TimeSpeed','fontweight','bold');
plot(EFDARes.TimeMean,EFDARes.WarpsDdt,'LineWidth',1,'Color',[.6 .6 .6]);
xlabel('Norm Time (a.u.)','fontweight','bold');

%% For demo and manuscript
SD = std(EFDARes.WarpsMinusUnity,[],2);
M = mean(EFDARes.WarpsMinusUnity,2);
XPatch = [EFDARes.TimeMean; flipud(EFDARes.TimeMean)];
YPatch = [M + SD; flipud(M-SD)];
patch(sp(2),'XData',XPatch,'YData',YPatch,'EdgeColor','none','FaceColor','r','FaceAlpha',.4);

curFigW2 = Create_Reuse_Figure([],'WarpedFn',figpos0 + [1000,300,300,200]);
sp(4) = subplot(1,1,1); hold on; title({'Temp2Spat: Aligned mean', 'warped with all warp fns'})
plot(EFDARes.TimeMean,EFDARes.FuncMeanWarped,'Color',[0.6 .6 .6],'LineWidth',1);
plot(EFDARes.TimeMean,EFDARes.FuncMeanWarpedMean,'k','LineWidth',3);
M = EFDARes.FuncMeanWarpedMean;
SD = EFDARes.FuncMeanWarpedStd;
YPatch = [M + SD; flipud(M-SD)];
patch(sp(4),'XData',XPatch,'YData',YPatch,'EdgeColor','none','FaceColor','r','FaceAlpha',.4);
ylim([0 8]);


end





function [curFigRes,curFigResVar] = PlotMetricsVsSigma(T,Pars2Plot,SigmasColForPlot,figpos0)
curFigRes = Create_Reuse_Figure([],'Error vs sigma',figpos0 + [30 100 900 660]);
%curFigRes = Create_Reuse_Figure([],'Error vs sigma',[-950 100 900 660]); %Home

os.SpacingsOut = [40 40 80 80];
os.SpacingsBetween = [40 80];
os.YLabels = 'All';
os.XLabels = 'First'; % May change to "First" or to "FirstEvery2nd" or "BotLeftCorn"
os.XTickMarks = 'First';
os.YTickMarks = 'All';
[sp, osTop] = TiledAxesMy(curFigRes,5,4,[0.0 0 1 1],os);
delete(sp(end));
for iicol = 1:length(Pars2Plot) %columns of parameters
    
    icol = Pars2Plot(iicol);
    parname = T.Properties.VariableNames{icol};
    ind0 = find(strcmpi(T.Stage,'TimePad Speed Aver'));
    ind1 = find(strcmpi(T.Stage,'TimeNorm Speed Aver'));
    %ind2 = find(strcmpi(T.Stage,'PathNorm Speed Aver'));
    ind3 = find(strcmpi(T.Stage,'EFDA Speed Aver'));
    Sigmas = T.SigmaScale(ind1);
%     if iscell(SigmasColForPlot)
%         Sigmas = T.SigmaScale(ind1); %T.Sigmas(ind1,SigmasColForPlot{2});
%     else
%         Sigmas = T{ind1,SigmasColForPlot};
%     end
    subplot(sp(iicol));
    if ismember(iicol,[17:20]), xlabel('Noise SD','fontweight','bold'); end
    ylabel(parname,'fontweight','bold','Interpreter','none');
    plot(Sigmas,T{ind0,icol},'Color',[0.3 0.7 0.2],'DisplayName','No rescaling');
    plot(Sigmas,T{ind1,icol},'Color',[0.3 0.2 0.7],'DisplayName','Time Norm');
    %plot(Sigmas,T{ind2,icol},'g','DisplayName','PathNormAver');
    plot(Sigmas,T{ind3,icol},'Color','k','DisplayName','EFDA');
    yline(0,'Color','k','LineStyle','--','HandleVisibility','off');
    if iicol == 1, lg = legend('location','best'); lg.Position(1:2) = [0.78 0.15]; end
end

linkaxes(sp,'x');

curFigResVar = Create_Reuse_Figure([],'Variances vs sigma',figpos0 + [330 100 900 500]);
[spVar, osTopVar] = TiledAxesMy(curFigResVar,3,4,[0.0 0 1 1],os);

VarsAll(:,:,1) = T.VarsOrigAll;
VarsAll(:,:,2) = T.VarsAlignedAll;
VarsAll(:,:,3) = T.VarsTShiftAll;
VarsAll(:,:,4) = T.VarsAlignedMeanWarpedAll;
VarsAll = VarsAll(ind3,:,:);

% Plot for different methods for the variability in each plot of the first row
VarMethods = {'Area of Var(t)','Area of SD(t)','MSE of L2-distances to the mean','RMSE of L2-distances to the mean'};
VarTypes = {'Orig spat','Aligned spat','TimeShift temp','Temp2Spat'};
for iicol = 1:4
    plot(spVar(iicol),Sigmas,VarsAll(:,1,iicol),'m','DisplayName',VarMethods{1});
    plot(spVar(iicol),Sigmas,VarsAll(:,2,iicol),'Color',[1 .66 0],'DisplayName',VarMethods{2});
    plot(spVar(iicol),Sigmas,VarsAll(:,3,iicol),'g','DisplayName',VarMethods{3});
    plot(spVar(iicol),Sigmas,VarsAll(:,4,iicol),'b','DisplayName',VarMethods{4});
    ylabel(spVar(iicol),VarTypes{iicol},'fontweight','bold');
end
lg = legend(spVar(iicol),'Orientation','horizontal','Location','best');
lg.Position(1:2) = [0.08 0.94];

% Plot different types of variability in each plot of the second row
for iicol = 1:4
    delete(spVar(4+iicol));
    spVar(8+iicol).Position(2) = 0.22;
end
for iicol = 1:4
    plot(spVar(8+iicol),Sigmas,VarsAll(:,iicol,1),'r','DisplayName',VarTypes{1});
    plot(spVar(8+iicol),Sigmas,VarsAll(:,iicol,2),'y','DisplayName',VarTypes{2});
    plot(spVar(8+iicol),Sigmas,VarsAll(:,iicol,3),'c','DisplayName',VarTypes{3});
    plot(spVar(8+iicol),Sigmas,VarsAll(:,iicol,4),'k','DisplayName',VarTypes{4});
    ylabel(spVar(8+iicol),VarMethods{iicol},'fontweight','bold');
end
lg2 = legend(spVar(8+iicol),'Orientation','horizontal','Location','best');
lg2.Position(1:2) = [0.19 0.59];

end

















