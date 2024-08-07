function [opts, EFDAResult] = EFDA_alignmentPublishing(X,opts,varargin)

%!!! Check with forums or GPT how to best verify function inputs... Presence, correct
%type, correct range.

% !! Integral/Derivative. Don't use der/int just for finding warps.
% Instead use them to find warps and warped fns, then differentiate/integrate those fns,
% then calculate the rest on them.


% The data will be resampled and time-normalized (mandatory prior to EFDA) according to 
% the value of opts.FResampleTarget (default 20)

% In case the input data doesn't have timestamps, but, for resampling and time-normalizing
% (mandatory prior to EFDA) you want to specify a target number of samples, 
% use option NSamplesTarget.


% April 5 -- Quick partially implemented. Upsample after quick alignment.
% Compute and display variabilities.
% Add filtering,
% check if nans work
%  WarpsWithDurs - think carefully. Maybe in that case I should not show normalized time
%  at all, but instead the mean duration?...

% Quick - change the method. Use a small uniform subsample (50-100 of signals) to estimate
% the mean. Then just compute the warpings.

% Inputs - as a structure or as varargin?
% X - a cell array Mx1 or 1xM with entries-arrays of N(i) x 2, i=1:M;
% Could be a csv of 1st column time (same for every signal), and signals of various
% durations, the endings padded with nans.
% Or could be a matrix where time is already resampled (equal durations assumed).

% Take inputs as matrix or cell or struct, with time and some X.
% Actually, if there's a button READ DATA, then the data should be able to read e.g. csv
%
% Preview quick stats and preview rough alignment
%     Durations histogram and stats
%     FFT power?
%     Lowpass filter via something quick and robust
%     Downsample to Nyquist 20 Hz, i.e. 40 Hz
%     If more than 200 signals, pick random 200
%     Nan-analyze
%     EFDA-align
%     Interpolate warps
%     Take original functions, warp them, plot the result and the warps
%     Preview. Probably aim for specific execution time (i.e. if better accuracy required, boost sampling and decrease lambda, but pick fewer functions)
%     Offer different nan processing: interpolate short ones or completely remove
%
% Do the warp alignment
%     Take Lambda, N samples,
%     Take optional filtering via Box xN, Savitsky, Or Butter.
%     Analyze Nans
%     Offer parallerl
%     EFDA-align
%     Plot. Compute variabilities. Display time-shift. Provide hint what time-shift means
%     Optionally compute warps and time-shifts scaled with original times
%     Optionally do EFDA-align w.r.t. derivatives or wrt integrals
%
% Return all the values as a structure.
%
% Do plotting in a single figure.
%
% Make a separete code where the same figure has callbacks, i.e. buttons, Load dialogue, sliders etc.


%% Parse inputs
% check types, check that within range
if nargin < 2
    opts = [];
end

%%% DEBUG ONLY - most frequently changed parameters
% opts.FSResampleTarget = 40;
% opts.EFDAParallel = opts.EFDAParallel;
% opts.EFDAWarpsWithDurs = 0;
% opts.EFDALambda = 0.01;
% opts.NSampPlot = 120;
% opts.NLinesPlotMax = 100;
% opts.EstimateExecutionTime = 0;



opts = EFDA_ParseOpts(opts,varargin{:});


% Return default vars values
if nargout == 1
    return
end

figpos0 = [1920 0 0 0]; % several monitors
figpos = [100 200 1280 720] + figpos0;
figname = 'Time-warping alignment';

% Check that input is cell!
% Check that each cell is an Nx2 column

%% Determine sampling rate (if FFT), determine Nresamp, check for nans, time-normalize
if iscell(X) % cell where y{i} = two-column vector of time and value
    [ynorm,tnorm,dur,opts,FSorig,Ntrials,~,nanTriallist,nanlocs] = EFDA_TimeNormSeries(X,opts);
else
    error('Input is not cell');
end
opts.nanTriallist = nanTriallist;
opts.Ntrials = Ntrials;
opts.DurMean = mean(dur,'omitnan');
trialnotnans = ones(size(ynorm,2),1) == 1;
trialnotnans(nanTriallist) = 0;

% Estimate time-warping execution time
if opts.EstimateExecutionTime || opts.EFDAQuick
    EFDA_tEst_1thread = EstimateEFDATime(ynorm,opts);
    opts.EFDA_tEst_1thread = EFDA_tEst_1thread;
    if EFDA_tEst_1thread < 1
        fprintf('<strong>Single-core execution is already faster than 1 s, removing the QuickEFDA flag</strong>n');
        opts.EFDAQuick = 0;
    end
end


if opts.EFDAGraphics
    EFDAFig = EFDA_PlotInput(ynorm,tnorm./tnorm(end),dur,opts,[],figname,figpos);
end
%% Consider filtering before EFDA

% function filter data EFDA

if ~isempty(opts.Filter)
    switch opts.Filter{1}
        case 'Butter'
            % In case time-normalization was made with too low frequency, redo it before filtering
            if opts.FSResampleTarget < opts.Filter{3}
                warning('Resampling frequency after normalization was increased to %.0f to match the double of the filter cutoff frequency',2*opts.Filter{3});
                opts.FSResampleTarget = opts.Filter{3};
                [ynorm,tnorm,dur,opts,FSorig,~,~,nanTriallist,nanlocs] = EFDA_TimeNormSeries(X,opts,'Quiet');
            end
            FilterParams.Type = 'ButterLow';
            FilterParams.Order = opts.Filter{2};
            FilterParams.FCutoff = opts.Filter{3};
            ynormF = EFDA_LowFilt(ynorm,FSorig,FilterParams,opts);
            % downsample again?
            %[ynorm,tnorm,dur,FSorig,[],[],nanTriallist] = TimeNormSeries(X,vars);
        case 'SavGol'
            FilterParams.Type = 'SavGol';
            FilterParams.Order = opts.Filter{2};
            FilterParams.FrameLength = opts.Filter{3};
            ynormF = EFDA_LowFilt(ynorm,FSorig,FilterParams);
        otherwise
            warning('Filter input not formatted properly. Skipping the filtering');
    end
    ynorm = ynormF;
end




%%
% Display GUI results?


%% Run EFDA
% if not gui.
tic
tEFDA0 = toc;
efdaResult = EFDA_align(ynorm(:,trialnotnans),tnorm./tnorm(end),opts,nanlocs);
efdaStruct = EFDA_fdawarp2struct(efdaResult);
EFDAResult = EFDA_PrepareResults(dur(trialnotnans),efdaStruct,nanTriallist,opts,nanlocs);
tEFDA1 = toc;
disp(sprintf('Alignment time %.1f s',tEFDA1-tEFDA0));
opts.LastEFDAExecutionTime = tEFDA1 - tEFDA0;


%% Plot the results if after EFDA
if opts.EFDAGraphics
    EFDAFig = EFDA_PlotResult(EFDAResult,opts,EFDAFig);
end


% DEBUG Display variabilities after Ramsay
% disp(sprintf('%.3f + %.3f - %.3f = %.3f;   R^2 = %.3f',EFDAResult.MSEAmp,EFDAResult.MSEPh,EFDAResult.MSETot,(EFDAResult.MSEAmp + EFDAResult.MSEPh - EFDAResult.MSETot),EFDAResult.MSEPh./EFDAResult.MSETot));


end





















%% Modules





%% EFDA settings


%% Convert fdawarp object to a regular matlab structure



%% Prepare a more convenient results structure




%% Plot the results after EFDA











%% Auxilliaries

function opts = EFDA_ParseOpts(optsIn,varargin)
% Have a table - ?cell or structure array? - of fieldnames and default values
% Read those that were supplied, update values
% Assign the updated values.

% Flags, such as EFDAQuick, Quiet, Graphics, PreserveIntegral, can be passed as varargin


% Default vars values
optsDef.EFDAQuick = 0;
optsDef.EFDAGraphics = 1;
optsDef.Gui = 0;
optsDef.PreserveIntegral = 0; % Experimental. After alignment, rescale amplitude of the signals to keep area under it invariant. Computed in EFDA_PrepareResults, may require heavier filtering or using larger Lambda.
optsDef.TimeStampsSupplied = 1;
optsDef.FSResampleTarget = 20; % [Hz] if time is in seconds; [1/timeunit] if otherwise (e.g. milliseconds or days or arclength-meters), about 1/2 of the resampling frequency in the longest signal, consider Nyquist for estimation.
optsDef.NSamplesTarget = nan;
optsDef.RemoveLongOutliers = 1; % Remove trials which are longer/shorter than Mean+3SD of durations.
optsDef.FSFFTUse = 0;
optsDef.FSFFTCutoffPower = 0.95;
optsDef.Filter = []; % {'Butter',ORDER [4],CUTOFF FREQ HZ [10]} OR {'SavGol',ORDER [3],FRAME LENGTH [5]}
optsDef.Quiet = 0;

optsDef.EFDAMethod = 'DP'; % This script was tested only for DP, but the others from the fdasrsf library might work as well.
optsDef.EFDALambda = 0.01; % Warp roughness penalty. Set 0 to match the math most rigorously. Set 0.001-0.05 for smoother warps.
optsDef.EFDAParallel = 1;
optsDef.EFDAPreserving = 'values'; % Experimental. Can 'derivative', or 'integral', of the original signals to obtain warps. The rest of computation is still done on the original signal.
optsDef.EFDAMaxIter = 2; % In my experience, the result does not change after 1-2 iterations.
optsDef.EFDAWarpsWithDurs = 0; % Resample normalized time to match the mean duration.
optsDef.MaxNanInterp = 0.1; % Completely ignore (and turn into nans) the time-series which have more than this fraction of nans. Otherwise, these nans will be interpolated (interp1, 'makima').
optsDef.NSampPlot = 150; % TBD
optsDef.NLinesPlotMax = 90; % 0 is all;  too large NSamp and NLines would prevent from exporting plot as vector.
optsDef.EstimateExecutionTime = 1; % Pick 5% of the data at desired sampling rate, try aligning, extrapolate the duration to the whole sample set.
% I found that in my i7-10700K, Time ~ 1.9e-4 * size(Func,2) * size(Func,1)^2;
% I.e. after TimeNorm, pick such number of functions >=2 that would take >= 1 second,
% if even 2 functions have too many timesamples, downsample them. Do optimum_reparam, get
% the time, and predict the total evaluation time.
% OR. to account for slower machines, pick any 2 functions, downsample or upsample to take
% (for my i7-10700K) 0.1 s, and run optimum_reparam repeatedly until total time reaches 1
% s. Then predict.


opts = optsDef;
if ~isempty(optsIn)
    assert(isstruct(optsIn),'The options variable needs to be a structure; Run the function with no inputs and a single output to see the default options.');
    flds = fieldnames(optsIn);
    fldsvalid = fieldnames(optsDef);
    for ifld = 1:length(flds)
        if any(ismember(flds{ifld},fldsvalid))
            opts.(flds{ifld}) = optsIn.(flds{ifld});
        else
            disp(sprintf('EFDA_Alignment Input: Ignoring the option field "%s". Re-run the function with no inputs and a single output to see valid field names.',flds{ifld}));
        end
    end

end
if ~isempty(varargin) && ~isempty(varargin{1,1})
    for iarg = 1:length(varargin)
        argg = varargin{iarg};
        fldsvalid = fieldnames(optsDef);
        if any(strcmpi(argg,fldsvalid))
            opts.(fldsvalid{(strcmpi(argg,fldsvalid))}) = 1;
        end
        if strcmpi(argg,'derivative')
            opts.EFDAPreserving = 'derivative';
        end
        if strcmpi(argg,'integral')
            opts.EFDAPreserving = 'integral';
        end

    end
end

end




function EFDA_tEst_1thread = EstimateEFDATime(ynorm,opts)
% Pick 2 arbitrary signals from the timenormed series
inds = randi(size(ynorm,2),1,2);
% If they have more than 200 time-samples, it might already challenge some older
% computers; if so, downsample.
Nmax = 200;
ytest0 = ynorm(:,inds);
N = size(ytest0,1);
ytest = ytest0;
if size(ytest0,1) > Nmax
    ytest = interp1(linspace(0,1,N)',ytest0,linspace(0,1,Nmax)');
    N = Nmax;
end
TimeNorm = linspace(0,1,N)';
q1 = f_to_srvf(ytest(:,1),TimeNorm);
q2 = f_to_srvf(ytest(:,2),TimeNorm);
tic
t1 = toc;
t2 = toc;
k = 0;
while (t2-t1) < 0.5 || k < 5
    gam_o = optimum_reparam(q1,q2,TimeNorm,opts.EFDALambda,opts.EFDAMethod);
    t2 = toc;
    k = k + 1;
end
coeffA = (t2 - t1) ./ k ./ N.^2;
EFDA_tEst_1thread = coeffA * size(ynorm,2) * size(ynorm,1).^2 * (2+opts.EFDAMaxIter);
fprintf('Estimated time %.1f s for time-warp alignment on one thread. Expect ~4-6 times faster for Parallel on 8 workers. \n',EFDA_tEst_1thread);
end


