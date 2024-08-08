function [opts, EFDAResult] = EFDA_alignmentPublishing(X,varargin)

% Performs time-warping alignment of multiple time series using the time-warping alignment (EFDA) algorithm.
% It offers various options for preprocessing, alignment, and post-processing, 
% including filtering, resampling, outlier removal, and visualization.

% Inputs
% X: A cell array of time series data. Each cell contains a Tx2 or Tx1 matrix, where 
% T is the number of time points and the second column (if present) represents the time stamps.
% varargin: Optional input arguments for specifying function parameters. Can be 
% passed as a structure or as individual arguments.

% Outputs
% EFDAResult: A structure containing the alignment results, including warped time series, 
% warping functions, and other relevant information.
% opts: the options structure detailing customizable options. (if no opts or keywords
% supplied, use that output to see the default settings).

% opts: A structure containing the function parameters and their values.
% Pass binary-flag options as keywords or pass any setting in second-argument structure opts.
    % EFDAQuick: If set to 1, performs alignment on a subset of the data, while also
    % downsampling it, to estimate the emerging mean; Then the remaining warping functions
    % are found. Expected time on a single thread 2 seconds. Use for large datasets where
    % M * N^2 > 10^6, for M signals of N samples.
    % EFDAGraphics: If set to 1, visualizes results. Default: 1.
    % PreserveIntegral: If set to 1, rescales aligned signals to preserve the integral. Experimental. Default: 0.
    % FSResampleTarget: Target sampling frequency for resampling. Default: 20 Hz.
    % NSamplesTarget: Target number of samples for resampling (overrides FSResampleTarget). Default: NaN.
    % RemoveLongOutliers: If set to 1, removes trials with duration significantly different from the mean. Default: 1.
    % Filter: Specifies the filter type (e.g., 'Butter', 'SavGol') and parameters. Default: empty.
    % Quiet: If set to 1, suppresses output messages and graphics. Default: 0.
    % EFDAMethod: EFDA method to use (e.g., 'DP'). Default: 'DP'.
    % EFDALambda: Warp roughness penalty for EFDA. Default: 0.01.
    % EFDAParallel: If set to 1, uses parallel computation. Default: 1.
    % EFDAPreserving: Specifies whether to preserve values, derivatives, or integrals during warping. Experimental. Default: 'values'.
    % EFDAMaxIter: Maximum number of EFDA iterations. Default: 2.
    % EFDAWarpsWithDurs: If set to 1, resamples normalized time to match the mean duration. Default: 0.
    % MaxNanInterp: Maximum fraction of NaN values allowed for interpolation. Default: 0.1 (10%) per signal
    % NSampPlot: Number of samples to plot. Default: 150.
    % NLinesPlotMax: Maximum number of lines to plot. Default: 90.
    % EstimateExecutionTime: If set to 1, estimates the execution time of EFDA. Default: 1.


% Aleksei Krotov
% Northeastern University
% 2024


%% Parse inputs
% check types, check that within range

[X, opts] = ParseInput(X,varargin{:});



%%% DEBUG ONLY - most frequently changed parameters
% opts.FSResampleTarget = 40;
% opts.EFDAParallel = opts.EFDAParallel;
% opts.EFDAWarpsWithDurs = 0;
% opts.EFDALambda = 0.01;
% opts.NSampPlot = 120;
% opts.NLinesPlotMax = 100;
% opts.EstimateExecutionTime = 0;


% 
% opts = EFDA_ParseOpts(opts,varargin{:});


% Return default vars values
if nargout == 1
    return
end

figpos0 = [0 0 0 0]; % change if plotting not on the main monitor 
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


%% Run EFDA
tt1 = tic;
tEFDA0 = toc(tt1);
efdaResult = EFDA_align(ynorm(:,trialnotnans),tnorm./tnorm(end),opts,nanlocs);
efdaStruct = EFDA_fdawarp2struct(efdaResult);
EFDAResult = EFDA_PrepareResults(dur(trialnotnans),efdaStruct,nanTriallist,opts,nanlocs);
tEFDA1 = toc(tt1);
disp(sprintf('Alignment time %.1f s',tEFDA1-tEFDA0));
opts.LastEFDAExecutionTime = tEFDA1 - tEFDA0;


%% Plot the results after EFDA if graphics enabled
if opts.EFDAGraphics
    EFDAFig = EFDA_PlotResult(EFDAResult,opts,EFDAFig);
end



end



















%% Auxilliaries

function [X, opts] = ParseInput(X,varargin)

N = nargin;

if N >= 2
    if isstruct(varargin{1})
        opts = varargin{1};
        Vins = varargin(2:end);
        opts = EFDA_ParseOpts(opts,Vins{:});
    else
        opts = [];
        opts = EFDA_ParseOpts(opts,varargin{:});
    end
else
    opts = [];
    opts = EFDA_ParseOpts(opts);
end

InputMsg = 'The first argument must be a Mx1 cell of Tx2 or Tx1 numbers or a TxM matrix';
Cond = (iscell(X) || ismatrix(X)) && (length(X) > 1);
assert(Cond,InputMsg);
if ismatrix(X)
    % Convert to a cell
    M = size(X,2);
    Xc = cell(M,1);
    for j = 1:M
        Xc{j} = X(:,j);
    end
end

end



function opts = EFDA_ParseOpts(optsIn,varargin)
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
        %fldsvalid = fieldnames(optsDef);
        fldsvalid = {'EFDAQuick','EFDAGraphics','PreserveIntegral','TimeStampsSupplied','RemoveLongOutliers','FSFFTUse','Quiet'};
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

% For several ensembles I have worked with on three different Windows PCs,
% single-core computation time has been fairly well predicted as ~ M * N^2,
% where M is number of signals to align, and N is their number of samples.

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


