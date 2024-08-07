function [ynorm,tnorm,dur,opts,FSorig,Ntrials,NTrialsUsed,nanTriallist,nanlocs] = EFDA_TimeNormSeries(src,opts,varargin)

% Trials with more nans than MaxNanTol, or with trials at the edges,
% will not be used for time-norm averaging or for the efda-alignment.
% Other trials with nans will be interpolated and used for the whole process - and only
% after efda the nans will be replaced at approximate locations of the original data's
% nans (assuming little-to-no warping).


FSTarget = opts.FSResampleTarget;
FFT_flag = opts.FSFFTUse;
FFT_CutoffPower = opts.FSFFTCutoffPower;
MaxNan = opts.MaxNanInterp;

flag_messages = 1;
if any(strcmpi(varargin,'Quiet')) || opts.Quiet == 1
    flag_messages = 0;
end


dur = nan(length(src),1);
L = nan(length(src),1);
Ntrials = length(src);
flag_Discard = zeros(Ntrials,1);

% If no time is provided (src is a cell where src{i} = [N x 1]), use No. of samples as
% time.
cellsz = nan(length(src),1);
for iitrial = 1:length(src)
    y0 = src{iitrial};
    if size(y0,2) > size(y0,1)
        y0 = y0';
    end
    cellsz(iitrial) = size(y0,2);
end
if ~all(cellsz == cellsz(1))
    error('data dimensions in the cell input must be consistent');
end
if cellsz(1) == 1
    opts.TimeStampsSupplied = 0;
    for iitrial = 1:length(src)
        y0 = src{iitrial};
        Nsamp = length(y0);
        y0(:,2) = y0;
        y0(:,1) = (1:Nsamp)';
        src{iitrial} = y0;
    end
end

% detect the duration of the longest signal and number of samples
for iitrial = 1:length(src)
    y0 = src{iitrial};
    if size(y0,1) < size(y0,2)
        y0 = y0'; %Vertical vector/matrix
        src{iitrial} = y0;
    end
    torig = y0(:,1);
    dur(iitrial,1) = torig(end,1);
    L(iitrial,1) = length(torig);

    if isnan(torig(end))
        flag_Discard(iitrial,1) = 1;
    end
end
FSorig = mean(L./dur,'omitnan');

% Remove too long trials? use 3*SD
indlong = find(abs(dur - mean(dur,'omitnan')) > 3 * std(dur,[],'omitnan')); % change to median and mad if the durations are largely non-normal
opts.NLongOutlier = 0;
if ~isempty(indlong)
    if opts.RemoveLongOutliers
        opts.NLongOutlier = length(indlong);
        flag_Discard(indlong) = 1;
        dur(indlong) = nan;
        if flag_messages
            disp(sprintf('%d trials were discarded due to their outlier-like duration',length(indlong)));
        end
    else
        if flag_messages
            warning('There are %d trials whose duration appears as an outlier (>3 SD) and might cause oversampling for the rest of the data. Consider discarding these trials (vars.RemoveLongOutliers = TRUE)',length(indlong));
        end
    end
end


%%% Optionally run FFT with power threshold detection
if FFT_flag
    disp(FFT_CutoffPower);
    error('No FFT set up yet');
end
if any(strcmpi(FSTarget,{'orig','original','source'}))
    FSTarget = FSorig / 2;
end


% if flag_messages
%     warning('Table with some shorter trial-columns not supported yet')
% end

if opts.TimeStampsSupplied == 1
    if FSorig < FSTarget
        error('Upsampling not supported, average sampling rate appears to be %.0f but requested sampling rate is %.0f',FSorig,FSTarget);
    end
    Nresamp = max(dur,[],'omitnan') .* FSTarget * 2; % Assuming FSTarget includes some Nyquist-related reserve, i.e. is larger than the expected cutoff frequency for low-pass processes.
    Nresamp = ceil(Nresamp / 10) * 10; % just make Nresamp a little nicer, with a bit of extra reserve.
elseif ~isnan(opts.NSamplesTarget)
    Nresamp = ceil(opts.NSamplesTarget);
else
    if flag_messages
        warning('No timestamps provided and no opts.NSamplesTarget specified. Data will be resampled to the length of the shortest signal - which might still be too much and cause excessively long alignment.')
    end

    Nresamp = max([1 min(L,[],'omitnan')]);
end


ynorm = nan(Nresamp,length(src));
nanlocs = cell(length(src),1);
warning('off','MATLAB:interp1:NaNstrip');
for iitrial = 1:length(src)

    y0 = src{iitrial};
    torig = y0(:,1);
    yorig = y0(:,2);

    % Check for large gaps
    nanFramelist = isnan(yorig(:,1));
    Ngaps = nnz(nanFramelist);
    if Ngaps > MaxNan * L(iitrial,1) | isempty(yorig)
        flag_Discard(iitrial,1) = 1;
        continue
    end

    % Check for gaps at the ends
    if nanFramelist(1) || nanFramelist(end)
        flag_Discard(iitrial,1) = 1;
    end

    % skip completely only if no data
    if all(nanFramelist)
        continue
    end

    % interpolate
    irangeorig = linspace(0,1,L(iitrial,1))';
    irangenew = linspace(0,1,Nresamp)';
    y1 = interp1(irangeorig,yorig(:,1),irangenew,'makima',NaN);

    % whatever nans were interpolated, re-add them. Even setting 'extrap' as nan does not
    % prevent interp1 from extrapolating unrealistic values.

    % scale nan locations
    inanorig = find(nanFramelist);
    inannew0 = unique(round(inanorig ./ L(iitrial,1) * Nresamp));

    % propagate by +-1 around that range(s)
    inannew1 = sort(unique([inannew0 - 1 inannew0 inannew0 + 1]));
    inannew1(inannew1>Nresamp) = [];
    inannew1(inannew1<1) = [];

    % replace these nans if they are at the edges, to make sure no-efda time-norm
    % averaging displays correctly.
    if nanFramelist(1) || nanFramelist(end)
        y1(inannew1,:) = nan;
    end

    % save the locations, to replace the interpolated nans after EFDA runs.
    nanlocs{iitrial,1} = inannew1;
    ynorm(:,iitrial)  = y1;
    

end
% With actual nan-trials from the Discard list, interpolate 


warning('on','MATLAB:interp1:NaNstrip');
NTrialsUsed = nnz(~flag_Discard);
nanTriallist = find(flag_Discard);

if opts.PreserveIntegral && ~strcmpi(opts.EFDAWarps,'orig')
warning('Area preserving warping is only supported when aligning original signals, not derivatives or integrals. Ignoring area preserving and doing a default shape-preserving warping.');
opts.PreserveIntegral = 0;
end

if opts.PreserveIntegral
    DurNormTemp = dur ./ mean(dur,'omitnan');
    ynorm = ynorm .*  DurNormTemp';
end

tnorm = linspace(0,1,size(ynorm,1))';
if opts.EFDAWarpsWithDurs
    tnorm = tnorm .* mean(dur,'omitnan');
end

end