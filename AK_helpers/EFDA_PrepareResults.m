function Result = EFDA_PrepareResults(durations, efdaStruct,nanlist,opts,nanlocs)
% Transfering most of the FDAWARP class fields to an output structure with more
% user-friendly names, plus computing a few extra variables, e.g. quantifying variabilities.
% Aleksei Krotov, Northeastern University <krotov.a@northeastern.edu>
% April 30, 2023

Result = struct;
Result.DurationsOrig = durations;
L = size(efdaStruct.f,1);

Result.TimeMean = efdaStruct.time; % this one is [0 1] as in the original package
Result.FuncOrig = efdaStruct.f;
Result.FuncAligned = efdaStruct.fn;
Result.FuncAlignedMean = mean(efdaStruct.fn,2,'omitnan');
Result.FuncAlignedStd = std(efdaStruct.fn,[],2,'omitnan');
%Result.
SRRFOrig = efdaStruct.q0; % Represent original SRRFs only if not aligning by derivatives or by integrals
%Result.
SRRFAligned = efdaStruct.qn; % Represent aligned SRRFs only if not aligning by derivatives or by integrals

Result.Warps = efdaStruct.gam;
%Result.
Psis = efdaStruct.psi;
% Plus compute "Time Speed" and "Time Shift"
Result.WarpsDdt = mydiff(efdaStruct.time,efdaStruct.gam); % derivative via central differences
Result.WarpsMinusUnity = efdaStruct.gam - linspace(efdaStruct.gam(1),efdaStruct.gam(end),L)';

% Prepare AlignedMean warped into original functions to reflect the contribution of
% temporal variability to the observed (raw) spatial variability.
for irec = 1:size(efdaStruct.f,2)
    gamI = invertGamma(efdaStruct.gam(:,irec));
    Result.FuncMeanWarped(:,irec) = warp_f_gamma(Result.FuncAlignedMean,gamI,Result.TimeMean,1);
end
Result.FuncMeanWarpedMean = mean(Result.FuncMeanWarped,2); %
Result.FuncMeanWarpedStd = std(Result.FuncMeanWarped,[],2); %

Result.FullDist = efdaStruct.FullDist; % distances to the mean in terms of original F. Only valid when 'orig' functions are aligned, not derivative or integral
Result.ElasticSpatDist = efdaStruct.ElasticAmpDist; % in terms of SRVF of F. Only valid when 'orig' functions are aligned, not derivative or integral
Result.ElasticPhDist = efdaStruct.ElasticPhDist; % in terms of SRVF of gammas. Only valid when 'orig' functions are aligned, not derivative or integral

% Elastic distances computed in the original method - I added it to the fdasrvf.timewarping,
% following the papers and the Srivastava's book.
% obj.ElasticAmpDist(ii) = sqrt(trapz(obj.time,(obj.mqn - obj.qn(:,ii)).^2));
% obj.ElasticPhDist(ii) = acos(trapz(obj.time,obj.psi(:,ii)));

% ADDING VARIABILITY METRICS - four different definitions to examine them.

%% An experimental feature: consider representing uniform scaling (time-norm) together with non-uniform scaling (time-warping)
% Time-warping results might become more interpretable, as all temporal scaling would be
% expressed in original time units. Furthermore, a landmark of interest in the aligned
% mean can be matched to the respective times of that landmark in the original ensemble.

% Careful, the variabilities may change too.

if opts.EFDAWarpsWithDurs
    Result.TimeMean = Result.TimeMean .* mean(Result.DurationsOrig,'omitnan');
    Result.Warps = Result.Warps .* Result.DurationsOrig';
    Result.WarpsMinusUnity = Result.Warps - Result.TimeMean;
    % Result.FuncMeanWarped change too?? No, just the timeMean changes. Up to the user
    % whether to use that time for plotting and/or computing variabilities, or to use
    % normalized time.
end


%% Another experimental feature, area-preservation
% In case researchers align a variable and assume that its integral should stay invariant
% (e.g. velocity and total distance travelled), this function will provide such
% invariance: intervals with shortened duration will have proportionally higher amplitudes.
if opts.PreserveIntegral 
    if strcmpi(opts.EFDAWarps,'orig')
        for irec = 1:size(efdaStruct.f,2)
            Result.FuncAligned(:,irec) = Result.FuncAligned(:,irec) .* Psis(:,irec).^2;
        end
        Result.FuncAlignedMean = mean(Result.FuncAligned,2,'omitnan');
        Result.FuncAlignedStd = std(Result.FuncAligned,[],2,'omitnan');
        for irec = 1:size(efdaStruct.f,2)
             gamI = invertGamma(efdaStruct.gam(:,irec));
             Result.FuncMeanWarped(:,irec) = warp_f_gamma(Result.FuncAlignedMean,gamI,Result.TimeMean,1);
        end
        
    end
end




%% LEAVE ONLY ONE VARIABILITY METRIC
%%
% Different ways to express variabilities.
[VarsOrig, VarNames] = Variabilities_EFDA(Result.FuncOrig, Result.TimeMean);
VarsAligned = Variabilities_EFDA(Result.FuncAligned, Result.TimeMean);
VarsTShift = Variabilities_EFDA(Result.Warps, Result.TimeMean); % Same ass WarpsMinusUnity
VarsAlignedMeanWarped = Variabilities_EFDA(Result.FuncMeanWarped, Result.TimeMean);

% Selected variability definition for ensemble of functions f [T x N] is 
% trapz(timemean,std(f)) ./ timemean(end);

% Expressing variability as the area of the SD(t) band (one-sided). Normalizing by average duration,
% in case opts.EFDAWarpsWithDurs is True.
% Assuming this will be the main, most intuitive, definition
Result.VarOrigSpat = VarsOrig(2);
Result.VarSpat = VarsAligned(2);
Result.VarTemp = VarsTShift(2);
Result.VarTemp2Spat = VarsAlignedMeanWarped(2);


% AGAIN ATTEMPT Kneip&RAMSAY (2008) or actually Ramsay Hooker Graves (2009)
N = size(efdaStruct.f,2);
x = Result.FuncOrig;
xhat = mean(Result.FuncOrig,2,'omitnan');
y = Result.FuncAligned;
yhat = Result.FuncAlignedMean;
h = Result.Warps;
Dh = Result.WarpsDdt;

MSETot = 1/N * sum(sum((x - xhat).^2,1)./L,2);
CRnum = mean(sum((Dh - mean(Dh,2,'omitnan')) .* (y.^2 - mean(y.^2,2,'omitnan')),1,'omitnan')./L,2,'omitnan');
CRdenom = mean(sum(y.^2,1,'omitnan')./L,2,'omitnan');
CR = 1 + CRnum / CRdenom;
MSEAmp = CR * mean(sum((y-yhat).^2,1,'omitnan')./L,2,'omitnan');
MSEPh = CR * sum(yhat.^2,1,'omitnan')./L - sum(xhat.^2,1,'omitnan')./L;

Result.MSETot = MSETot;
Result.MSEAmp = MSEAmp;
Result.MSEPh = MSEPh;

% % The next one is From Kneip & Ramsay (2008) Combining Registration and fitting....
% % determining independence of aligned functions and gammas.
% % It might remain unused.
% % Why are these two definitions produce different numeric results? Verify.....
% C1 = mean(trapz(Result.FuncOrig.^2)) ./ mean(trapz(Result.FuncAligned.^2));
% C2 = mean(trapz(Result.WarpsDdt .* Result.FuncAligned.^2)) ./ mean(trapz(Result.FuncAligned.^2));
% C = C1;
% Result.VarAmpMSEw = C * mean(trapz((Result.FuncAligned - Result.FuncAlignedMean).^2)./L);
% %C = mean(trapz(Result.WarpsDdt .* Result.FuncAligned.^2)) ./ mean(trapz(Result.FuncAligned.^2));
% %Result.VarAmpMSEw = C * trapz(Result.VarAmpF); Slightly different from the below
% MSTotal = mean(trapz((Result.FuncOrig - mean(Result.FuncOrig)).^2) ./ L);
%
%
% %%
% %%%%%% Careful here, verify expressions for variance and order of operations!
% % Check with Kneip&Ramsay; Check with Srivastava book too;
% %%
%
%
% % The next one is also from Kneip & Ramsay (2008). Weighted.
% % Supposedly, VarOrigAmpF = VarAmpMSEw + VarTempMSEw
% Result.VarTempMSEw = C * trapz(Result.FuncAlignedMean.^2) ./ L - trapz(mean(Result.FuncOrig,2,'omitnan')) ./ L;
%
% % From Srivastava: variance across the mean-aligned warped with all warps.
% Result.VarTempMSESri = trapz(Result.FuncMeanWarpedStd.^2) ./ L;%; ./ mean(durations);
%
% Result.VarTempTimeSpeed = trapz(var(Result.WarpsDdt,[],2)) ./ L; % ./ mean(durations);
% Result.VarTempTimeShift = trapz(var(Result.WarpsMinusUnity,[],2)) ./ L; % ./ mean(durations);



Result.NannedTrials = nanlist;
Result.IterCostFn = efdaStruct.qun;
Result.lambda = efdaStruct.lambda;
Result.method = efdaStruct.method;
Result.type = efdaStruct.type;
Result.TimeElapsed = efdaStruct.TimeEl;
Result.Iterations = efdaStruct.Iterations;
switch opts.EFDAPreserving
    case 'value'
        Result.AlignmentPreserving = 'signal values';
        Result.OrigFdaStruct = [];
    case 'integral'
        Result.AlignmentPreserving = 'integral values';
        Result.OrigFdaStruct = efdaStruct.origstruct;
    case 'derivative'
        Result.AlignmentPreserving = 'derivative values';
        Result.OrigFdaStruct = efdaStruct.origstruct;
end
    

% Process nan-trials
% If there were discarded time-series, add them to this structure as completely nan-series, 
% to retain the original trial numeration.
if ~isempty(nanlist)
    Flds = {'FuncOrig','FuncAligned','Warps','WarpsDdt','WarpsMinusUnity','FuncMeanWarped'}; % fields in Result structure that have Nsamp x Ntrials layout
    for ifld = 1:length(Flds)
        A = Result.(Flds{ifld});
        for inan = 1:length(nanlist)
            if nanlist(inan) > size(A,2)
                A(:,nanlist(inan)) = nan(size(A,1),1);
            else
                Temp = A(:,nanlist(inan):end);
                A(:,nanlist(inan)) = nan(size(A,1),1);
                A(:,nanlist(inan)+1:end+1) = Temp;
            end
        end
        Result.(Flds{ifld}) = A;
    end
    % A field with Ntrial x 1 layout
    FldsN1 = {'DurationsOrig','FullDist','ElasticSpatDist','ElasticPhDist'};
    for ifld = 1:length(FldsN1)
    
    A = Result.(FldsN1{ifld});
    for inan = 1:length(nanlist)
        if nanlist(inan) > length(A)
            A(end+1) = nan;
        else
        Temp = A(nanlist(inan):end);
        A(nanlist(inan)) = nan;
        A((nanlist(inan)+1):(end+1)) = Temp;
        end
    end
    Result.(FldsN1{ifld}) = A;
    end
end

% Add nans, that were previously interpolated.
for irec = 1:size(Result.FuncOrig,2)
    if ~isempty(nanlocs{irec})
        Result.FuncOrig(nanlocs{irec},irec) = nan;
        Result.FuncAligned(nanlocs{irec},irec) = nan;
        Result.Warps(nanlocs{irec},irec) = nan;
        Result.WarpsDdt(nanlocs{irec},irec) = nan;
        Result.WarpsMinusUnity(nanlocs{irec},irec) = nan;
        Result.FuncMeanWarped(nanlocs{irec},irec) = nan;
        
        %error('EFDA_PrepareResult::Replace nans if any');
    end
end
% y1(inannew1,:) = nan;


end





function [Vars, VarNames] = Variabilities_EFDA(F, t)
% units of F for example [m] and normalized time is [a.u.]
VarNames = {'Area of Var(t)','Area of SD(t)','MSE of L2-distances to the mean','RMSE of L2-distances to the mean'};

L = size(F,1); N = size(F,2);


Vars(1,1) = trapz(t, var(F,[],2)) / t(end); % [m^2]
Vars(1,2) = trapz(t,std(F,[],2)) / t(end); % [m]
Vars(1,3) = mean(trapz(t,(F - mean(F,2)).^2) / t(end)); % [m^2]   trapz(t, var(F,2)) / t(end)
Vars(1,4) = sqrt(mean(trapz(t,(F - mean(F,2)).^2) / t(end))); % [m]

end


function y = mydiff(t,x,varargin)
% Takes vector of time (may be not equidistantly spaced) Nx1 or 1xN,
% Takes vector of values or a matrix of values N x M
% Preferably, N is a column for both. Automatic reshaping works if M < N.

% Aleksei Krotov
% Northeastern University, 2023

warning('off','MATLAB:interp1:NaNstrip');
s = size(x);
ndof = s(2);
if any(strcmpi(varargin,'Method'))
    InterpMethod = varargin{find(strcmpi(varargin,'Method'))+1};
else
    InterpMethod = 'makima';  %Could be changed (e.g. to makima) to avoid undulations
end
flag_Speed = 0;
if any(strcmpi(varargin,'Speed'))
    flag_Speed = 1;
end

y = nan(size(x,1),size(x,2));
if flag_Speed
    y = nan(size(x,1),1);
end
if (nnz(isnan(x)) > 0.8 * numel (x))
    % Too many nans, returning a nan structure.
    %disp(sprintf('<strong>More than 80%% nans</strong>'));
    return
end
st = size(t);
if st(1) < st(2)
    t = t';
end
assert(iscolumn(t));
assert(length(t) == s(1),'Time and data lengths are not equal or data need to be transposed!');

% Remember nan locations
ind0 = isnan(x(:,1));


dt = t(3:end) - t(1:end-2); 
dx = x(3:end,:) - x(1:end-2,:);
% Computing double-forward differences
yy = dx ./ dt;
% They correspond to central differences at sample points 2:end-1
% We will substitute derivative values at these points with those at 2 and end-1
t1 = linspace(t(1),t(end),length(t)-2)';
for idof = 1:ndof
    y(:,idof) = interp1My(t1,yy(:,idof),t,InterpMethod); % ,'extrap' % Extrap sometimes shoots the result in infinity.
end


warning('on','MATLAB:interp1:NaNstrip');

y(ind0,:) = nan;

if flag_Speed
    y = sqrt(sum(y.^2,2));
end

end


function Vout = interp1My(X,V,Xq,method,extrapval)
% A wrapper of interp1 that returns nan in a few cases when the original function seems to
% throw errors. The checkups are not exhaustive.
% Aleksei Krotov
% Northeastern University, 2023

if all(isnan(V))
    Vout = nan(length(Xq),size(V,2));
    return
end

if all(isnan(Xq))
    Vout = nan(length(Xq),size(V,2));
    return
end

if all(isnan(X))
    Vout = nan(length(Xq),size(V,2));
    return
end

NNNans = nnz(~isnan(V(:,1)));
if NNNans < 3
    Vout = nan(length(Xq),size(V,2));
    return
end

if nargin > 4 && isstring(method) && isstring(extrapval)
    if strcmpi(method,'linear') && strcmpi(extrapval,'none')
        Vout = interp1(X,V,Xq);
        return
    elseif strcmpi(extrapval,'none')
        Vout = interp1(X,V,Xq,method);
        return
    else
        Vout = interp1(X,V,Xq,method,extrapval);
        return
    end
end

% Eventually pass to the regular interp1
if nargin < 4
    Vout = interp1(X,V,Xq);
elseif nargin < 5
    Vout = interp1(X,V,Xq,method);
else
    Vout = interp1(X,V,Xq,method,extrapval);
end

end

