function objwarp = EFDA_align(FuncOrig,TimeOrig,opts,nanlocs)

option.parallel = opts.EFDAParallel; % turns on MATLAB parallel processing - speeds up a lot!
option.closepool = 0; % determines whether to close matlabpool
option.smooth = 0; % smooth data using standard box filter
option.sparam = 25; % number of times to run filter if enabled
option.method = opts.EFDAMethod; % optimization method (DP, DP2, SIMUL, RBFGS).
% DP reacts well to different lambdas, others seem to be insensitive to it.
option.w = 0.0; % BFGS weight, only for RBFGS
option.spl = true; % use spline interpolation (doesn't matter for our signals)
option.MaxItr = opts.EFDAMaxIter;  % maximum iterations (if result is suboptimal, increase).
% In my experience, little to no improvement results from more than 2 iterations.
option.IterationHistory_flag = 0; % turn on to examine all iterations (< MaxItr); increases computation time by up to 2 times
option.Quiet = opts.Quiet;


if opts.EFDAQuick % Aim for 1-2-second performance on a single-thread
    Func2Align = FuncOrig;
    if strcmpi(opts.EFDAPreserving,'derivative')
        Func2Align = mydiff(TimeOrig,FuncOrig,'NoReshape');
    end
    if strcmpi(opts.EFDAPreserving,'integral')
        Func2Align = cumtrapz(TimeOrig,FuncOrig);
    end

    
    % Now I have the full-blown version working. For the quick one, perform alignment on
    % whatever Func2Align - only to find the mean - then compute remaining warps to that mean.

    % Then, just copy the approach from the full version in case of derivative and
    % integral


    % Computation burden
    tFull = opts.EFDA_tEst_1thread; % = t_single_optimum_reparam .* (Niter + 2) .* N^2 * M
    N = size(FuncOrig,1);
    M = size(FuncOrig,2);

    % is it sufficient to just decrease NSamples to 20? 20 provides relatively smooth
    % curve in a regular 
    %t20 = tFull ./ N^2 .* 20.^2;

    % So to decrease time to 1 seconds, take between 3 to M functions, resampled to 20
    % to NFramesTarget frames.
    MinS = [20 3];
    %mintime = tFull ./ N^2 ./ M .* MinS(1).^2 .* MinS(2);
    % Choose MinS(1) < n < N of samples per function and select a subset of MinS(2) < m < M  functions
    % To estimate EFDA mean from them. Then find optimum warps for the rest (M-m)
    % downsampled functions. Then upsample everything back to N.

    
    % How to choose n and m?
    % While keeping one-thread time under 1 s, maximize such a combination of n and m that each of
    % them is as high ratio of [N M] as possible. Empirically, the expression of "scale"
    % below appears appropriate.
    %time = nan(N,1);
    scale = nan(N,1);
    m = nan(N,1);
    %n = nan(N,1);
    for i = MinS(1):N
        %n(i) = i;
        m(i) = min([M floor((N^2 * M / tFull) ./ i.^2)]);
        %time(i) = tFull ./ N^2 ./ M .* i.^2 .* m(i);
        scale(i) = log(sum([N M] ./ [i.^0.75 m(i)]));
    end
    [~,n] = min(scale);
    n = max([10 n]);
    m = max([min([M floor((N^2 * M / tFull) ./ n.^2)]) 3]);

    NM = [length(TimeOrig) size(FuncOrig,2)];
    Msubquick = m;%ceil(NM(2) .* 0.2);
    indsubquick = sort(unique(randperm(NM(2),Msubquick))); % A subsample of original functions

    NResampQuick = n;
    FuncQuick2Align = nan(NResampQuick,M);
    TimeQuick = linspace(0,1,NResampQuick)';
    for itrial = 1:M
        FuncQuick2Align(:,itrial) = interp1My(TimeOrig,Func2Align(:,itrial),TimeQuick);
    end

    tt2 = tic;
    t1 = toc(tt2);
    
    % Find the aligned mean using the reduced functions
    objwarp0 = fdawarp(FuncQuick2Align,TimeQuick); % Creating object
    objwarp = fdawarp(FuncQuick2Align(:,indsubquick),TimeQuick); % Creating object
    option.parallel = 0;
    option.Quiet = 1;
    option.MaxItr = 1;
    objwarp = objwarp.time_warping(opts.EFDALambda,option); % Alignment

    % Then find alignment to that mean from all the original functions, same reduced number of samples
    objwarp0.fmean = objwarp.fmean;
    objwarp0.mqn = objwarp.mqn;
    objwarp0.stats = objwarp.stats;
    objwarp0.qun = objwarp.qun;
    objwarp0.lambda = objwarp.lambda;
    objwarp0.method = objwarp.method;
    objwarp0.gamI = objwarp.gamI;
    objwarp0.rsamps = objwarp.rsamps;
    objwarp0.type = objwarp.type;
    objwarp0.FullDist = [];
    objwarp0.ElasticAmpDist = [];
    objwarp0.ElasticPhDist = [];
    isub = 1;
    for itrial = 1:M
        if ismember(itrial,indsubquick)
            objwarp0.q0(:,itrial) = objwarp.q0(:,isub);
            objwarp0.gam(:,itrial) = objwarp.gam(:,isub);
            objwarp0.fn(:,itrial) = objwarp.fn(:,isub);
            objwarp0.qn(:,itrial) = objwarp.qn(:,isub);
            objwarp0.psi(:,itrial) = objwarp.psi(:,isub);
            isub = isub + 1;
        else
            objwarp0.q0(:,itrial) = f_to_srvf(FuncQuick2Align(:,itrial),TimeQuick);
            objwarp0.gam(:,itrial) = optimum_reparam(objwarp0.mqn,objwarp0.q0(:,itrial),TimeQuick,opts.EFDALambda,opts.EFDAMethod);
            objwarp0.fn(:,itrial) = warp_f_gamma(objwarp0.f(:,itrial),objwarp0.gam(:,itrial),TimeQuick);
            objwarp0.qn(:,itrial) = f_to_srvf(objwarp0.fn(:,itrial),TimeQuick);
            objwarp0.psi(:,itrial) = f_to_srvf(objwarp0.gam(:,itrial),TimeQuick);
        end
    end
    
    % Then interpolate to original samples and we're done here.
    flds = {'f','fn','qn','q0','fmean','mqn','gam','psi','gamI'};
    for ifld = 1:length(flds)
        a = interp1(TimeQuick,objwarp0.(flds{ifld}),TimeOrig,'makima',nan);
        objwarp0.(flds{ifld}) = a;
    end
    objwarp0.time = TimeOrig;

    % Find the distances and stats
    isub = 1;
    for itrial = 1:M
        if ismember(itrial,indsubquick)
            objwarp0.FullDist(itrial) = objwarp.FullDist(isub);
            objwarp0.ElasticAmpDist(itrial) = objwarp.ElasticAmpDist(isub);
            objwarp0.ElasticPhDist(itrial) = objwarp.ElasticPhDist(isub);
            isub = isub + 1;
        else
            objwarp0.FullDist(itrial) = sqrt(trapz(objwarp0.time,(objwarp0.fn(:,itrial) - objwarp0.fmean).^2));
            objwarp0.ElasticAmpDist(itrial) = sqrt(trapz(objwarp0.time,(objwarp0.mqn - objwarp0.qn(:,itrial)).^2));
            objwarp0.ElasticPhDist(itrial) = acos(trapz(objwarp0.time,objwarp0.psi(:,itrial)));
        end
    end
    objwarp = objwarp0;
    std_f0 = std(objwarp.f, 0, 2); 
    std_fn = std(objwarp.fn, 0, 2); 
    fgam = zeros(length(TimeOrig),size(FuncOrig,2));
    for ii = 1:size(FuncOrig,2)
        fgam(:,ii) = warp_f_gamma(objwarp.fmean,objwarp.gam(:,ii)',objwarp.time);
    end
    std_fgam = std(fgam, 0, 2);
    objwarp.stats.orig_var = trapz(objwarp.time,std_f0.^2);
    objwarp.stats.amp_var = trapz(objwarp.time,std_fn.^2); 
    objwarp.stats.phase_var = trapz(objwarp.time,std_fgam.^2); 

    obj0 = fdawarp2struct(objwarp);
    objwarp = obj0;


    % Re-convert derivative/integral if needed
    if strcmpi(opts.EFDAPreserving,'derivative')
        % Integrate fn
        objwarp.f = FuncOrig;
        objwarp.fn = cumtrapz(TimeOrig,obj0.fn);
        objwarp.fmean = cumtrapz(TimeOrig,obj0.fmean);
        %objwarp.alignment = 'from the derivative alignment';
    elseif strcmpi(opts.EFDAPreserving,'integral')
        objwarp.f = FuncOrig;
        objwarp.fn = mydiff(TimeOrig,obj0.fn,'NoReshape');
        objwarp.fmean = mydiff(TimeOrig,obj0.fmean);
        %objwarp.alignment = 'from the integral alignment';
    end
    if any(strcmpi(opts.EFDAPreserving,{'derivative','integral'}))
        % recompute q0, qn, mqn, gam, psi, stats, Dist, gamI
        % zero qun
        for itrial = 1:size(FuncOrig,2)       
            objwarp.q0(:,itrial) = f_to_srvf(objwarp.f(:,itrial),TimeOrig);
            objwarp.qn(:,itrial) = f_to_srvf(objwarp.fn(:,itrial),TimeOrig);
            objwarp.gam(:,itrial) = optimum_reparam(objwarp.qn(:,itrial),objwarp.q0(:,itrial),TimeOrig,opts.EFDALambda,opts.EFDAMethod);
            objwarp.psi(:,itrial) = f_to_srvf(objwarp.gam(:,itrial),TimeOrig);
        end
        objwarp.mqn = mean(objwarp.qn,2,'omitnan');
        objwarp.gamI = SqrtMeanInverse(objwarp.gam');
        for itrial = 1:size(FuncOrig,2)
            objwarp.FullDist(itrial) = sqrt(trapz(objwarp.time,(objwarp.fn(:,itrial) - objwarp.fmean).^2));
            objwarp.ElasticAmpDist(itrial) = sqrt(trapz(objwarp.time,(objwarp.mqn - objwarp.qn(:,itrial)).^2));
            objwarp.ElasticPhDist(itrial) = acos(trapz(objwarp.time,objwarp.psi(:,itrial)));
        end
        std_f0 = std(objwarp.f, 0, 2);
        std_fn = std(objwarp.fn, 0, 2);
        fgam = zeros(length(TimeOrig),size(FuncOrig,2));
        for itrial = 1:size(FuncOrig,2)
            fgam(:,itrial) = warp_f_gamma(objwarp.fmean,objwarp.gam(:,itrial)',objwarp.time);
        end
        std_fgam = std(fgam, 0, 2);
        objwarp.stats.orig_var = trapz(objwarp.time,std_f0.^2);
        objwarp.stats.amp_var = trapz(objwarp.time,std_fn.^2);
        objwarp.stats.phase_var = trapz(objwarp.time,std_fgam.^2);
        objwarp.qun = [];
        objwarp.origstruct = obj0;
    end


    t2 = toc(tt2);
    objwarp.TimeEl = t2 - t1;
    
else
    Func2Align = FuncOrig;
    
    % Use derivative or integral for alignment. Then, respectively, integrate or
    % differentiate. If the function is velocity, "integral" would preserve total distance
    % traveled, while "derivative" would preserve acceleration values. "Integral" is
    % expected to be less susceptible to noise, while "derivative" is expected to be more
    % accurate but susceptible to noise.
    if strcmpi(opts.EFDAPreserving,'derivative')
        Func2Align = mydiff(TimeOrig,FuncOrig,'NoReshape');
    elseif strcmpi(opts.EFDAPreserving,'integral')
        Func2Align = cumtrapz(TimeOrig,FuncOrig);
    end
    
    % What to do with distance metrics, warps, etc.?
    % Compute them manually and store in the main structure. 
    % But create a substructure with all the info about actual aligned functions (der or
    % int)


    objwarp = fdawarp(Func2Align,TimeOrig); % Creating object
    objwarp = objwarp.time_warping(opts.EFDALambda,option); % Alignment
    
    obj0 = fdawarp2struct(objwarp);
    objwarp = obj0;

    if strcmpi(opts.EFDAPreserving,'derivative')
        % Integrate fn
        objwarp.f = FuncOrig;
        objwarp.fn = cumtrapz(TimeOrig,obj0.fn);
        objwarp.fmean = cumtrapz(TimeOrig,obj0.fmean);
        %objwarp.alignment = 'from the derivative alignment';
    elseif strcmpi(opts.EFDAPreserving,'integral')
        objwarp.f = FuncOrig;
        objwarp.fn = mydiff(TimeOrig,obj0.fn,'NoReshape');
        objwarp.fmean = mydiff(TimeOrig,obj0.fmean);
        %objwarp.alignment = 'from the integral alignment';
    end
    if any(strcmpi(opts.EFDAPreserving,{'derivative','integral'}))
        % recompute q0, qn, mqn, gam, psi, stats, Dist, gamI
        % zero qun
        for itrial = 1:size(FuncOrig,2)       
            objwarp.q0(:,itrial) = f_to_srvf(objwarp.f(:,itrial),TimeOrig);
            objwarp.qn(:,itrial) = f_to_srvf(objwarp.fn(:,itrial),TimeOrig);
            objwarp.gam(:,itrial) = optimum_reparam(objwarp.qn(:,itrial),objwarp.q0(:,itrial),TimeOrig,opts.EFDALambda,opts.EFDAMethod);
            objwarp.psi(:,itrial) = f_to_srvf(objwarp.gam(:,itrial),TimeOrig);
        end
        objwarp.mqn = mean(objwarp.qn,2,'omitnan');
        objwarp.gamI = SqrtMeanInverse(objwarp.gam');
        for itrial = 1:size(FuncOrig,2)
            objwarp.FullDist(itrial) = sqrt(trapz(objwarp.time,(objwarp.fn(:,itrial) - objwarp.fmean).^2));
            objwarp.ElasticAmpDist(itrial) = sqrt(trapz(objwarp.time,(objwarp.mqn - objwarp.qn(:,itrial)).^2));
            objwarp.ElasticPhDist(itrial) = acos(trapz(objwarp.time,objwarp.psi(:,itrial)));
        end
        std_f0 = std(objwarp.f, 0, 2);
        std_fn = std(objwarp.fn, 0, 2);
        fgam = zeros(length(TimeOrig),size(FuncOrig,2));
        for itrial = 1:size(FuncOrig,2)
            fgam(:,itrial) = warp_f_gamma(objwarp.fmean,objwarp.gam(:,itrial)',objwarp.time);
        end
        std_fgam = std(fgam, 0, 2);
        objwarp.stats.orig_var = trapz(objwarp.time,std_f0.^2);
        objwarp.stats.amp_var = trapz(objwarp.time,std_fn.^2);
        objwarp.stats.phase_var = trapz(objwarp.time,std_fgam.^2);
        objwarp.qun = [];
        objwarp.origstruct = obj0;
    end


end




end



function warpstruct = fdawarp2struct(warpobj)
warpstruct = struct;
fldns = fieldnames(warpobj);
for ifldn = 1:length(fldns)
    fldn = fldns{ifldn};
    warpstruct.(fldn) = warpobj.(fldn);
end
end
