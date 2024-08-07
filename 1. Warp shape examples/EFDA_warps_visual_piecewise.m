function EFDA_warps_visual_piecewise(varargin)
% Demonstration of effect that a few examples of piecewise warping functions have on a
% one-dimensional signal.
% Warps are generated with three-to-four knots at predefined locations and piece-wise
% linear segments between them.
% This visualization helped me digest the concept of warping when I started studying this
% approach.

% The original signal is similar to the one in the manuscript, a two-peaked Gaussian. In
% this function, only warping functions "symmetric" w.r.t. the identity line (t=t) are
% plotted, i.e. these warps oscillate around identity. 
% See also the _continuous_asymm.m and _continuous_symm.m functions for the continuous warps

% Aleksei Krotov
% Northeastern University, 2024.


CloseFig = any(strcmpi(varargin,'CloseFig'));
SaveFig = any(strcmpi(varargin,'SaveFig'));
figpos0 = [0 0 0 0]; % Change if plotting not on the main monitor

Nsamp = 200;
time = linspace(0,1,Nsamp);

%a = 10;
% The function is two gaussian peaks having same sigma, and the mean symmetric about 0.5
d_mu = 0.15; %must be <=0.5
sigm = 0.07;
fsource = normpdf(time,0.5-d_mu,sigm) + normpdf(time,0.5+d_mu,sigm);
%figure, plot(time, fsource);

% Use the GaussModel function with the same parameters
Pars = [4;... %A1
    3;... %A2
    0.16;... % B1 Full-width at half-maximum
    0.16;... % B2
    0.4;    % Delta-T time between the two peaks
    0.3]; % Time of the first peak
%Pars = setParams1();
fH = @(tt,Pars) Pars(1) * exp(-4.*log(2) ./ (Pars(3)).^2 .* (tt-Pars(6)).^2) +...
    Pars(2) * exp(-4.*log(2) ./ (Pars(4)).^2 .* (tt-Pars(6)-Pars(5)).^2);
fsource = fH(time,Pars);
fdersource = mydiff(time',fsource);
fintsource = cumtrapz(time',fsource);



%% Warping functions

% 1. compress
alinear = 0.8; % That would mean compress
gamm(:,1) = alinear * time;

% 2. stretch
alinear = 1.2;
gamm(:,2) = alinear * time;

% 3. piecewise slower then faster (convex); t0 is on the Y axis, because f(gamm(t))
t0 = 0.35; % node position
a1 = 0.7; % slope on the first part
a2 = (1 - t0) / (1 - t0/a1); %slope on the second part
b2 = 1 - a2; % intercept of the second part
t0x = t0 / a1;
gamm(:,3) = a1 * time .* (time < t0x) + (a2 * time + b2) .* (time >= t0x);
% a2 = (1 - a1 * t0) / (1-t0); % if t0 was on the X axis.
% gamm(:,3) = a1 * time .* (time < t0) + (a2 * time + b2) .* (time >= t0); % if t0 was on the X


% Just once
% q1 = f_to_srvf(fwarped(:,3),time');
% q2 = f_to_srvf(fsource',time');
% gam = optimum_reparam(q2,q1,time',0.01);
% fw = warp_f_gamma(fwarped(:,3),gam,time.');
% fw = warp_f_gamma(f2,gam,time.');
% [dy, dx] = elastic_distance(fwarped(:,3), fsource', time', 0.01);
% binsize = mean(diff(time));
% psi = sqrt(gradient(gam,binsize));

% 4. piecewise faster then slower (concave)
t0 = 0.65;
a1 = 1.3;
a2 = (1 - t0) / (1 - t0/a1);
b2 = 1 - a2;
t0x = t0 / a1;
gamm(:,4) = a1 * time .* (time < t0x) + (a2 * time + b2) .* (time >= t0x);
% a2 = (1 - a1 * t0) / (1-t0);  % if t0 was on the X axis;
% gamm(:,4) = a1 * time .* (time < t0) + (a2 * time + b2) .* (time >= t0);

% 5. Shift a part towards later time (convex) without warping that part
t0 = 0.05; % shift begins (Y axis)
t1 = 0.61; % shift ends
d = - 0.17; % value of shift (vertical difference with identity); Constraint!  d < t1 - 1;
t0x = t0 - d; % beginning of parallel no-warp shift
a0 = t0 / t0x; % warping of the early segment
a1 = 1; % no warp of the central part
b1 = d;
t1x = t1 - d;
a2 = (1 - t1) / (1 - t1x);
b2 = 1 - a2;
gamm(:,5) = a0 * time .* (time < t0x) + (a1 * time + b1) .* (time >= t0x & time < t1x) + (a2 * time + b2) .* (time >= t1x);

% 6. Shift a part towards earlier time (concave) without warping that part
% the same?
t0 = 0.39; % shift begins (Y axis)
t1 = 0.95; % shift ends
d = 0.17; % value of shift (vertical difference with identity); Constraint!  d < t0;
t0x = t0 - d; % beginning of parallel no-warp shift
a0 = t0 / t0x; % warping of the early segment
a1 = 1; % no warp of the central part
b1 = d;
t1x = t1 - d;
a2 = (1 - t1) / (1 - t1x);
b2 = 1 - a2;
gamm(:,6) = a0 * time .* (time < t0x) + (a1 * time + b1) .* (time >= t0x & time < t1x) + (a2 * time + b2) .* (time >= t1x);

% Compute timeshift and timespeed
gammSpeed = mydiff(time,gamm);
gammShift = gamm - time';


%% Warps

for i = 1:size(gamm,2)
    fwarped(:,i) = warp_f_gamma(fsource,gamm(:,i),time);
end

% Special treatment for compress and stretch
fwarped(:,1) = interp1(time, fsource, (time(end)-time(1)).* alinear .*time + time(1))';%compress
fwarped(:,2) = interp1(time, fsource, (time(end)-time(1)).* (1/alinear) .*time + time(1))';%compress

fder = nan(length(time),6);
fint = nan(length(time),6);
fder(:,3:end) = mydiff(time',fwarped(:,3:end));
fint(:,3:end) = cumtrapz(time',fwarped(:,3:end));

%% Plotting

figname = 'EFDA warpings demo piecewise';
Create_Reuse_Figure([],figname,[200 20 1800 600] + figpos0);

% Original Fn
posOrig = [0.02 0.15 0.06 0.3];
sp.Original = subplot('Position',posOrig); hold on; title('Original','fontsize',14);
%xticklabels({}); xticks('auto');
ylabel('$F(t)$','interpreter','latex');
xlabel('$t$','interpreter','latex');
view(-90,90);
sp.Original.XAxisLocation = 'top';
plot(time, fsource,'k');



if 1
    posgam = [0.0 0.15 0.1 0.3];
    poswarped = [0.0 0.5 0.1 0.1];
    posder = [0.0 0.66 0.1 0.1];
    posint = [0.0 0.82 0.1 0.1];
    
    % all other warps. [0,1] -> [0,1];
    
    for i = 3:size(gamm,2)
        posgam = posgam + [0.13 0 0 0];
        if i>3, posgam = posgam + [0.07 0 0 0]; end
        sp.gam(i) = subplot('Position',posgam); hold on; box on;
        
        xticks([0 1]);yticks ([0 1]);
        if i == 3
            sp.gam(i).YLabel.Interpreter = 'latex';
            PlotMyFormat('Warp Fn',...
                '','','','','','','TitleSize',10);
            ylabel('\gamma','interpreter','tex');
        end
        xlabel('$t$','interpreter','latex');
        plot(time, gamm(:,i),'r');
        plot(time, time, '--k','LineWidth',1);
        
       
        
        % Timespeed and TimeShift
        posgamSp = [posgam(1) 0 0 0] + [0.12 0.2 0.06 0.1];
        sp.gamSp(i) = subplot('Position',posgamSp); hold on; box on;
        plot(time,gammSpeed(:,i),'r');
        yline(1,'LineStyle','--','Color','k','LineWidth',1);
        xticks([0 1]); xlabel('$t$','interpreter','latex');
        if i == 3
            PlotMyFormat('TSpeed',...
                '','','','','','','TitleSize',10);
        end
        
        posgamSh = [posgam(1) 0 0 0] + [0.12 0.33 0.06 0.12];
        sp.gamSh(i) = subplot('Position',posgamSh); hold on; box on;
        plot(time,gammShift(:,i),'r');
        yline(0,'LineStyle','--','Color','k','LineWidth',1);
        xticklabels({}); xticks([0 1]);
        if i == 3
            PlotMyFormat('TShift',...
                '','','','','','','TitleSize',10);
        end
        
        % Warped fn
        poswarped = poswarped + [0.13 0 0 0];
        if i>3, poswarped = poswarped + [0.07 0 0 0]; end
        sp.warped(i) = subplot('Position',poswarped); hold on;
        if i == 3
            sp.warped(i).YLabel.Interpreter = 'latex';
            PlotMyFormat('Orig & Warped',...
                '','$F(\gamma(t))$','','','','','TitleSize',10);
        else
            PlotMyFormat('',...
                '','','','','','');
        end
        xticks([0 1]);
        plot(time, fsource,'k');
        plot(time, fwarped(:,i),'r');
        
         % Add derivative
        posder = posder + [0.11 0 0 0];
        posder(1) = posgam(1);
        sp.der(i) = subplot('Position',posder); hold on;
        if i == 3
            sp.der(i).YLabel.Interpreter = 'latex';
            PlotMyFormat('Derivatives',...
                '','','','','','','TitleSize',10);
        end
        xticklabels({});
        xticks([0 1]);
        plot(time, fdersource,'k');
        plot(time, fder(:,i),'r');
        
        
        % Add integral
        posint = posint + [0.11 0 0 0];
        posint(1) = posgam(1);
        sp.int(i) = subplot('Position',posint); hold on;
        if i == 3
            sp.der(i).YLabel.Interpreter = 'latex';
            PlotMyFormat('Integrals',...
                '','','','','','','TitleSize',10);
        end
        xticklabels({});
        xticks([0 1]);
        plot(time, fintsource,'k');
        plot(time, fint(:,i),'r');
        
        
    end
    
    linkaxes(sp.gamSp,'xy');
    linkaxes(sp.gamSh,'xy');
    linkaxes(sp.der,'xy');
    linkaxes(sp.int,'xy');
    
end




if SaveFig
    savefig(gcf,sprintf('EFDA Analyses/Generic algorithm with examples/%s.fig',figname));
    saveas(gcf,sprintf('EFDA Analyses/Generic algorithm with examples/%s.png',figname));
end

if CloseFig
    close(gcf);
end

end

