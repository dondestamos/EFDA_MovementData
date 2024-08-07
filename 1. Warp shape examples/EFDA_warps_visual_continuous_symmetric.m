function EFDA_warps_visual_continuous_symmetric(varargin)

% When amplitude gets too large (examine for each warp), the pattern may change in terms
% of the "fixed" (gamma = t) points and in terms of points with largest d_gamma/dt.

CloseFig = any(strcmpi(varargin,'CloseFig'));
SaveFig = any(strcmpi(varargin,'SaveFig'));

Nsamp = 200;
time = linspace(0,1,Nsamp);
figpos0 = [0 0 0 0]; % Change if plotting not on the main monitor
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
fsource = fH(time,Pars)';

fdersource = mydiff(time',fsource);
fintsource = cumtrapz(time',fsource);


%% Warping functions

% symmetric about identity, .5 fixed, sigmoid
i = 1; %convex first
[psi(:,i), gamm(:,i), shvec(:,i)] = genwarp(0.3,1,time,'InverseFlag',1,'Shift',pi/2);

i = 2; %concave first
[psi(:,i), gamm(:,i), shvec(:,i)] = genwarp(0.3,1,time,'InverseFlag',0,'Shift',pi/2);

% symmetric about identity, .33 and .67 fixed
i = 3;
[psi(:,i), gamm(:,i), shvec(:,i)] = genwarp(0.3,1.5,time,'InverseFlag',1,'Shift',pi/2);

i = 4;
[psi(:,i), gamm(:,i), shvec(:,i)] = genwarp(0.3,1.5,time,'InverseFlag',0,'Shift',pi/2);

% symmetric about identity, .25, .5, .75 fixed
i = 5;
[psi(:,i), gamm(:,i), shvec(:,i)] = genwarp(0.3,2,time,'InverseFlag',1,'Shift',pi/2);

i = 6;
[psi(:,i), gamm(:,i), shvec(:,i)] = genwarp(0.3,2,time,'InverseFlag',0,'Shift',pi/2);

% symmetric about identity, .2, .4, .6, .8 fixed
i = 7;
[psi(:,i), gamm(:,i), shvec(:,i)] = genwarp(0.3,2.5,time,'InverseFlag',1,'Shift',pi/2);

i = 8;
[psi(:,i), gamm(:,i), shvec(:,i)] = genwarp(0.3,2.5,time,'InverseFlag',0,'Shift',pi/2);

% Compute timeshift and timespeed
gammSpeed = mydiff(time,gamm);
gammShift = gamm - time';

%% Warps

for i = 1:size(gamm,2)
    fwarped(:,i) = warp_f_gamma(fsource,gamm(:,i),time')';
    
    % Norm preservving?
    % fwarped(:,i) = fwarped(:,i) .* psi(:,i);
    
    % Area preserving
    % fwarped(:,i) = fwarped(:,i) .* psi(:,i).^2;
    
end

fder = mydiff(time',fwarped);
fint = cumtrapz(time',fwarped);

%% Plotting

figname = 'EFDA warpings demo continuous symmetric';
Create_Reuse_Figure([],figname,figpos0 + [100 20 1800 600]);

% Old option, just a panel of original, warps, warped, and psi-functions
if 0
    % Original Fn
    posOrig = [0.02 0.4 0.04 0.32];
    sp.Original = subplot('Position',posOrig); hold on; title('Original','fontsize',14);
    %xticklabels({}); xticks('auto');
    ylabel('$F(t)$','interpreter','latex');
    xlabel('$t$','interpreter','latex');
    view(-90,90);
    sp.Original.XAxisLocation = 'top';
    plot(time, fsource,'k');
    
    posgam = [0.0 0.4 0.08 0.32];
    poswarped = [0.0 0.8 0.08 0.15];
    pospsi = [0.0 0.15 0.08 0.15];
    % all other warps. [0,1] -> [0,1];
    
    for i = 1:size(gamm,2)
        posgam = posgam + [0.1 0 0 0];
        sp.gam(i) = subplot('Position',posgam); hold on; box on;
        %     xlabel('$t$','interpreter','latex');
        xticks([0 1]);yticks ([0 1]);
        if i == 1, ylabel('\gamma','interpreter','tex'); end
        plot(time, gamm(:,i),'r');
        plot(time, time, '--k','LineWidth',1);
        
        
        poswarped = poswarped + [0.1 0 0 0];
        sp.warped(i) = subplot('Position',poswarped); hold on;
        if i == 1
            sp.warped(i).YLabel.Interpreter = 'latex';
            PlotMyFormat('',...
                '','$F(\gamma(t))$','','','','');
        else
            PlotMyFormat('',...
                '','','','','','');
        end
        xticks([0 1]);
        plot(time, fsource,'k');
        plot(time, fwarped(:,i),'r');
        
        
        pospsi = pospsi + [.1 0 0 0];
        sp.psi(i) = subplot('Position',pospsi); hold on;
        %xlabel('$t$','interpreter','latex');
        sp.psi(i).XLabel.Interpreter = 'latex';
        if i == 1
            sp.psi(i).YLabel.Interpreter = 'latex';
            PlotMyFormat('',...
                '$t$','$\psi$','','','','');
        else
            PlotMyFormat('',...
                '$t$','','','','','');
        end
        xticks([0 1]);
        plot(time, psi(:,i),'r');
        yline(1,'Color','k','LineStyle','--','LineWidth',1);
    end
end




% Smaller panels. But adding TimeShift and TimeSpeed
if 1
    % Original Fn
    posOrig = [0.02 0.3 0.04 0.142];
    sp.Original = subplot('Position',posOrig); hold on; title('Original','fontsize',14);
    %xticklabels({}); xticks('auto');
    ylabel('$F(t)$','interpreter','latex');
    xlabel('$t$','interpreter','latex');
    view(-90,90);
    sp.Original.XAxisLocation = 'top';
    plot(time, fsource,'k');
    
    posgam = [0.0 0.3 0.0533 0.142];
    poswarped = [0.0 0.5 0.0533 0.1];
    posder = [0.0 0.65 0.0533 0.1];
    posint = [0.0 0.8 0.0533 0.1];
    
    pospsi = [0.0 0.1 0.0533 0.07];
    % all other warps. [0,1] -> [0,1];
    
    for i = 1:size(gamm,2)
        posgam = posgam + [0.11 0 0 0];
        sp.gam(i) = subplot('Position',posgam); hold on; box on;
        %     xlabel('$t$','interpreter','latex');
        xticks([0 1]);yticks ([0 1]);
        plot(time, gamm(:,i),'r');
        plot(time, time, '--k','LineWidth',1);
        if i == 1
            PlotMyFormat('Warp Fn',...
                '','$F(\gamma(t))$','','','','','TitleSize',10);
        end
        if i == 1, ylabel('\gamma','interpreter','tex'); end
        
        
        
        poswarped = poswarped + [0.11 0 0 0];
        sp.warped(i) = subplot('Position',poswarped); hold on;
        if i == 1
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
        sp.der(i) = subplot('Position',posder); hold on;
        if i == 1
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
        sp.int(i) = subplot('Position',posint); hold on;
        if i == 1
            sp.der(i).YLabel.Interpreter = 'latex';
            PlotMyFormat('Integrals',...
                '','','','','','','TitleSize',10);
        end
        xticklabels({});
        xticks([0 1]);
        plot(time, fintsource,'k');
        plot(time, fint(:,i),'r');
        
        
        
        
        
        
        pospsi = pospsi + [.11 0 0 0];
        sp.psi(i) = subplot('Position',pospsi); hold on;
        %xlabel('$t$','interpreter','latex');
        sp.psi(i).XLabel.Interpreter = 'latex';
        if i == 1
            sp.psi(i).YLabel.Interpreter = 'latex';
            PlotMyFormat('Shooting vector',...
                '$t$','$\psi$','','','','');
            sp.psi(i).Title.FontSize = 10;
        else
            PlotMyFormat('',...
                '$t$','','','','','');
        end
        xticks([0 1]);
        plot(time, shvec(:,i),'r');
        %yline(1,'Color','k','LineStyle','--','LineWidth',1);
        
        
        % Timespeed and TimeShift
        posgamSp = [posgam(1) 0 0 0] + [0.07 0.32 0.028 0.05];
        sp.gamSp(i) = subplot('Position',posgamSp); hold on; box on;
        plot(time,gammSpeed(:,i),'r');
        yline(1,'LineStyle','--','Color','k','LineWidth',1);
        if i == 1
            PlotMyFormat('TSpeed',...
                '','','','','','','TitleSize',8,'FontSize',8);
        end
        xticks([0 1]); xlabel('$t$','interpreter','latex');
        
        posgamSh = [posgam(1) 0 0 0] + [0.07 0.4 0.028 0.05];
        sp.gamSh(i) = subplot('Position',posgamSh); hold on; box on;
        plot(time,gammShift(:,i),'r');
        yline(0,'LineStyle','--','Color','k','LineWidth',1);
        xticklabels({}); xticks([0 1]);
        if i == 1
            sp.psi(i).YLabel.Interpreter = 'latex';
            PlotMyFormat('TShift',...
                '','','','','','','TitleSize',8,'FontSize',8);
            %sp.psi(i).Title.FontSize = 10;
        end
        
        fprintf('%.2f  ',L2norm(fwarped(:,i)));
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