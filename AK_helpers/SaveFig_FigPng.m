function SaveFig_FigPng(varargin)

%example
%SaveFig_FigPng('FigName',figname,'Path','Analysis vars/Curvature related');

if any(strcmpi(varargin,'Handle'))
    fighand = varargin{find(strcmpi(varargin,'Handle'))+1};
else
    fighand = [];
end

if any(strcmpi(varargin,'FigName'))
    figname = varargin{find(strcmpi(varargin,'FigName'))+1};
else
    figname = 'SavedFig';
end

if any(strcmpi(varargin,'Path'))
    pathname = varargin{find(strcmpi(varargin,'Path'))+1};
else
    pathname = '';
end

if any(strcmpi(varargin,'Replace')) %Do not rename, just replace existing (if any)
    flag_replace = 1;
else
    flag_replace = 0;
end

flag_FigSave = 1;
if any(strcmpi(varargin,'NoFig')) %Do not rename, just replace existing (if any)
    flag_FigSave = 0;
end


if length(findobj('type','figure')) < 1
    disp('No Open Figure');
    return
end

if isempty(fighand)
    f = gcf;
else
    f = fighand;
end

c = clock;
c = sprintf('%d-%d-%d-%d',c(2),c(3),c(4),c(5));

if strcmp(figname,'0')
    figname = f.Name;
end

if isfile(sprintf('%s/%s',pathname,figname)) && ~flag_replace
    figname = sprintf('%s-%s',figname,c);
end

if flag_FigSave
savefig(f,sprintf('%s/%s.fig',pathname,figname));
end
saveas(f,sprintf('%s/%s.png',pathname,figname));

end