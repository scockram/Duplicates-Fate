%% Loads helpers from this one file into the global namespace. Typically
% called as per the function line:
%   [ cdf ] = helpers()
% No inputs required
function [ cdf, tidyplot, histnorm, heatmap, normv, pdfp ] = helpers()
  cdf = @cdf_helper;
  tidyplot = @tidyplot_helper;
  histnorm = @histnorm_helper;
  heatmap = @heat_helper;
  normv = @normv_helper;
  pdfp = @pdf_helper;
end

%% Returns the x and y axis for a pdf plot
function [x,y] = pdf_helper(data, n)
  x = linspace(min(data), max(data), n);
  y = zeros(1,n);
  data = normv_helper(data);
  
  for i = 1:length(data)
    j = 1+floor((n-1)*data(i));
    y(j) = y(j) +  1;
  end
  y = y / length(data);
end

%% Returns the normalized vector
function normalized = normv_helper(x)
  normalized = (x-min(x))/(max(x)-min(x));
end

%% Heatmap
function [X,Y,Z] = heat_helper(x,y,xrange,yrange)
  % Preprocessing
  n = length(xrange)-1;
  [ X, Y ] = meshgrid(xrange, yrange);
  Z = X - X;
  
  % Categorise all the things
  for i = 1:length(x)
    xr = 1 + floor(n*x(i));
    yr = 1 + floor(n*y(i));
    Z(xr, yr) = Z(xr, yr) + 1;
  end
end


%% Computes the cumulative distribution function of some distribution.
% Output is suitable for plotting
% Alex Podgaetsky, September 2003 - alex@wavion.co.il
function [Xplot,Yplot] = cdf_helper(X)
  tmp = sort(reshape(X,prod(size(X)),1));
  Xplot = reshape([tmp tmp].',2*length(tmp),1);

  tmp = [1:length(X)].'/length(X);
  Yplot = reshape([tmp tmp].',2*length(tmp),1);
  Yplot = [0; Yplot(1:(end-1))];
end


%% Tidies the currently focused figure, and sets labels if variables
% x_l and y_l are set, otherwise prompts for input
function tidyplot_helper(x_l, y_l)
  set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1         );

  %if ~exist('x_l','var')
  %    x_l = input('X-Label: ', 's');
  %    y_l = input('Y-Label: ', 's');
  %end

  xlabel( x_l );
  ylabel( y_l );

  set(gcf, 'PaperPositionMode', 'auto');
end


%% Normalised histogram generator
function varargout = histnorm_helper(varargin)
  doPlot = 0;
  if ischar (varargin{end}) && strcmpi ('plot', varargin{end})
      doPlot = 1;
      varargin = varargin(1:end-1);
  elseif nargout == 0
      doPlot = 1;    
  end

  % normalize so the "area" is 1
  [xo,no] = hist (varargin{:});
  binwidths = diff ([ min(varargin{1}),
                      no(1:end-1)+diff(no)/2,
                      max(varargin{1}) ]);
  xonorm = xo/sum (xo .* binwidths);
  varargout = {xonorm, no};
  varargout = varargout(1:nargout);

  % do plot
  if doPlot
    cax = axescheck(varargin{:});

    % modify vertices of bar plot
    hist (varargin{:});
    if isempty (cax)
        cax = gca;
    end
    ch = findobj (get (cax, 'children'), 'type', 'patch'); ch = ch(1);
    vertices = get (ch, 'vertices');
    for idx = 1:numel (xonorm)
        vertices((idx-1)*5+[3 4],2) = xonorm(idx);     % hope it works :)
    end
    set (ch, 'vertices', vertices);
  end
end
