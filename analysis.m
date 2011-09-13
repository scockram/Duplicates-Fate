%% Load the data
clear;
datafile = 'Data/Dynamics/0.10plus/data.mat';
load(datafile);
[ cdf, tidyplot, histnorm, heatmap, normv, pdfp ] = helpers();
n2s = @num2str;

%% Data Cleaning
% Selecting data which is interesting
v = find( (swt ~= 1) & (sdup ~=1) ) ;
v2 = find( (swt ~= 1) & (sdup == 1) );
v3 = find( (swt == 1) & (sdup ~= 1) );
v4 = find( (swt == 1) & (sdup == 1) );

disp('------');
disp([ 'Datafile: ' datafile ]);
disp([ 'Data set size (Number of mutations): ' n2s(length(swt)) ]);
disp([ 'Cases where instability in WT and DUP: ' n2s(length(v4)) ]);
disp([ 'Cases where instability ONLY in WT: ' n2s(length(v3)) ]);
disp([ 'Cases where instability ONLY in DUP: ' n2s(length(v2)) ]);
disp([ 'Other cases: ' n2s(length(v)) ]);
disp('------');

%% Multiple CDF's

plot(x,y,x1,y1)
tidyplot('Network Fitness Effect')
??? Input argument "y_l" is undefined.

Error in ==> helpers>tidyplot_helper at 78
  ylabel( y_l );
 
tidyplot('Network Fitness Effect', '')
tidyplot('Network Fitness Effect', 'Proportion of Data')

%% PDF of the actual mutation effects
[X,Y] = pdfp(muteff,51);
plot(X,Y);

%% Scatter Plot (Not too revealing)
plot(swt(v), sdup(v), '.', 'MarkerSize', 4);
tidyplot('s_{wt}', 's_{dup}')

%% Heat Map of swt vs sdup
r = linspace(0,1,101);
[x,y,z] = heatmap(swt(v), sdup(v), r, r) ;
z = log(z+1);
contourf(x,y,z)
tidyplot('s_{wt}', 's_{dup}')

%% PDF plot
%sdiff = sdup - swt;
sdiff = sdup(v) - swt(v);
eps = .001;
sdiff2 = normv(sdiff(abs(sdiff) < eps));
n=101;
x = linspace(-eps,eps,n);
y = zeros(1,n);
for i=1:length(sdiff2)
  j = 1+floor((n-1)*sdiff2(i));
  y(j) = y(j) +  1;
end
y = y/length(sdiff2);
y2 = zeros(1,n);
for i=1:length(y);
  if i==1
    y2(i) = 0.9*y(i) + 0.1*y(i+1);
  elseif i==length(y)
    y2(i) = 0.9*y(i) + 0.1*y(i-1);
  else
    y2(i) = 0.8*y(i) + 0.1*y(i-1) + 0.1*y(i+1);
  end
end
plot(x,y,'x')
tidyplot('s_{dup} - s_{wt}', 'Data Density');

%% Qualitative Numbers for the PDF plot at different "zooms"
bp = [];
deltalabels = {};
deltas = (10.^(-(1:6)));
for delta=deltas
  sdiffr = sdiff(abs(sdiff) < delta);
  epsilon = delta/(10^3);
  a = [sum(sdiffr<-epsilon) sum(abs(sdiffr<=epsilon)) sum(sdiffr>epsilon)];
  bp = [ bp; a ];
  deltalabels = [ deltalabels {[ num2str(delta) ]} ];
end
normdata = (bp ./ [ sum(bp,2) sum(bp,2) sum(bp,2) ]);
bar(normdata, 'stacked')
legend('s_{diff} < -\epsilon', '|s_{diff}| < \epsilon', 's_{diff} > \epsilon')
set(gca, 'XTickLabel', deltalabels)
tidyplot('\delta', 'Proportion of Data')
title('\epsilon = \delta / 10^3, looking at data in [-\delta, \delta]')

%% auxillary plot to the previous
plot(normdata, '--.')
set(gca, 'xtick', 1:length(deltas))
set(gca, 'xticklabel', deltalabels)
legend('s_{diff} < -\epsilon', '|s_{diff}| < \epsilon', 's_{diff} > \epsilon')
tidyplot('\delta', 'Proportion of Data')

%% Muteff vs sdiff heat map
[x,y,z] = heatmap(normv(sdiff), normv(muteff), ...
  linspace(-1, 1, 151), linspace(min(muteff), max(muteff), 151));
zz = log(z+1);
plot(sdiff, muteff)
tidyplot('s_{dup} - s_{wt}', 'Mutation Effect');
