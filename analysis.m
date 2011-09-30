%% Load the data
clear; clc
for fg=0:4
datafile = 'Data/Dynamics/0.10plus/data.mat';
load(datafile);
[ cdf, tidyplot, histnorm, heatmap, normv, pdfp ] = helpers();
n2s = @num2str;
groups = {'All Data', 'Receptors', 'Kinases', 'Phosphotases', 'Effectors'};

% Filter functional groups
if fg~=0
  r = find(class==fg);
  vars = {'class','muteff','s','sd','sdup','swt'};
  for i=1:length(vars)
    eval([ vars{i} '=' vars{i} '(r);' ])
  end
end

%% Data Cleaning
% Selecting data which is interesting
v  = find( (swt ~= 1) & (sdup ~=1) ) ;
v2 = find( (swt ~= 1) & (sdup == 1) );
v3 = find( (swt == 1) & (sdup ~= 1) );
v4 = find( (swt == 1) & (sdup == 1) );

% disp('------');
% disp([ 'Datafile: ' datafile ]);
% disp([ 'Data set size (Number of mutations): ' n2s(length(swt)) ]);
% disp([ 'Cases where instability in WT and DUP: ' n2s(length(v4)) ]);
% disp([ 'Cases where instability ONLY in WT: ' n2s(length(v3)) ]);
% disp([ 'Cases where instability ONLY in DUP: ' n2s(length(v2)) ]);
% disp([ 'Other cases: ' n2s(length(v)) ]);
% disp('------');

%% PDF plot
%sdiff = sdup - swt;
rangel = [ 0  .1 .2 .4 .6 .8 0 ];
rangeu = [ .1 .2 .4 .6 .8  1 1 ];
outer = [];
for i=1:length(rangel)
  %subplot(6,1,i)
  rb = rangel(i);
  ru = rangeu(i);
  m = (swt(v) >= rb) & (swt(v) < ru);
  sdiff = sdup(v) - swt(v);
  %[x,y] = pdfp(sdup(m), 101);
  %plot(x,y)
  %axis([-1 1 0 0.7]);
  %tidyplot('','')
  dataintotable = sdiff(m);
  dataintotable = sd(v) - swt(v);
  dataintotable = dataintotable(m);
  outer = [ outer; [ rb ru mean(dataintotable) median(dataintotable) ] ];
end
disp(groups{fg+1})
disp(outer')
%tidyplot('s_{dup} - s_{wt}', 'Data Density');


end

%% Replicate Original Analysis
clf; hold on;
colors = [ 1 0 0; 0 0 1; 0 .7 .5; 0 0 0; ];
for c=1:4
  [x,y] = cdf(sd(class==c));
  plot(x,y,'color',colors(c,:))
end
legend('Receptors', 'Kinases', 'Phosphotases', 'Effectors', 'Location','NorthWest')
tidyplot('Network Fitness Effect','Cumulative Frequency')
title('Distribution of Network Fitness Effect of Duplication and Mutation')

%% Multiple CDF's (wt and dup)
[x,y] = cdf(swt);
[x1,y1] = cdf(sdup);
plot(x,y,'b',x1,y1,'r')
tidyplot('Network Fitness Effect', 'Proportion of Data')
legend('s_{wt}', 's_{dup}')

%% PDF of the actual mutation effects
[X,Y] = pdfp(muteff,51);
plot(X,Y);

%% Scatter Plot (Not too revealing)
plot(swt(v), sdup(v), '.', 'MarkerSize', 4);
tidyplot('s_{wt}', 's_{dup}')

%% Heat Map of swt vs sdup
r = linspace(0,1,101);
[x,y,z] = heatmap(swt(v), sdup(v), r, r) ;
z = log(z+1)';
contourf(x,y,z)
tidyplot('s_{wt}', 's_{dup}')

%% Qualitative Numbers for the PDF plot at different "zooms"
bp = [];
deltalabels = {};
deltas = (10.^(-(0:0.5:9)));
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
[x,y,z] = heatmap(normv(sdiff), normv(muteff(v)), ...
  linspace(-1, 1, 151), linspace(min(muteff), max(muteff), 151));
zz = log(z+1);
plot(sdiff, muteff)
tidyplot('s_{dup} - s_{wt}', 'Mutation Effect');
