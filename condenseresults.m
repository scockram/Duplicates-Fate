% From Data/Dynamics/*/*.mat to one data file containing six equal length
% vectors for creating figures and plots from.

% This can be changed
dir = 'Data/Dynamics/0.10plus/';

% Finds what .mat files exist in the dir
stuffindir = what(dir);
datafiles = stuffindir.mat;
nd = length(datafiles);

% initialise storage media
s = zeros(1,250000);
sd = zeros(1,250000);
swt = zeros(1,250000);
sdup = zeros(1,250000);
class = zeros(1,250000);
muteff = zeros(1,250000);
p = 1;

% For calculating s values
sc =  @(s1, s2) 1 - exp(-abs(s1-s2)/0.1); % sigma defined as 0.1

% Loop through datafiles doing some manipulation
for i=1:nd
  tic
  df = cell2mat(datafiles(i));
  file = strcat(dir, df);
  d = load(file);
  
  for j=1:length(d.D_dup)
    for k=1:size(d.D_dup_m, 1)
      s(p) = sc(d.D_wt, d.D_dup(j));
      sd(p) = sc(d.D_wt, d.D_dup_m(k,j));
      swt(p) = sc(d.D_wt, d.D_wt_m(k,j));
      sdup(p) = sc(d.D_dup(j), d.D_dup_m(k,j));
      class(p) = d.class(j);
      muteff(p) = d.mut_eff(k,j);
      
      p = p+1;
    end
  end
  
  disp(i);
end

s = s(1:(p-1));
sd = sd(1:(p-1));
swt = swt(1:(p-1));
sdup = sdup(1:(p-1));
class = class(1:(p-1));
muteff = muteff(1:(p-1));

save(strcat(dir,'data.mat'), 's', 'sd', 'swt', 'sdup', 'class', 'muteff');





