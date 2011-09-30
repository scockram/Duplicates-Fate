% Script to handle the processing (duplication and mutation analysis)
% of existing networks. Designed to run in parallel with other processors.
%
% Expects Data/Networks/queue to have be set up by ./preprocessing
%   - to_process
%   - processing
%   - processed
%
% Outline: finds a network to process from the list, and perform some file
% manipulation to keep this list in check. THEN, does analysis on this network.
% THEN, does more file manipulation on the lists in order to tidy them up.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%  Initially: Set simulation and mutation settings here
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Settings
settings_mutation = struct(...
  'number', 20, ...
  'perturb_style', 'relative', ...
  'perturb', 'row', ...
  'perturb_severity', 0.1, ...
  'big_change', 'row', ...
  'big_change_selection', 0, ...
  'delete_or_regen', 'delete');
settings_simulation = struct();

% Seed random number generator
[dummy_variable, host] = system('hostname');
hex = DataHash(strcat(host,datestr(clock)));
seed = mod(hex2dec(hex(1:16)), 2^32);
s = RandStream('mt19937ar','Seed', seed);
RandStream.setDefaultStream(s)

% Directory layout
from_root = 'Data/Networks/';
to_root = 'Data/Dynamics_alt/';
to_process = strcat(from_root, 'to_process.txt');
processing = strcat(from_root, 'processing');
processed = strcat(from_root, 'processed');

% Create mutator object
mut = Mutator(settings_mutation);
mut.mutate_the_duplicate = 0; % Are we doing alternative stuff or not?

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%  Processing control
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Hacky method of do { } while (); in matlab
while true
  % Load list of organisms
  list = textread(to_process, '%s');
  % If there exists nothing that needs processing
   if length(list) == 0
     break
   else
     % Pick random element from list
     el = randi(length(list));
     ele = cell2mat(list(el));
 
     % Remove the element from to_process
     % Note: probably not efficient, but is not relatively slow
     a = []; b = [];
     if el ~= 1
       a = list(1:(el-1));
     end
     if el ~= length(list)
       b = list((el+1):length(list));
     end
     list_new = [a; b];
     fid = fopen(to_process, 'w');
     for i=1:length(list_new)
       fprintf(fid, '%s\n', cell2mat(list_new(i)));
     end
     fclose(fid);
 
     % Add element to processing queue
     fid = fopen(processing, 'a');
     fprintf(fid, '%s\n', ele);
     fclose(fid);
   end

  % Load the organism
  data_file = strcat(from_root, ele, '.mat');
  data = load(data_file);
  wt = data.organism;

  % Useful
  n = wt.size;
  m = settings_mutation.number;

  % Set up variables useful for processing
  D_wt = data.D;
  D_dup = zeros(1, n-1);
  class = zeros(1,n-1);
  D_dup_m = zeros(m,n-1);
  D_wt_m = zeros(m,n-1);
  mut_eff = zeros(m,n-1);
  
  % Loop through all genes
  for i=1:(n-1)
    % Create duplicate, and see what the dynamics of the duplicate is
    dup = wt.duplicategene(i);
    sim_dup = Simulation(dup);

    % Dynamics = ?
    if sim_dup.completed == 1
      D_dup(i) = sim_dup.dynamics;
    else
      D_dup(i) = Inf;
    end

    % Class of protein = ?
    class(i) = wt.protein_categorisations(i);

    % Now perform mutations -- we do mutations even if duplication creates an
    % unstable system -- mutation may make stable again
    for j=1:m
      % Create a mutated duplicate
      [ mut, dup_m ] = mut.mutate(dup);
      sim_dup_m = Simulation(dup_m);

      % Dynamics = ?
      if sim_dup_m.completed == 1
        D_dup_m(j,i) = sim_dup_m.dynamics;
      else
        D_dup_m(j,i) = Inf;
      end

      % Create a mutated WT
      % First: find out which gene the mutant is (WHY?: NOT ALWAYS i)
      g = mut.gene - (mut.gene>dup.duplicated);
      wt_m = wt;
      a = dup_m.k(g,:);
      b = [ a(1:(g-1)) , a((g+1):length(a)) ];
      wt_m.k(g,:) = b;
      sim_wt_m = Simulation(wt_m);

      % Dynamics = ?
      if sim_wt_m.completed == 1
        D_wt_m(j,i) = sim_wt_m.dynamics;
      else
        D_wt_m(j,i) = Inf;
      end

      % Effects of mutation = ?
      eff = sum(sum(dup_m.k - dup.k));
      mut_eff(j,i) = eff;
    end
  end

  % Save data
  save(strcat(to_root,ele,'.mat'), ...
    'D_wt', 'D_dup', 'class', 'D_dup_m', 'D_wt_m', 'mut_eff');

  % add to processed
  fid = fopen(processed, 'a');
  fprintf(fid, '%s\n', ele);
  fclose(fid);

  % remove from processing queue
  p = textread(processing, '%s');
  fid = fopen(processing, 'w');
  for i=1:length(p)
    if ~strcmp(cell2mat(p(i)) , ele)
      fprintf(fid, '%s\n', cell2mat(p(i)));
    end
  end
  fclose(fid);
end
