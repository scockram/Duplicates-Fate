%% In order to pregenerate a large number of networks with some specific
% settings. Just runs a Controller with some specific settings that
% prevents processing beyond the initial network analysis.

function generate_networks()
  % DataHash() gets caught in a fickle when classes are considered
  warning off

  % Create a controller
  controller = Controller();
  
  % Define the callbacks that are required for the Controller to act as
  % a network generator only.
  null_callback = @(varargin) 0;
  controller.fail_callback = @save_failed_network;
  controller.prot_callback = @save_network;
  controller.dup_callback = null_callback;
  controller.mut_callback = null_callback;

  % Sets settings for the simulation
  controller.mutation_settings = controller.mutation_settings;
  controller.simulation_settings.do_duplication = 0;
  controller.organism_settings = struct(...
    'ideal_connectivity', 0.7,...
    'multiple_receptors', 1,...
    'action_on_receptor', 1,...
    'min_response', 0.01); % min response here gets further categorised
                           % in @save_network.

  controller.number_of_tests = 100;

  % Must seed the random number generator with some suitably unique
  % string (Matlab random numbers are not too random)
  [dummy_variable, host] = system('hostname');
  hex = DataHash(strcat(host,datestr(clock)));
  % RandStream accepts numbers < 2^32, and we reduce the number of digits
  % from the hextstream in order to not lose accuracy, large accuracy in
  % matlab is a bit iffy.
  seed = mod(hex2dec(hex(1:16)), 2^32);
  s = RandStream('mt19937ar','Seed', seed);
  RandStream.setDefaultStream(s)

  % Start running the experiment
  controller.start();
end


%% Saves the network in multiple forms:
% Adds #hash of network to list of network hashes
% Saves Network#hash.mat - matlab data file for the network
% Saves Network#hash.txt - plaintext data file for the network
function save_network(organism, r, D)
  % Computes 'unique' hash of the organism
  hash = DataHash(struct(organism));

  % Categorise the organism based on signal response
  if r < 0.01
    return
  elseif r < 0.05
    subdir = '0.01to0.05';
  elseif r < 0.1
    subdir = '0.05to0.10';
  else
    subdir = '0.10plus';
  end
  root = strcat('Networks/',subdir,'/');

  % Save organism as raw matlab data
  save(strcat(root,hash,'.mat'), 'organism', 'r', 'D');

  % Save organism as more readable text
  % variables
  k = regexprep(mat2str(organism.k), ';', ';\n');
  a = mat2str(organism.a);
  d = mat2str(organism.d);
  alpha = mat2str(organism.alpha);
  % Write to file
  fid = fopen(strcat(root,hash,'.txt'), 'w');
  fprintf(fid, '%5f %5f\n%s\n%s\n\n%s\n\n%s\n',r,D,a,d,alpha,k);
  fclose(fid);

  % Append hash to list of valid networks
  fid = fopen(strcat(root,'valids.txt'), 'a');
  fprintf(fid, '%s %s %5f\n', datestr(clock), hash, r);
  fclose(fid);
end


%% Makes a record of failed networks.
% Keep no record, not needed
function save_failed_network(organism, r, D)
end
