%% In order to pregenerate a large number of networks with some specific
% settings. Just runs a Controller with some specific settings that
% prevents processing beyond the initial network analysis.

function generate_networks()
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
  % -- str2uniquenum

  % Start running the experiment
  controller.start();
end


%% Saves the network in multiple forms:
% Adds #hash of network to list of network hashes
% Saves Network#hash.mat - matlab data file for the network
% Saves Network#hash.txt - plaintext data file for the network
function save_network(organism, r, D)
  % --
  disp([ r D ])
end


%% Makes a record of failed networks.
% Only stores a hash of the network - just to keep a record.
function save_failed_network(organism, r, D)
end
