classdef Controller
  %CONTROLLER Controls the whole experiment (all organisms)
  
  properties
    % Because we run the experiment multi-threaded
    thread

    % For generating organisms, hashes of properties
    mutation_settings
    organism_settings
    simulation_settings

    % Callbacks! - these can be modified by the controller controller (...)
    % to allow for certain things to be captured.
    fail_callback
    prot_callback
    dup_callback
    mut_callback

    % Actually related to the whole experiment
    number_of_tests

    % Persistant storage
    D_store
    protein_classification_store
    size_of_organism_store
    connectivity_of_organism_store
    connectivity_of_gene_store
    Ddash_store
    Ddashdash_store
    Ddashdashdash_store
    muteffect_store
  end
    
  methods
    %  set default values here to allow class clearing
    function self = Controller()
      self.mutation_settings = struct(...
        'quantity', 20,...
        'perturb_style', 'relative',...
        'perturb', 'none',...
        'perturb_severity', 0.1,...
        'big_change', 'col',...
        'big_change_selection', 1,...
        'delete_or_regen', 'delete');
      self.organism_settings = struct(...
        'ideal_connectivity', 0.7,...
        'multiple_receptors', 1,...
        'action_on_receptor', 1,...
        'min_response', 0.008);
      self.simulation_settings = struct(...
        'do_duplication',  1);
      self.fail_callback = @(org, D) 0;
      self.prot_callback = @(org, D) 0;
      self.dup_callback  = @(org, gene, org_dup, D, Dd, s) 0;
      self.mut_callback  = @(org, gene, org_dup, org_mut, D, Dd, Ddd, s, sd) 0;
      self.number_of_tests = 100;
    end
    
    function self = start(self)
      % Initialise the storage of output variables
      self.protein_classification_store = zeros(self.number_of_tests, 1);
      self.size_of_organism_store = zeros(self.number_of_tests, 1);
      self.D_store = zeros(self.number_of_tests, 1);
      self.connectivity_of_organism_store = zeros(self.number_of_tests, 1);
      self.Ddash_store = zeros(self.number_of_tests, 15);
      self.connectivity_of_gene_store = zeros(self.number_of_tests, 30);
      self.protein_classification_store = zeros(self.number_of_tests, 30);
      self.Ddashdash_store = zeros(self.number_of_tests, 30,...
          self.mutation_settings.quantity);
      self.Ddashdashdash_store = zeros(self.number_of_tests, 30,...
          self.mutation_settings.quantity);
      self.muteffect_store = zeros(self.number_of_tests, 30,...
          self.mutation_settings.quantity);
      
      % Run simulations until enough organisms have been tested
      test = 1;
      while test <= self.number_of_tests
        org = Organism(self.organism_settings);
        sim = Simulation(org);
        
        % If and only if this is a "valid" simulation, proceed with further
        % analysis
        if sim.completed == 1
          % Store the values relevant to the main network
          self.D_store(test) = sim.dynamics;
          self.size_of_organism_store(test) = org.size;
          self.connectivity_of_organism_store(test) = org.actual_connectivity;
          
          % A bit of matrix manipulation to store protein connectivities.
          a1 = org.protein_connectivities;
          a2 = org.protein_categorisations;
          b = self.connectivity_of_gene_store(test,:);
          filler1 = zeros(1, length(b)-length(a1));
          filler2 = zeros(1, length(b)-length(a2));
          self.connectivity_of_gene_store(test,:) = [ a1 filler1 ];
          self.protein_classification_store(test,:) = [ a2 filler2 ];
          
          % callback as required
          self.prot_callback(org, sim.response, sim.dynamics);
          
          % If we are doing duplications
          if self.simulation_settings.do_duplication
            % Must do a duplication on all genes (not output!)
            for i=1:(org.size-1)
              org_dup = org.duplicategene(i);
              sim_dup = Simulation(org_dup);
              
              % Only care about some things if the duplication simulation
              % completed, some things = mutate it and actually use the
              % dynamics the system claims to have
              if sim_dup.completed ~= 1
                self.Ddash_store(test,i) = Inf;
              else
                self.Ddash_store(test,i) = sim_dup.dynamics;
              end
              
              % callback as required
              self.dup_callback(org, i, org_dup, sim.dynamics, ...
                self.Ddash_store(test,i), ...
                self.fitness(sim.dynamics, self.Ddash_store(test,i)));
                
              % Mutation related things
              for o=1:self.mutation_settings.quantity
                % Create a mutator and set its mutation properties.
                mut = Mutator(self.mutation_settings);
                org_mut = mut.mutate(org_dup);
                sim_mut = Simulation(org_mut);

                % Store the values if relevant!
                if sim_mut.completed == 1
                  self.Ddashdash_store(test,i,o) = sim_mut.dynamics;
                else
                  self.Ddashdash_store(test,i,o) = Inf;
                end

                % Create a mutated wildtype (not mutated duplicate)
                wt_mut = org;
                a = org_mut.k(i,:);
                b = [ a(1:(i-1)) , a((i+1):length(a)) ];
                wt_mut.k(i,:) = b;
                sim_wtm = Simulation(wt_mut);

                % Store the values if relevant!
                if sim_wtm.completed == 1
                  self.Ddashdashdash_store(test,i,o) = sim_wtm.dynamics;
                else
                  self.Ddashdashdash_store(test,i,o) = Inf;
                end

                % Record the effect of the mutation
                eff = org_mut.k(i,:) - org_mut.k(i+1,:);
                self.muteffect_store(test,i,o) = sum(eff);
                
                % callback as required
                self.mut_callback(org, i, org_dup, org_mut, sim.dynamics, ...
                  self.Ddash_store(test,i), ...
                  self.Ddashdash_store(test,i,o), ...
                  self.fitness(sim.dynamics, self.Ddash_store(test,i)), ...
                  self.fitness(sim.dynamics, self.Ddashdash_store(test,i,o)));
              end % End mutation
            end % End duplication
          end
          
          % Was a successful test, counts towards the required quota
          test = test + 1;
        else
          self.fail_callback(org, sim.dynamics);
        end % End of if original test was stable
      end
    end
    
    function s = fitness(self, D,Dd)
      sigma = 0.1;
      s = 1 - exp(- abs(Dd-D)/sigma);
    end
  end
end

