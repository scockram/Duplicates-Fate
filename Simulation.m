classdef Simulation
  %SIMULATION Simulates the response of an organism to a ligand input
  
  properties
    % Output variables
    time
    concentrations
    ligand
    response
    dynamics
    completed = 0;
    organism
    
    end_time % when the system reaches stability again
    
    min_response = 0.008; % As an absolute change
    max_time_to_stable = 50;
    
    % Useful accessor
    odeoptions = odeset('AbsTol',1e-4,'RelTol',1e-4);
  end
  
  methods
    % Simulation initializer
    function self = Simulation(organism)
      % Variables needed for processing the simulation
      self.organism = organism;
      self.min_response = organism.init_settings.min_response;
      
      try
        % Initial steady state reaching
        self.ligand = 0;
        trange = 0:600;
        active0 = 0.5 * ones(1, organism.size);
        % Solve and then make variables available for stability analysis
        [t1,y1] = ode23(@organism.odes, trange, active0, self.odeoptions, self.ligand);
        self.time = [ t1 ];
        self.concentrations = [ y1 ];

        % Check if stable now - if not stop
        if check_stability(self) == 0
          error('not stable')
        end

        % Now introduce ligand and solve for atleast 400 steps
        self.ligand = 1;
        trange = 600:1000;
        active0 = y1(size(y1,1),:);
        % Solve and then make variables available for stability analysis
        [t2,y2] = ode23(@organism.odes, trange, active0, self.odeoptions, self.ligand);
        self.time = [ t1; t2 ];
        self.concentrations = [ y1; y2 ];
        
        % remove ligand, run for some time (to stabilise!)
        self.ligand = 0;
        trange = 1000:1050;
        active0 = y2(size(y2,1),:);
        [t3,y3] = ode23(@organism.odes, trange, active0, self.odeoptions, self.ligand);
        self.time = [ t1; t2; t3 ];
        self.concentrations = [ y1; y2; y3 ];

        % Check if stable now - if not stop
        if check_stability(self) == 0
          error('not stable')
        end
        
        % keep going (to ensure it doesn't stop being symbol)
        trange = 1050:1200;
        active0 = y3(size(y3,1),:);
        % Solve and then make variables available for stability analysis
        [t4,y4] = ode23(@organism.odes, trange, active0, self.odeoptions, self.ligand);
        self.time = [ t1; t2; t3; t4 ];
        self.concentrations = [ y1; y2; y3; y4 ];

        % Check if stable now - if not stop
        if check_stability(self) == 0
          error('not stable')
        end
        
        % Assume if we got here it is - so far - completed
        self.completed = 1;
        
      % Here, we catch any times that the code errored - including when it
      % turned out the system was not stable when it should have been. If
      % this is the case the simulation did not complete.
      catch
        self.completed = 0;
      end
      
      if self.completed == 1
        % Diagnostic variables
        [ self.response, self.dynamics ] = self.calcdynamics();

        % If the response was too small, then we say nothing actually
        % happened!
        if self.response < self.min_response
          self.completed = 0;
        end
      end
    end
    
    % This checks the stability of the system AT THIS MOMENT
    function stable = check_stability(self)
      organism = self.organism;
      
      % Find the values of the ODEs at this point.
      concs = self.concentrations(length(self.time),:);
      avg_rate = mean(abs(organism.odes(0, concs, self.ligand)));
      
      % Stable if less than some number
      if avg_rate < 1e-5
        stable = 1;
      end
    end
    
    function [ r, D ] = calcdynamics(self)
      % Reference variables
      p=self.concentrations(:,size(self.concentrations,2));
      t_pre = 600;
      p_pre = p(t_pre);
      min_time_to_stabilise = 20;
      max_time_to_stabilise = 50;
      
      summation = 0; u=0;
      for s=t_pre:(1000 + min_time_to_stabilise)
        % Relating to the integral
        summation = summation + abs(p_pre - p(s));
      end
      
      % More useful reference variables
      t_post = s;
      p_post = p(t_post);

      % relevant outputs
      r = summation / ( t_post - t_pre );
      D = r + p_pre + p_post;
    end
  end
end
