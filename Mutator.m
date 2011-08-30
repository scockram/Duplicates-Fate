classdef Mutator
  %MUTATOR Only changes the interaction matrix.

  properties
    %% Unique Mutation Settings
    % 1 - absolute
    % 0 - relative
    perturb_style = 1;
    
    % Clear
    perturb_severity = 0.1
    
    % 1 - perturb
    % 0 - don't
    perturb_row = 0;
    perturb_col = 0;
    
    %  1 - do one only
    %  0 - don't
    % -1 - randomly select big changes
    big_change_row = 0;
    big_change_col = 0;
    
    % how it chooses which interaction to play with
    %  1 - purely random choice 
    %  0 - 50/50 from absent/existent things, then random from each of those
    % -1 - only existing
    % -2 - only non-existing
    big_change_selection_method = -1
    
    % 1 - regen
    % 0 - delete
    big_change_delete_regen = 0;

    % 1 - mutation is on the duplicated gene
    % 0 - is not
    mutate_the_duplicate = 1;
    
    %% Interals
    gene
  end
  
  methods
    %% Helper functions
    
    % Creates the mutator, sets the settings based on the input struct.
    function self = Mutator(settings)
      % Perturbations
      switch settings.perturb
        case 'both'
          self.perturb_row = 1;
          self.perturb_col = 1;
        case 'row'
          self.perturb_row = 1;
        case 'col'
          self.perturb_col = 1;
      end
      self.perturb_severity = settings.perturb_severity;
      switch settings.perturb_style
        case 'absolute'
          self.perturb_style = 1;
        case 'relative'
          self.perturb_style = 0;
      end
      
      % Random big change
      switch settings.big_change
        case 'both'
          self.big_change_row = 1;
          self.big_change_col = 1;
        case 'row'
          self.big_change_row = 1;
        case 'col'
          self.big_change_col = 1;
      end
      self.big_change_selection_method = settings.big_change_selection;
      switch settings.delete_or_regen
        case 'regen'
          self.big_change_delete_regen = 1;
        case 'delete'
          self.big_change_delete_regen = 0;
      end
    end
    
    % Mutates the organism - calls the small mutation methods individually
    % based on the settings that were given at creation
    function org = mutate(self, org)
      % useful
      if mutate_the_duplicate
        self.gene = org.duplicated;
      else
        y = randi(org.size-1);
        self.gene = y + (y>=org.duplicated);
      end
      
      % Repetitive code, but does the job fine.
      if self.perturb_row == 1
        org.k = self.do_perturb_row(org);
      end
      
      if self.perturb_col == 1
        org.k = self.do_perturb_col(org);
      end
      
      if self.big_change_row == 1
        org.k = self.do_big_change(org, 'row');
      end
      
      if self.big_change_col == 1
        org.k = self.do_big_change(org, 'col');
      end
    end
    
    %% Actual mutation functions 
    function k = do_perturb_row(self, org)
      k = org.k;
      
      % Choose the "base" - absolute or relative perturbation
      if self.perturb_style == 1
        base = k(self.gene, :) ~= 0;
      else
        base = k(self.gene, :);
      end
      
      % apply it
      k(self.gene, :) = k(self.gene, :) + ...
        base * self.perturb_severity .* (2*rand(1,org.size) - 1);
    end
    
    function k = do_perturb_col(self, org)
      k = org.k;
      
      % Choose the "base" - absolute or relative perturbation
      if self.perturb_style == 1
        base = k(:, self.gene) ~= 0;
      else
        base = k(:, self.gene);
      end
      
      % apply it
      k(:, self.gene) = k(:, self.gene) + ...
        base * self.perturb_severity .* (2*rand(org.size,1) - 1);
    end
    
    % DRY - transpose the matrix and apply a random row mutation if a
    % column mutation is what was required.
    function k = do_big_change(self, org, orientation)
      k = org.k;
      if orientation == 'col'
        k = k';
      end
      
      % Create a list of possibly mutatable indices of the row
      switch self.big_change_selection_method
        case 1
          possible = find(k(self.gene, :) < 100);
        case 0
          % If there exists 1 or 0 possible interactions then DO THE OTHER THING
          % Note: while loop is a pretty lazy way of doing this, but is easiest,
          % and doesn't cause any massive inefficiencies.
          possible = [];
          while length(possible) < 2
            if rand() < 0.5
              possible = find(k(self.gene, :) ~= 0);
            else
              possible = find(k(self.gene, :) == 0);
            end
          end
        case -1
          possible = find(k(self.gene, :) ~= 0);
        case -2
          possible = find(k(self.gene, :) == 0);
      end
      
      % Pick something random from the possible options
      chosen = randi(length(possible), 1);
      
      % The class of teh /acting/ protein is based upon whether we are
      % acting on rows or close
      if orientation == 'col'
        % This is horrible and hacky, the -1 is needed because this vector
        % isn't updated upon duplication - but this is the only reason it
        % woudl be required, so it's done here. TODO: refactor if
        % necessary.
        if chosen > self.gene
          class = org.protein_categorisations(chosen-1);
        else
          class = org.protein_categorisations(chosen);
        end
      else
        class = org.protein_categorisations(self.gene);
      end
      
      % If it exists, regen or delete
      if k(self.gene, chosen) ~= 0
        % Regen
        if self.big_change_delete_regen == 1
          switch class
            case 1
              k(self.gene, chosen) = rand();
            case 2
              k(self.gene, chosen) = rand();
            case 3
              k(self.gene, chosen) = - rand();
            case 4
              k(self.gene, chosen) = 2*rand()-1;
          end
        % Delete
        else
          k(self.gene, chosen) = 0;
        end
      % If it doesn't, create
      else
          switch class
            case 1
              k(self.gene, chosen) = rand();
            case 2
              k(self.gene, chosen) = rand();
            case 3
              k(self.gene, chosen) = - rand();
            case 4
              k(self.gene, chosen) = 2*rand()-1;
          end
      end % End of if ==0 or not thing
      
      % Transpose it again if required
      if orientation == 'col'
        k = k';
      end
    end
  end
  
end

