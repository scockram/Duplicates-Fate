classdef Organism
  %ORGANISM Describes the interactions within the organism
  
  properties
    % Interaction coefficients
    k
    a
    d
    alpha
    
    % Useful to save
    init_settings
    actual_connectivity
    protein_connectivities
    protein_categorisations
    
    % convenience variables
    rs
    ks
    es
    ps
    size
    
    % useful
    duplicated
    duplicated_type
  end
  
  methods
    % Organism initializer
    function self = Organism(settings)
      % Store this
      self.init_settings = settings;
      
      % Convenience variables
      if settings.multiple_receptors == 1
        self.rs = 1 + floor(rand()*6);  % No. of receptor
      else
        self.rs = 1;
      end
      self.ks = floor(rand()*7);      % kinase
      self.ps = floor(rand()*7);      % phosphotase
      self.es = 1 + floor(rand()*6);  % effector
      self.size = self.rs+self.ks+self.ps+self.es+1; % +1 = P_out
      self.protein_connectivities = zeros(1,self.size);
      self.protein_categorisations = zeros(1,self.size);
      
      % Set tge initial atrributes as 0
      self.k = zeros(self.size, self.size);
      self.a = zeros(1, self.size);
      self.d = zeros(1, self.size);
      self.alpha = zeros(1, self.size);
      
      % Interaction matrix And alpha
      for i=1:self.rs % receptors
          self.k(i,:) = self.genrow(i);
          self.alpha(i) = rand()*2 - 1;
          self.protein_connectivities(i) = sum(self.k(i,:)~=0)/(length(self.k(i,:))-1);
          self.protein_categorisations(i) = 1;
      end
      for i=(self.rs+1):(self.rs+self.ks) % kinases
          self.k(i,:) = self.genrow(i);
          self.protein_connectivities(i) = sum(self.k(i,:)~=0)/(length(self.k(i,:))-1);
          self.protein_categorisations(i) = 2;
      end
      for i=(self.rs+self.ks+1):(self.rs+self.ks+self.ps) % phosphotases
          self.k(i,:) = -1 * self.genrow(i);
          self.protein_connectivities(i) = sum(self.k(i,:)~=0)/(length(self.k(i,:))-1);
          self.protein_categorisations(i) = 3;
      end
      for i=(self.rs+self.ks+self.ps+1):(self.rs+self.ks+self.ps+self.es) % effectors
          self.k(i, self.size) = rand()*2 - 1;
          self.protein_connectivities(i) = sum(self.k(i,:)~=0)/(length(self.k(i,:)));
          self.protein_categorisations(i) = 4;
      end
      
      % Relaxation rates
      for i=1:(self.size-1)
          rate = (rand() * 2 - 1) * 0.1;
          if rate<0
              self.d(i) = rate;
          else
              self.a(i) = rate;
          end
      end
      
      % After generating k, find out what the connectivity of the network
      % actually is (so ignore effectors here)
      s = self.k(1:(self.size-1-self.es),1:(self.size-1));
      self.actual_connectivity = 1 - ( sum(sum(s==0))-size(s,1) ) / ( size(s,1)*size(s,2) );
    end
    
    % For generating non-effector values
    function row = genrow(self, rowid)
        row = zeros(1, self.size);
        for i=1:(self.size-1)
            if (i ~= rowid) && (rand() < self.init_settings.ideal_connectivity)
                row(i) = rand();
            end
        end
    end
    
    % Duplicates one of the genes of the organism
    function organism = duplicategene(organism, gene)
      % Housekeeping for future mutations
      organism.duplicated = gene;
      if gene <= organism.rs
        organism.duplicated_type = 1;
      elseif gene <= organism.rs + organism.ks
        organism.duplicated_type = 2;
      elseif gene <= organism.rs + organism.ks + organism.ps
        organism.duplicated_type = 3;
      elseif gene <= organism.rs + organism.ks + organism.ps + organism.es
        organism.duplicated_type = 4;
      end
      
      % Useful
      n = organism.size;
      a = organism.k;
      
      % Matrix manipulation, grab the "four corners" of the interaction
      a1 = a(1:(gene-1), 1:(gene-1));
      a2 = a(1:(gene-1), (gene+1):n);
      a3 = a((gene+1):n, 1:(gene-1));
      a4 = a((gene+1):n, (gene+1):n);
      
      % Create a new vector for row i, with an extra 0 in it
      rowi = a(gene, 1:n);
      rowi = [ rowi(1:gene), 0, rowi((gene+1):n) ];
      
      % Get the column and split it into two parts
      coli = a(1:n, gene);
      colitop = coli(1:(gene-1));
      colibot = coli((gene+1):n);
      
      % Patch this all together
      B = [ a1 colitop colitop a2; 
            rowi; rowi;
            a3 colibot colibot a4];
      organism.k = B;
      
      % Now patch the other (not k) stuff together
      organism.a = [ organism.a(1:gene), 0, organism.a((gene+1):n) ];
      organism.d = [ organism.d(1:gene), 0, organism.d((gene+1):n) ];
      organism.alpha = [ organism.alpha(1:gene), 0, organism.alpha((gene+1):n) ];
      organism.protein_categorisations = [ organism.protein_categorisations(1:gene), ...
                                           organism.protein_categorisations(gene:n) ];
      
      % Update size
      organism.size = organism.size + 1;
      % TODO: update the other rs,ks,ps,es stores
    end

    % Mutates a gene (that has been duplicated
    function organism = mutateduplicate(organism, settings)
      % Convenience
      rs = organism.rs;
      ks = organism.ks;
      ps = organism.ps;
      es = organism.es;
      gene = organism.duplicated;
      
      % Loopy loop
      for i=1:organism.size
        % DO NOT ACT ON SELF
        if gene ~= i
        
          % manipulate the row
          if strcmp(settings.act_on, 'row') || strcmp(settings.act_on, 'both')
            % Regenerate as required
            if rand() < settings.severity
              if gene <= organism.rs
                organism.k(gene,i) = rand();
              elseif gene <= organism.rs + organism.ks
                organism.k(gene,i) = rand();
              elseif gene <= organism.rs + organism.ks + organism.ps
                organism.k(gene,i) = - rand();
              elseif gene <= organism.rs + organism.ks + organism.ps + organism.es
                organism.k(gene,i) = 2*rand() - 1;
              end

            % Perturb the data if we don't regenerate it
            else
              organism.k(gene,i) = ...
                organism.k(gene,i) + (organism.k(gene,i) ~= 0) * settings.severity * (1 - 2*rand());
            end
          end  % Stop acting on rows

          % maniupate the column
          if strcmp(settings.act_on, 'col') || strcmp(settings.act_on, 'both')
            % Regenerate as required
            if rand() < settings.severity
              if i <= organism.rs
                organism.k(i,gene) = rand();
              elseif i <= organism.rs + organism.ks
                organism.k(i,gene) = rand();
              elseif i <= organism.rs + organism.ks + organism.ps
                organism.k(i,gene) = - rand();
              elseif i <= organism.rs + organism.ks + organism.ps + organism.es
                organism.k(i,gene) = 2*rand() - 1;
              end

            % Perturb the data if we don't regenerate it
            else
              organism.k(i,gene) = ...
                organism.k(i,gene) + (organism.k(gene,i) ~= 0) * settings.severity * (1 - 2*rand());
            end
          end % stop acting on the cols
        
        end % end checking i~=gene
      end % End loop
    end % end function definition
    
    % Returns the ODEs that have to be integrated to represent network
    function ydot = odes(self, t, active, ligand)
        % Helpful
        k = self.k;
        a = self.a;
        d = self.d;
        alpha = self.alpha;
        inactive = 1 - active;
        
        % Actually set up the ODES
        ydot = zeros(self.size, 1);
        for i=1:self.size
            % Relaxation
            ydot(i) = ydot(i) + (a(i) * inactive(i));
            ydot(i) = ydot(i) + (d(i) * active(i));
            
            % Ligand
            if alpha(i) > 0
                ydot(i) = ydot(i) + (inactive(i) * ligand * alpha(i));
            elseif alpha(i) < 0
                ydot(i) = ydot(i) + (active(i) * ligand * alpha(i));
            end
            
            % Protein INTERactions
            for j = 1:self.size
                if k(j,i) > 0
                    ydot(i) = ydot(i) + (inactive(i) * k(j,i) * active(j));
                elseif k(j,i) < 0
                    ydot(i) = ydot(i) + (active(i) * k(j,i) * active(j));
                end
            end
        end
    end % End odes()
    
    % Returns teh jacobian matrix for the system
    function jac = jacobian(self, t, active, ligand)
    end
  end
  
end

