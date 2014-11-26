classdef Rxn  < handle
    %Rxn class representing a chemical reaction in the glycosylation
    % reaction network
    %
    % A Rxn object is a generic representation of a glycosylation reaction.
    %  It consists of an enzyme, a reactant and a product.
    %
    % Rxn properties:
    %  reac               -  the reactant in the reaction, i.e, acceptor
    %  prod               -  the product in the reaction
    %  donor              -  the glycosyl donor
    %  enz                -  the catalytic enzyme
    %  name               -  the name of the reaction
    %  mark               -  if the reaction is visited
    %  markState          -  if visited, the status of the mark
    %  cost               -  a cost to form the product from the reactant
    %
    % Rxn methods:
    %  RXN                -  create a Rxn object
    %  getReactant        -  retrieve the "reactant" property
    %  getProduct         -  retrieve the "product" property
    %  getDonor           -  retrieve the "donor" property
    %  getEnzyme          -  retrieve the "enzyme" property
    %  getName            -  retrieve the "name" property
    %  getCost            -  retrieve the "cost" property
    %  setReactant        -  set the "reactant" property
    %  setProduct         -  set the "product" property
    %  setDonor           -  set the "donor" property
    %  setEnzyme          -  set the "enzyme" property
    %  setName            -  set the "name" property
    %  setCost            -  set the "cost" property
    %  setMark            -  set the "mark" property as true
    %  clearMark          -  set the "mark" property as false
    %  isMarked           -  return the "mark" property
    %  clone              - copy a Rxn object
    %
    % See also Pathway,GlycanSpecies,Compt,Enz.
    
    % Author: Gang Liu
    % Date Lastly Updated 9/20/13
    
    properties
       rxntype;
    end   
    
    properties
        %REAC a reactant in the glycosylation reaction
        % REAC property is a GlycanSpecies object
        %
        %
        %See also Rxn
        reac
    end
    
    properties
        %PROD a product in the glycosylaiton reaction
        % PROD property is a GlycanSpecies object
        %
        %See also Rxn
        prod
    end
    
    properties
        %DONOR a glycosyl donor in the reaction
        % DONOR property is a GlycanSpecies object
        %
        %See also Rxn
        donor
    end
    
    properties
        %ENZ a glycosyltransferase in the reaction
        % ENZ property is a ENZ object
        %
        %See also Rxn
        enz
    end
    
    properties
        %COST a cost to form the product from the reactant
        % COST property is a numeric scalar.
        %
        %See also Rxn
        cost   % for graph structure purpose
    end
    
    properties
        %MARK the status if the reaction is visited
        % MARK property is a logical scalar
        %
        %See also Rxn
        mark  % for graph structure purpose
    end
    
    properties
        %MARKSTATE the state of the mark if visisted
        % MARKSTATE property is a character array
        %
        %See also Rxn
        markState  % for graph structure purpose
    end
    
    properties
        %NAME name of the reaction
        % NAME property is a character arrray
        %
        %See also Rxn
        name
    end
    
    properties
        %ID name of the reaction
        % ID property is a character arrray
        %
        %See also Rxn
        id
    end
    
    
    properties
        % RXNKINETICS reaction kinetics
        %   RXNKINETICS property is a MATLAB structure with 3 fields
        %   including "mathexpr','speciesnames','parameter'
        %
        % See also EnzKinetics.
        rxnkinetics
    end
    
    methods
        function setrxnkinetics(obj)
           if (~isempty(obj.enz))
                if(isa(obj.enz.enzkinetics,...
                        'MMenKinetics'))
                    obj.rxnkinetics = MMenRxnKinetics(...
                        obj.enz.enzkinetics,obj.reac);
                elseif(isa(obj.enz.enzkinetics,...
                        'BiBiKinetics'))
                    obj.rxnkinetics = BiBiRxnKinetics(...
                        obj.enz.enzkinetics,obj.reac,obj.prod);
                elseif(isa(obj.enz.enzkinetics,...
                        'RevMMenKinetics'))
                    obj.rxnkinetics = MMenRxnKinetics(...
                        obj.enz.enzkinetics,obj.reac,obj.prod);
                else
                
                end
           else % transport kinetics
               try 
                if(isempty(obj.reac))  % flux in reaction
                    qin             = obj.prod.compartment.reactortype.qflowin;
                    if(obj.prod.compartment.reactortype.speciesidinletdb.isKey(obj.prod.id))
                       initconc = obj.prod.compartment.reactortype.speciesidinletdb(obj.prod.id);
                    else
                       error('MATLAB:GNAT:ERRORSPECIESID','SPECIES IS NOT FOUND IN THE PATHWAY');
                       % initconc = 0;
                    end
                    
                    obj.rxnkinetics = TransportKinetics(TransportKinetics.zeroorder,...
                            qin,obj.prod.id,initconc);        
                else
                    obj.rxnkinetics = TransportKinetics(TransportKinetics.firstorder,...
                       obj.reac.compartment.reactortype.ktransport,obj.reac.id);                            
                end
               catch err
                   disp('Error in Rxn');
               end
           end
        end
    end
    
    methods (Static)
        function rxnobj=loadmat(matfilename)
            rxnstruct=load(matfilename);
            p = fieldnames(rxnstruct);
            if(length(p)~=1)
                error('MATLAB:GNAT:WRONGINPUT','wrong variable stored');
            end
            rxnobj = rxnstruct.(p{1});
            if(isa(rxnobj,'Rxn'))
                rxnobj.resetjava;
            else
                error('MATLAB:GNAT:WRONGINPUT','wrong variable stored');
            end
        end
    end
    
    methods
        function resetjava(obj)
            
            if ~(isempty(obj.getReactant))
                obj.getReactant.glycanStruct.resetjava;
            end
            
            if ~(isempty(obj.getProduct))
                obj.getProduct.glycanStruct.resetjava;
            end
            
        end
    end
    
    methods
        function obj = Rxn(varargin)
            %Rxn create a Rxn object.
            %
            % GR  = Rxn(REACTANT,PRODUCT) creates a Rxn object with its
            %  acceptor REACTANT and the product PRODUCT.
            %
            % GR  = Rxn(REACTANT,PRODUCT,ENZYME) creates a Rxn object with its
            %  acceptor REACTANT, the product PRODUCT and its catalytic
            %  enzyme ENZYME.
            %
            % GR  = Rxn(REACTANT,PRODUCT,ENZYME,COST) specifies the cost of PRODUCT formation
            %  from REACTANT using Enzyme.
            %
            % See also GlycanStruct,Compt,Pathway,GlycanNetModel.
            
            if (nargin == 0)
                obj.name    = 'empty';
            elseif(nargin==1)
                errorReport(mfilename,'IncorrectInputNumber');
            elseif(nargin==2)
                if(isa(varargin{1},'GlycanSpecies') && isa( varargin{2},'GlycanSpecies'))
                    obj.reac  =  varargin{1};
                    obj.prod =  varargin{2};
                    obj.cost =  0;
                elseif(isempty(varargin{1}) && isa( varargin{2},'GlycanSpecies'))
                    obj.reac  = [];
                    obj.prod =  varargin{2};
                    obj.cost =  0;
                elseif(isempty(varargin{2}) && isa( varargin{1},'GlycanSpecies'))
                    obj.reac  = varargin{1};
                    obj.prod =  [];
                    obj.cost =  0;
                else
                    error('MATLAB:GNAT:WRONGTINPUTTYPE','Incorrect Input Type');
                end
            elseif(nargin==3)
                if(isa(varargin{2},'GlycanSpecies')...
                && isa( varargin{2},'GlycanSpecies')...
                        &&isa(varargin{3},'Enz'))
                    obj.reac   = varargin{1};
                    obj.prod  =  varargin{2};
                    obj.enz    =  varargin{3};
                    obj.cost  =  0;
                elseif((isa(varargin{1},'GlycanSpecies') &&...
                        isa(varargin{2},'GlycanSpecies') )...
                        &&isa(varargin{3},'char'))
                    obj.reac      =  varargin{1};
                    obj.prod      =  varargin{2};
                    obj.name      =  varargin{3};
                    obj.cost      =  0;
                elseif(isempty(varargin{1})&& isa( varargin{2},'GlycanSpecies')...
                       &&isempty(varargin{3}))
                    obj.reac      = [];
                    obj.enz       = [];
                    obj.prod      = varargin{2};
                    obj.cost      =  0;
                elseif(isempty(varargin{2})&& isa( varargin{1},'GlycanSpecies')...
                       &&isempty(varargin{3}))
                    obj.prod      = [];
                    obj.enz       = [];
                    obj.reac      =  varargin{1};    
                    obj.cost      =  0;
                elseif isa(varargin{1},'GlycanSpecies')...
                        && isa(varargin{1},'GlycanSpecies')...
                        && isempty(varargin{3})
                    obj.prod      =  varargin{2};
                    obj.enz       = [];
                    obj.reac      =  varargin{1};    
                    obj.cost      =  0;
                elseif  isempty(varargin{1})...
                        && isa(varargin{2},'GlycanSpecies')...
                        && ischar(varargin{3})
                    obj.prod      =  varargin{2};
                    obj.enz       =  [];
                    obj.reac      =  [];    
                    obj.cost      =  0;
                    obj.name      = varargin{3};
               elseif  isempty(varargin{2})...
                        && isa(varargin{1},'GlycanSpecies')...
                        && ischar(varargin{3})
                    obj.prod      =  varargin{1};
                    obj.enz       =  [];
                    obj.reac      =  [];    
                    obj.cost      =  0;  
                    obj.name      = varargin{3};
                else
                    fprintf(1,'%i %i %i ',isa(varargin{1},'GlycanSpecies'),...
                        isa(varargin{2},'GlycanSpecies'),isa(varargin{3},'char'));
                    error('MATLAB:GNAT:WRONGTINPUTTYPE','Incorrect Input Type');
                end
            elseif(nargin==4)
                if(isa( varargin{1},'GlycanSpecies') && isa( varargin{2},'GlycanSpecies')...
                        &&isa( varargin{3},'Enz') &&isa( varargin{4},'double'))
                    obj.reac  =  varargin{1};
                    obj.prod  =  varargin{2};
                    obj.enz    =  varargin{3};
                    obj.cost   =  varargin{4};               
                else
                    error('MATLAB:GNAT:INCORRECTINPUTTYPE','Incorrect Input Type');
                end
            end
        end % end constructor
    end
    
    methods
        function prodToForm = getProduct(obj)
            % getProduct get the prodcut of the reaction
            %
            % See also Rxn
            prodToForm = obj.prod;
        end
        
        function setProduct(obj,prodToForm)
            % setProduct set the prodcut of the reaction
            %
            % See also Rxn
            obj.prod = prodToForm ;
        end
        
        function reacToAct = getReactant(obj)
            % getReactant get the reactant of the reaction
            %
            % See also Rxn
            reacToAct = obj.reac;
        end
        
        function setReactant(obj,reacToAct )
            % setReactant set the reactant of the reaction
            %
            % See also Rxn
            obj.reac = reacToAct;
        end
        
        function enzToFind = getEnzyme(obj)
            % getEnzyme get the enzyme of the reaction
            %
            % See also Rxn
            enzToFind = obj.enz;
        end
        
        function setEnzyme(obj,enzToSet)
            % setEnzyme set the enzyme of the reaction
            %
            % See also Rxn
            obj.enz = enzToSet;
        end
        
        function donor = getDonor(obj)
            % getDonor get the donor of the reaction
            %
            % See also Rxn
            donor = obj.donor;
        end
        
        function setDonor(obj,donor)
            % setDonor set the donor of the reaction
            %
            % See also Rxn
            obj.donor = donor;
        end
        
        function name = getName(obj)
            % getName get the name of the reaction
            %
            % See also Rxn
            name = obj.name;
        end
        
        function setName(obj,name)
            % setName set the name of the reaction
            %
            % See also Rxn
            obj.name = name;
        end
        
        function costOfRxn = getCost(obj)
            % getCost get the cost of the reaction
            %
            % See also Rxn
            costOfRxn = obj.cost;
        end
        
        function  setCost(obj,costOfRxn)
            % setCost set the cost of the reaction
            %
            % See also Rxn
            obj.cost = costOfRxn;
        end
        
        function setMark(obj)
            % setMark set the reaction "visited"
            %
            % See also Rxn
            obj.mark = true;
        end
        
        function clearMark(obj)
            % clearMark set the reaction "NOT visited"
            %
            % See also Rxn
            obj.mark = false;
        end
        
        function markStatus = isMarked(obj)
            % isMarked return true if the reaction is "visted"
            %
            % See also Rxn
            markStatus=obj.mark;
        end
    end
    
    methods
        toNameStr = toString(obj);
    end
    
    methods
        sbmlexpr  = tosbmlrxn(obj);
    end
    
    methods
        function new = clone(obj)
            %CLONE copy a Rxn object
            %  clone(RXNobj) creates a new RXN object
            %
            %  See also RXN
            
            % Instantiate new object of the same class.
            new = feval(class(obj));
            
            % Copy all non-hidden properties.
            p = properties(obj);
            for i = 1:length(p)
                new.(p{i}) = obj.(p{i});
            end
        end
    end
    
end


