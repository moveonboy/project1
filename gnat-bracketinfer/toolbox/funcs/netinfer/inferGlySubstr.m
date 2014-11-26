function varargout=inferGlySubstr(prodObj,enzObj,varargin)
% inferGlySubstr infer the substrate based on the enzyme and product.
%
%   [substrSpecies] = inferGlySubstr(prodObj,enz) infers
%    the substrate if the enzyme acts on the substrate to form the product
%    . If no substrate is  found, substrSpecies is empty.
%    If more  than one substrates are found, they are stored in
%    a CellArrayList object substrSpecies.
%
%   [numSubstr,substrSpecies] = inferGlySubstr(prodObj,enz) returns the
%    number of substrates inferred.
%
%   [numSubstr,substrSpecies,rxns] = inferGlySubstr(prodObj,enz)
%   returns a list of reactions (a CellArrayList object rxns) if the enzyme acts
%   to form the product. If no substrate is found, rxns returns as an empty
%   CellArrayList object.
%
%  [numSubstr,substrSpecies,rxns,pathway] = inferGlySubstr(prodObj,enz)
%   returns the pathway  if the enzyme acts to form the product prodObj. If
%   no substrate  is found, pathway is set empty.
%
%      Example 1:
%            mani  =GHEnz.loadmat('mani.mat');
%            m8species  = GlycanSpecies(glycanMLread('M8.glycoct_xml'));
%            [nsubstr, m9species] = inferGlySubstr(m8species,mani);
%            options  = displayset('showmass',true,'showLinkage',true,...
%                             'showRedEnd',true);
%            for i = 1: nsubstr
%                glycanViewer(m9species.get(i).glycanStruct,options);
%            end
%
%      Example 2:
%             mani  =GHEnz.loadmat('mani.mat');
%            m8species  = GlycanSpecies(glycanMLread('M8.glycoct_xml'));
%            [nsubstr, m9species,m8rxns] = inferGlySubstr(m8species,mani);
%            for i = 1: nsubstr
%               glycanRxnViewer(m8rxns.get(i));
%            end
%
%      Example 3:
%            mani  =GHEnz.loadmat('mani.mat');
%            m8species  = GlycanSpecies(glycanMLread('M8.glycoct_xml'));
%            [nsubstr, m9species,m8rxns,m8pathway] = inferGlySubstr(m8species,mani);
%            glycanPathViewer(m8pathway);
%
% See also inferGlyProd,inferGlyRevrPath.

% Author: Gang Liu
% Date Last Updated: 11/18/14

if(~isa(prodObj,'GlycanSpecies'))
    error('MATLAB:GNAT:INCORRECINPUTTYPE','WRONG INPUT TYPE');
end

if(length(varargin)==1)
    if(strcmpi(varargin{1},'bracket'))
      usebracket=1;
    elseif(strcmpi(varargin{1},'nobracket'))
      usebracket=0;
    else
      error('MATLAB:GNAT:INCORRECTINPUTSTRING','WRONG INPUT STRING');  
    end
else
    error('MATLAB:GNAT:INCORRECTINPUTNUMBER','WRONG NUMBER OF INPUTS');  
end

if(isa(enzObj,'GTEnz'));
    substrSpecies = inferGTGlySubstr(prodObj,enzObj,usebracket);    
elseif(isa(enzObj,'GHEnz'))
    substrSpecies = inferGHGlySubstr(prodObj,enzObj,usebracket);    
else
    error('MATLAB:GNAT:NOTSUPPORTEDENZYME','NOT SUPPORTED ENZYME IN GNAT');
end

if(nargout>=3)
   rxns=CellArrayList;
   if(substrSpecies.length>0)
       for i = 1 : length(substrSpecies)
          rxns.add(Rxn(substrSpecies.get(i),prodObj,enzObj));
       end
   end 
    
   if(nargout==4)   
         path=Pathway;
        if(substrSpecies.length>0)
           path.addGlycans(substrSpecies);
           path.addGlycan(prodObj);
           path.addRxns(rxns);
        end
   end
end

if(nargout==1)
    varargout{1}=substrSpecies;
elseif(nargout==2)
    varargout{1}=length(substrSpecies);
    varargout{2}=substrSpecies;
elseif(nargout==3)
    varargout{1}=length(substrSpecies);
    varargout{2}=substrSpecies;
    varargout{3}=rxns;
elseif(nargout==4)
    varargout{1}=length(substrSpecies);
    varargout{2}=substrSpecies;
    varargout{3}=rxns;
    varargout{4}=path;
end

end