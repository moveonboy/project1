function isstructcontained = bracketContain(obj,obj2)
%BRACKETCONTAIN: To check if glycan struct contain another glycan structure
%  with bracket. 
% 
%  Syntax:  
%      issubset=glycanstructobj1.bracketContain(glycanstructobj2)
%      issubset=bracketContain(glycanstructobj1,glycanstructobj2)
%  
%  Input: 
%      glycanstructobj1: GlycanStruct object 1
%      glycanstructobj2: GlycanStruct object 2
%   
%  Output:
%      issubset indicating GlycanStruct object 2 is a subset of
%      GlycanStruct object 1
%
%  Example:
%      
%
% See also compContain

%Author: Yusen Zhou and Gang Liu
%Date: 11/19/2014

ismonocompsame=obj.compContain(obj2);
if(ismonocompsame==0)
    isstructcontained = 0;
    return    
end

nonbracketcomp    = isempty(obj.bracket)  && isempty(obj2.bracket);
obj1bracketcomp   = ~isempty(obj.bracket) && isempty(obj2.bracket);
obj2bracketcomp   = isempty(obj.bracket)  && ~isempty(obj2.bracket);
obj12bracketcomp = ~isempty(obj.bracket) && ~isempty(obj2.bracket);

isstructcontained = 1;
if(nonbracketcomp)
   isstructcontained = obj.contains(obj2);
   return
end

if(obj1bracketcomp)
    objtocompare=obj.clone;
    objtocompare.removeBracket;
    if(~objtocompare.contains(obj2))
      isstructcontained = 0;               
    end
    return
end

if(obj2bracketcomp)
    obj2tocompare = obj2.clone;
    obj2tocompare.removeBracket;   
    if(~obj.contains(obj2tocompare))
      isstructcontained = 0;
    end
    return
end

if(obj12bracketcomp)    
    if(~obj.glycanjava.getRoot.subtreeContains(...
            obj2.glycanjava.getRoot))
        isstructcontained = 0;
    end    
end

end

