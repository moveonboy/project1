function isstructcontained = bracketContain2(obj,obj2)
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
% See also compContain

%Author: Yusen Zhou and Gang Liu
%Date: 11/19/2014

ismonocompsame=compcontain(obj,obj2);
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
    obj1list = bracket2exact(obj1);
    for i = 1 : length(obj1list)
       if(~obj1list{i}.glycanjava.root.contains(obj2.glycanjava.root))
           isstructcontained = 0;
           return
       end
    end
    return
end

if(obj2bracketcomp)    
    obj2list = bracket2exact(obj2);
    for i = 1 : length(obj2list)
       if(~obj1.glycanjava.root.contains(obj2list{i}.glycanjava.root))
           isstructcontained = 0;
           return
       end
    end
    return
end

if(obj12bracketcomp)
    obj1list = bracket2exact(obj1);
    obj2list = bracket2exact(obj2);
    for i = 1 : length(obj1list)
        for j = 1 : length(obj2list)
            if(~obj1list{i}.glycanjava.root.contains(obj2list{j}.glycanjava.root))
                isstructcontained = 0;
                return
            end
        end
    end    
end

end

