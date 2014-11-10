function Buildingblocks = getBlocks(glycanspecies)
%isSubsetStructure returns true if one glyanstruct(obj) contains another
%  structure(obj2) and the number of obj2 structure that obj structure
%  contains.
% 
% isSubsetStructure(obj,obj2) compares the two structures, and check if
%    obj2 is fragmentation of obj.
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 9/23/14

load('GlyBlocksDB.mat')
glycanStruct   = glycanspecies.glycanStruct;
Buildingblocks = '';
counter        = 0;
Mapkeys        = keys(blockdatabase);

for i = 1 : length(Mapkeys)
    ithkey = Mapkeys{i};
    ithblockstruct = blockdatabase{ithkey};
    if(~isequal(ithkey,'LeX'))&&(~isequal(ithkey,'LacNAcI'))&&(~isequal(ithkey,'LacNAcII'))
        issubsetstruct = glycanStruct.contains(ithblockstruct);
        isexist        = checkexistance(ithblockstruct,Buildingblocks);
        if(~isexist)&&(issubsetstruct)
            Buildingblocks{counter+1} = ithblockstruct;
            counter = counter+1;
        end
    elseif(isequal(ithkey,'LeX'))
        [~,num1] = isSubsetStructure(glycanStruct,blockdatabase{'SLeX'});
        [~,num2] = isSubsetStructure(glycanStruct,blockdatabase{'a2_6SSLeX'});
        [~,num3] = isSubsetStructure(glycanStruct,blockdatabase{ithkey});
        isexist                = checkexistance(ithblockstruct,Buildingblocks);
        if(num3>num2+num1)&&(~isexist)
            Buildingblocks{counter+1} = ithblockstruct;
            counter = counter+1;
        end
    elseif(isequal(ithkey,'LacNAcI'))
        [~,num1] = isSubsetStructure(glycanStruct,blockdatabase{'SLacNAcI'});
        [~,num2] = isSubsetStructure(glycanStruct,blockdatabase{'a2_6SLacNAcI'});
        [~,num3] = isSubsetStructure(glycanStruct,blockdatabase{ithkey});
        isexist                = checkexistance(ithblockstruct,Buildingblocks);
        if(num3>num2+num1)&&(~isexist)
            Buildingblocks{counter+1} = ithblockstruct;
            counter = counter+1;
        end
    elseif(isequal(ithkey,'LacNAcII'))
        [~,num1] = isSubsetStructure(glycanStruct,blockdatabase{'SLacNAcII'});
        [~,num2] = isSubsetStructure(glycanStruct,blockdatabase{'a2_6SLacNAcII'});
        [~,num3] = isSubsetStructure(glycanStruct,blockdatabase{ithkey});
        isexist                = checkexistance(ithblockstruct,Buildingblocks);
        if(num3>num2+num1)&&(~isexist)
            Buildingblocks{counter+1} = ithblockstruct;
            counter = counter+1;
        end
    end
end
end

function isexist = checkexistance(Blocks,Buildingblocks)
isexist = 0;
for i = 1 : length(Buildingblocks)
    ithBlocks = Buildingblocks{i};
    if(isequal(ithBlocks.name,Blocks.name))
        isexist = 1;
        break
    end
end
end