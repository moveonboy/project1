function GlyEnz = GlyEnzdatabase(MSglycandatafile,enzDB)
%GlyEnzdatabase returns the infromation of enzymes acted in the generation
% of glycans and store it in a database.
%
%
%GlyEnz = GlyEnzdatabase(MSglycandatafile,enzDB) using local glycan
% database and enzyme database as inputs to generate the database of
% combination of glycans and enzyme.
%
%
%
% Author: Yusen Zhou
% Date Lastly Updated: 10/1/14
GlyEnz = containers.Map;
load(MSglycandatafile);
load(enzDB);
for i = 1 : length(glycaninMS)
    ithglycan         = glycaninMS.get(i);
    ithglycanEnz      = createEnz(ithglycan.glycanStruct,enzDB);
    GlyEnz(ithglycan.glycanStruct.name) = ithglycanEnz;
end
end

function glycanEnz = createEnz(obj,enzDB)
glycanEnz   = '';
allResidues = getAllResidues(obj);
counter     = 0;
Mapkeys     = keys(enzDB);
for i = 1 : length(Mapkeys)
    ithkey = Mapkeys{i};
    ithenz = enzDB(ithkey);
    for j = 1 : length(allResidues)
        jthresidue    = allResidues{1,j};
        if(~isequal(jthresidue.residueType.name,ithenz.resfuncgroup.name))
            continue;
        end
        
        parentresidue = jthresidue.getParent;
        if(~isequal(parentresidue.residueType.name,ithenz.resAtt2FG.name))
            continue;
        end
        
        isresfuncvalid = 0;
        if(isequal(ithenz.linkFG.bond,jthresidue.linkageParent.bonds))&&...
                (isequal(ithenz.linkFG.anomer,jthresidue.anomer.symbol))    
            isresfuncvalid = 1;
        end
        
        isresAttvalid = 0;
        if(isequal(ithenz.linkresAtt2FG.bond,parentresidue.linkageParent.bonds))&&...
            (isequal(ithenz.linkresAtt2FG.anomer,parentresidue.anomer.symbol))
            isresAttvalid = 1;
        end
        
        isexist = checkexistance(ithenz,glycanEnz);
        
        if(~isexist)&& isresfuncvalid && isresAttvalid
           glycanEnz{counter+1,1} = ithenz;
           counter = counter+1;
        end
        
    end
end

end

function isexist = checkexistance(Enz,glycanEnz)
isexist = 0;
for i = 1 : length(glycanEnz)
    ithGlyEnz = glycanEnz{i};
    if(isequal(ithGlyEnz.name,Enz.name))
        isexist = 1;
        break
    end
end
end
