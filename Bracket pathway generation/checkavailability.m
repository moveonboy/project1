function isProdValid = checkavailability(nreResidues,bracketresidue,enzObj,enzObjArray)
isProdValid = 1;
numrequiredresidue = 0;
for i = 1 : length(nreResidues)
    ithresidue = nreResidues{1,i};
    if(isequal(ithresidue.residueType.name,enzObj.resfuncgroup.name))
        numrequiredresidue = numrequiredresidue+1;
    end
end

acceptorname = '';
counter      = 0;
for i = 1 : length(enzObjArray)
    ithenzObj = enzObjArray.get(i);
    if(isequal(ithenzObj.resAtt2FG.name,enzObj.resfuncgroup.name))
        acceptorname{counter+1,1} = ithenzObj.resfuncgroup.name;
        counter = counter + 1;
    end
end

numacceptdresidue = 0;
for i = 1 : length(bracketresidue)
    ithresidue = bracketresidue{1,i};
    if(isequal(ithresidue.residueType.name,'#bracket'))
        ithchildren = ithresidue.getChildren;
        for j = 1 : length(ithchildren)
            jthchild = ithchildren(j);
            if(isequal(jthchild.residueType.name,enzObj.resfuncgroup.name))
                continue
            end
            
            for jj = 1 : length(acceptorname)
                jjthacceptorname = acceptorname{jj,1};
                if(isequal(jthchild.residueType.name,jjthacceptorname))
                    numacceptdresidue = numacceptdresidue+1;
                end
            end
        end
    end
end

if(numacceptdresidue>=numrequiredresidue)
    isProdValid = 0;
end
end