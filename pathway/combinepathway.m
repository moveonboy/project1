function nlinkedpath = combinepathway(varargin)
if(nargin<2)
    error('MATLAB:GNAT:ERRORNONCOMPLEX','WRONG NUMBER OF INPUT');
end

nlinkedpath = Pathway;
filename = [varargin{1} '.mat'];
pathway1 = load(filename);

% Combine the Rxns
totallistOfRxns = pathway1.nlinkedpath.theRxns;
for i = 2 : length(varargin)
    filename = [varargin{i} '.mat'];
    ithpathway = load(filename);
    ithlistOfRxn = ithpathway.nlinkedpath.theRxns;
    for j = 1 : length(ithlistOfRxn)
        jthrxn   = ithlistOfRxn.get(j);
        isnewrxn = 1;
        for jj = 1 :length(totallistOfRxns)
            jjthrxn   = totallistOfRxns.get(jj);
            isSameRxn = comparerxn(jthrxn,jjthrxn);
            if(isSameRxn);
                isnewrxn = 0;
                break;
            end
        end
        if(isnewrxn)
           totallistOfRxns.add(jthrxn)
        end
    end
end
nlinkedpath.theRxns = totallistOfRxns;

% Combine the species
totallistOfSpecies = CellArrayList;
for i = 1 : length(totallistOfRxns)
    ithrxn     = totallistOfRxns.get(i);
    ithrxnreac = ithrxn.reac;
    ithrxnprod = ithrxn.prod;
    isnewreac = 1;
    isnewprod = 1;
    isbreak   = 0;
    for j = 1 :length(totallistOfSpecies)
        jthspecies = totallistOfSpecies.get(j);
        if(isequal(ithrxnreac.glycanStruct.name,jthspecies.glycanStruct.name))
            nlinkedpath.theRxns.get(i).reac = jthspecies;
            isnewreac = 0;
            isbreak   = isbreak+1;
        elseif(isequal(ithrxnprod.glycanStruct.name,jthspecies.glycanStruct.name))
            nlinkedpath.theRxns.get(i).prod = jthspecies;
            isnewprod = 0;
            isbreak   = isbreak+1;
        end
        if(isbreak==2)
            break
        end
    end
    if(isnewreac)
        totallistOfSpecies.add(ithrxnreac)
    end
    if(isnewprod)
        totallistOfSpecies.add(ithrxnprod)
    end
end
nlinkedpath.theSpecies = totallistOfSpecies;

% Combine the enzymes
totallistOfEnzs = pathway1.nlinkedpath.theEnzs;
for i = 2 : length(varargin)
    filename = [varargin{i} '.mat'];
    ithpathway = load(filename);
    ithlistOfEnzs = ithpathway.nlinkedpath.theEnzs;
    for j = 1 : length(ithlistOfEnzs)
        jthenz   = ithlistOfEnzs.get(j);
        isnewenz = 1;
        for jj = 1 :length(totallistOfEnzs)
            jjthenz = totallistOfEnzs.get(jj);
            if(isequal(jthenz.name,jjthenz.name));
                isnewenz = 0;
                break;
            end
        end
        if(isnewenz)
           totallistOfEnzs.add(jthenz)
        end
    end
end
nlinkedpath.theEnzs = totallistOfEnzs;

end

function isSameRxn = comparerxn(rxn1,rxn2)
isreac = isequal(rxn1.reac.glycanStruct.name,rxn2.reac.glycanStruct.name);
isprod = isequal(rxn1.prod.glycanStruct.name,rxn2.prod.glycanStruct.name);
isenz  = isequal(rxn1.enz.name,rxn2.enz.name);
if(isreac)&&(isprod)&&(isenz)
    isSameRxn = 1;
else
    isSameRxn = 0;
end
end