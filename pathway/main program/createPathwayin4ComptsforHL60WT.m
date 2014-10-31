function createPathwayin4ComptsforHL60WT()
%% load pathway from *.mat file
matfilename         = 'HL60WTN_GlycanPathwayReduced.mat';
fprintf(1,'%s\n','load pathway file');
HL60wildtypepathway  = Pathway.loadmat(matfilename);

%% set up four compts
theGolgi4Compts = createGolgi4Compts();

%% Set up species and rxns in each compartment
fprintf(1,'%s\n','clone pathway to 4 compartments');
HL60wildtypepathway.cloneNetworkinAllCompts(theGolgi4Compts);

% set up enzyme name
HL60wildtypepathway = setenzname(HL60wildtypepathway);

%% Set up species ID/Init Conc for each species, reaction ID for each reaction
fprintf(1,'%s\n','set up IDs in the pathway');
HL60wildtypepathway.setSpeciesIDIndex;
HL60wildtypepathway.setRxnsIDIndex;
HL60wildtypepathway.setSpeciesInitConcZero;
HL60wildtypepathway.setEnzIDIndex;

%% Set up rxn kinetics in each compartment
fprintf(1,'%s\n','set up Enzyme kinetics in the pathway');
BiBikineticsdb                = CreateBiBiSeqDBforHL60();  % BiBikinetic database
sugarnucldb                   = setSugarDB();
substrComptOption.substrCompt = 1;

nglycanenzdistrib      = setEnzDist();
HL60wildtypepathway     = pathwaySetKinetics(HL60wildtypepathway,...
    BiBikineticsdb,substrComptOption,sugarnucldb,nglycanenzdistrib);

%% set up compartment reactor model
fprintf(1,'%s\n','set up Reactor model');

speciesInletConc      = setSpeciesinletconcForHL60(HL60wildtypepathway);
reactorpara           = struct('qflow',1,'kt',0.18*60,...
    'inletConc',speciesInletConc);
glycanNetReactorModel = pathway2ReactorModel(HL60wildtypepathway,...
    TypeStringConst.CSTR,reactorpara);
glycanNetViewer(glycanNetReactorModel);

% Output SBML file
fprintf(1,'%s\n','output SBML file');
testsbmlstruct  = glycanNetReactorModel.toSBMLStruct(true);
glycanNetReactorModel.glycanNet_sbmlmodel = testsbmlstruct;
save('HL60Model.mat','glycanNetReactorModel');
OutputSBML(testsbmlstruct);
end

function nglycanenzdistrib=setEnzDist()
nglycanenzdistrib = containers.Map;
enznames = {'mgat1','mgat2','mgat4','mgat3',...
    'mgat5','b4GalI','manii','mani','Futa23',...
    'Fut8','iGnt','STGalTs'};
for i = 1 : length(enznames)
    if(~strcmpi(enznames{i},'b4GalI'))&&(~strcmpi(enznames{i},'iGnt'))&&...
            (~strcmpi(enznames{i},'STGalTs'))&&(~strcmpi(enznames{i},'Fut8'))&&...
            (~strcmpi(enznames{i},'Futa23'))
        nglycanenzdistrib(enznames{i})=struct('cis_golgi',0.15,...
            'media_golgi',0.45,'trans_golgi',0.30,'trans_golgi_network',0.30);
    elseif(strcmpi(enznames{i},'b4GalI'))||(strcmpi(enznames{i},'STGalTs'))
        nglycanenzdistrib(enznames{i})=struct('cis_golgi',1,...
            'media_golgi',0.05,'trans_golgi',0.20,'trans_golgi_network',0.70);
    elseif(strcmpi(enznames{i},'iGnt'))
        nglycanenzdistrib(enznames{i})=struct('cis_golgi',1,...
            'media_golgi',0.60,'trans_golgi',0.30,'trans_golgi_network',0.10);
    elseif(strcmpi(enznames{i},'Fut8'))||(strcmpi(enznames{i},'Futa23'))
        nglycanenzdistrib(enznames{i})=struct('cis_golgi',0.1,...
            'media_golgi',0.15,'trans_golgi',0.45,'trans_golgi_network',0.30);
    end
end
end

function testpathway=setenzname(testpathway)
enzdbmatfilename  =  'HL60enzDB.mat';
enzdb             =  enzdbmatLoad(enzdbmatfilename);
mgat1             = enzdb('mgat1');
mgat2             = enzdb('mgat2');
mgat3             = enzdb('mgat3');
mgat4             = enzdb('mgat4');
mgat5             = enzdb('mgat5');
b4GalI            = enzdb('b4GalI');
manii             = enzdb('manii');
Fut8              = enzdb('Fut8');
iGnt              = enzdb('iGnt');
mani              = enzdb('mani');
Futa23            = enzdb('Futa23');
STGalTs          = enzdb('STGalTs');
for i = 1 : length(testpathway.theEnzs)
    if(isequal(testpathway.theEnzs.get(i).ecno,mgat1.ecno))
        testpathway.theEnzs.get(i).name = 'mgat1';
    elseif(isequal(testpathway.theEnzs.get(i).ecno,mgat2.ecno))
        testpathway.theEnzs.get(i).name = 'mgat2';
    elseif(isequal(testpathway.theEnzs.get(i).ecno,mgat3.ecno))
        testpathway.theEnzs.get(i).name = 'mgat3';
    elseif(isequal(testpathway.theEnzs.get(i).ecno,mgat4.ecno))
        testpathway.theEnzs.get(i).name = 'mgat4';
    elseif(isequal(testpathway.theEnzs.get(i).ecno,mgat5.ecno))
        testpathway.theEnzs.get(i).name = 'mgat5';
    elseif(isequal(testpathway.theEnzs.get(i).ecno,b4GalI.ecno))
        testpathway.theEnzs.get(i).name = 'b4GalI';
    elseif(isequal(testpathway.theEnzs.get(i).ecno,manii.ecno))
        testpathway.theEnzs.get(i).name = 'manii';
    elseif(isequal(testpathway.theEnzs.get(i).ecno,Fut8.ecno))
        testpathway.theEnzs.get(i).name = 'Fut8';
    elseif(isequal(testpathway.theEnzs.get(i).ecno,iGnt.ecno))
        testpathway.theEnzs.get(i).name = 'iGnt';
    elseif(isequal(testpathway.theEnzs.get(i).ecno,mani.ecno))
        testpathway.theEnzs.get(i).name = 'mani';
    elseif(isequal(testpathway.theEnzs.get(i).ecno,Futa23.ecno))
        testpathway.theEnzs.get(i).name = 'Futa23';
    elseif(isequal(testpathway.theEnzs.get(i).ecno,STGalTs.ecno))
        testpathway.theEnzs.get(i).name = 'STGalTs';
    end
end

end

function speciesInletConc = setSpeciesinletconcForHL60(testpathway)
speciesInletConc = CellArrayList;
firstcomptname = testpathway.findFirstCompt.name;
species       = GlycanSpecies(glycanMLread('M5.glycoct_xml'));
inletspecies  = testpathway.findSpeciesByStruct(...
    species,firstcomptname);
speciesInletConc.add(struct('species',inletspecies.id ,'inletconc',1));
% speciesInletConc.add(struct('species','G2' ,'inletconc',0.5));
for i = 1 :  length(testpathway.theSpecies)
    speciesid = ['G' num2str(i)];
    if(strcmpi(speciesid,inletspecies.id))
        continue
    else
        Inletstruct = struct('species',speciesid ,'inletconc',0);
        speciesInletConc.add(Inletstruct);
    end
end
end

function theCompt = createGolgi4Compts()
    theCompt   = CellArrayList;
    cisCompt   = Compt('cis_golgi');
    mediaCompt = Compt('media_golgi');
    transCompt = Compt('trans_golgi');
    tgnCompt   = Compt('trans_golgi_network');
    cisCompt.posteriorcompt   = mediaCompt;
    mediaCompt.priorcompt     = cisCompt;
    mediaCompt.posteriorcompt = transCompt;
    transCompt.priorcompt     = mediaCompt;
    transCompt.posteriorcompt = tgnCompt;
    tgnCompt.priorcompt       = transCompt;
    theCompt.add(cisCompt);
    theCompt.add(mediaCompt);
    theCompt.add(transCompt);
    theCompt.add(tgnCompt);
end

function sugarnucldb = setSugarDB()
    sugarnucldb = containers.Map;
    sugarnucldb('UDP_GlcNAc') = 9.6*1e3;
    sugarnucldb('UDP_Gal')    = 4*1e3;
    sugarnucldb('GDP_Fuc')    = 4*1e3;
    sugarnucldb('CMP_NeuAc')  = 2.5*1e3;
    sugarnucldb('UDP')        = 0;
    sugarnucldb('GDP')        = 0;
    sugarnucldb('CMP')        = 0;
end