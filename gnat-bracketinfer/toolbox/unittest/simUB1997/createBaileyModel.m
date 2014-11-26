% test example
function createBaileyModel()
clc;clear

%% set four compt in UB1997
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

%% Load Enzyme Database
enzdbmatfilename   = 'glyenzDB.mat';
enzdb              = enzdbmatLoad(enzdbmatfilename);

%% Retrieve 33 Glycan Structures
glycanstructdb     = setGlyStruct();

%% Set Up 33 species and 33 rxns in each compartment
allspecies = CellArrayList;
allrxns    = CellArrayList;
for i = 1 : length(theCompt);
    [allspecies,allrxns] = setGlySpeciesRxninCompt(allspecies,allrxns,...
     glycanstructdb,theCompt.get(i),enzdb);
end

%% Set up species ID for each species
allspecies = setSpeciesIDAndInitConc(allspecies);
allspecies.get(1).initConc   = 0.0019;
allspecies.get(1).initAmount = 0.0019;
allspecies.get(2).initConc   = 0.0019;
allspecies.get(2).initAmount = 0.0019;

%% set up enzyme kinetics
allenzs = setupEnzKinetics(enzdb,glycanstructdb);

%% set up pathway
testpathway = Pathway;
testpathway.addGlycans(allspecies);
testpathway.addRxns(allrxns)
testpathway.addEnzs(allenzs);
testpathway.setRxnsIDIndex;
testpathway.setEnzIDIndex;
testpathway.setEnzRxnsKinetics;

%% Set up rxn kinetics in each compartment
theRxns             = testpathway.theRxns;
testpathway.theRxns = setupRxnKinetics(theRxns,enzdb);

%% set up compartment reactor
testpathway.setComptReactorType(TypeStringConst.CSTR);

% set qflow and k in each compts
qflow  = 1;
k      = 0.18*60;
testpathway.setComptsQflow(qflow);
testpathway.setComptskt(k);

%first compt
setSpeciesinletconcForUB1997(cisCompt);
testpathway.setTransporRxnKinetics;

modelname        = 'testSBML';
testmodel        = GlycanNetModel(theCompt,testpathway,modelname);
testsbmlstruct   = testmodel.toSBMLStruct;
glycanNetViewer(testmodel);
OutputSBML(testsbmlstruct);
end

function setSpeciesinletconcForUB1997(cisCompt)
speciesInletConc = CellArrayList;
speciesInletConc.add(struct('species','m9species_cis' ,'inletconc',1));
speciesInletConc.add(struct('species','m8species_cis' ,'inletconc',0));
speciesInletConc.add(struct('species','m7species_cis' ,'inletconc',0));
speciesInletConc.add(struct('species','m6species_cis' ,'inletconc',0));
speciesInletConc.add(struct('species','m5species_cis','inletconc',0));
speciesInletConc.add(struct('species','m5gnspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m4gnspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gnspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn2species_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn3species_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn3bspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn4species_cis','inletconc',0));
speciesInletConc.add(struct('species','m5gngspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m4gngspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gngspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn2gspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn3gspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn3bgspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn4bgspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m5gngnbspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m4gngnbspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gngnbspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn2gnbspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn3gnbspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn3bgnbspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn4gnbspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m5gngnbgspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m4gngnbgspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gngnbgspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn2gnbgspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn3gnbgspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn3bgnbgspecies_cis','inletconc',0));
speciesInletConc.add(struct('species','m3gn4gnbgspecies_cis','inletconc',0));
cisCompt.reactortype.setSpeciesInletDB(speciesInletConc)
end

function allglycans = setSpeciesIDAndInitConc(allglycans)
speciesname ={'m9species_cis','m8species_cis','m7species_cis','m6species_cis',...
    'm5species_cis','m5gnspecies_cis','m4gnspecies_cis','m3gnspecies_cis',...
    'm3gn2species_cis','m3gn3species_cis','m3gn3bspecies_cis','m3gn4species_cis',...
    'm5gngspecies_cis','m4gngspecies_cis','m3gngspecies_cis','m3gn2gspecies_cis',...
    'm3gn3gspecies_cis','m3gn3bgspecies_cis','m3gn4bgspecies_cis','m5gngnbspecies_cis',...
    'm4gngnbspecies_cis','m3gngnbspecies_cis','m3gn2gnbspecies_cis','m3gn3gnbspecies_cis',...
    'm3gn3bgnbspecies_cis','m3gn4gnbspecies_cis','m5gngnbgspecies_cis','m4gngnbgspecies_cis',...
    'm3gngnbgspecies_cis','m3gn2gnbgspecies_cis','m3gn3gnbgspecies_cis','m3gn3bgnbgspecies_cis',...
    'm3gn4gnbgspecies_cis',...
    'm9species_media','m8species_media','m7species_media','m6species_media'...
    'm5species_media','m5gnspecies_media','m4gnspecies_media','m3gnspecies_media',...
    'm3gn2species_media','m3gn3species_media','m3gn3bspecies_media','m3gn4species_media',...
    'm5gngspecies_media','m4gngspecies_media','m3gngspecies_media','m3gn2gspecies_media',...
    'm3gn3gspecies_media','m3gn3bgspecies_media','m3gn4bgspecies_media','m5gngnbspecies_media',...
    'm4gngnbspecies_media','m3gngnbspecies_media','m3gn2gnbspecies_media','m3gn3gnbspecies_media',...
    'm3gn3bgnbspecies_media','m3gn4gnbspecies_media','m5gngnbgspecies_media','m4gngnbgspecies_media',...
    'm3gngnbgspecies_media','m3gn2gnbgspecies_media','m3gn3gnbgspecies_media','m3gn3bgnbgspecies_media',...
    'm3gn4gnbgspecies_media',...
    'm9species_trans','m8species_trans','m7species_trans','m6species_trans',...
    'm5species_trans','m5gnspecies_trans','m4gnspecies_trans','m3gnspecies_trans',...
    'm3gn2species_trans','m3gn3species_trans','m3gn3bspecies_trans','m3gn4species_trans',...
    'm5gngspecies_trans','m4gngspecies_trans','m3gngspecies_trans','m3gn2gspecies_trans',...
    'm3gn3gspecies_trans','m3gn3bgspecies_trans','m3gn4bgspecies_trans','m5gngnbspecies_trans',...
    'm4gngnbspecies_trans','m3gngnbspecies_trans','m3gn2gnbspecies_trans','m3gn3gnbspecies_trans',...
    'm3gn3bgnbspecies_trans','m3gn4gnbspecies_trans','m5gngnbgspecies_trans','m4gngnbgspecies_trans',...
    'm3gngnbgspecies_trans','m3gn2gnbgspecies_trans','m3gn3gnbgspecies_trans','m3gn3bgnbgspecies_trans',...
    'm3gn4gnbgspecies_trans',...
    'm9species_tgn','m8species_tgn','m7species_tgn','m6species_tgn',...
    'm5species_tgn','m5gnspecies_tgn','m4gnspecies_tgn','m3gnspecies_tgn',...
    'm3gn2species_tgn','m3gn3species_tgn','m3gn3bspecies_tgn','m3gn4species_tgn',...
    'm5gngspecies_tgn','m4gngspecies_tgn','m3gngspecies_tgn','m3gn2gspecies_tgn',...
    'm3gn3gspecies_tgn','m3gn3bgspecies_tgn','m3gn4bgspecies_tgn','m5gngnbspecies_tgn',...
    'm4gngnbspecies_tgn','m3gngnbspecies_tgn','m3gn2gnbspecies_tgn','m3gn3gnbspecies_tgn',...
    'm3gn3bgnbspecies_tgn','m3gn4gnbspecies_tgn','m5gngnbgspecies_tgn','m4gngnbgspecies_tgn',...
    'm3gngnbgspecies_tgn','m3gn2gnbgspecies_tgn','m3gn3gnbgspecies_tgn','m3gn3bgnbgspecies_tgn',...
    'm3gn4gnbgspecies_tgn'};

for i=1:length(allglycans)
    allglycans.get(i).id = speciesname{i};
end

for i=1:length(allglycans)
    allglycans.get(i).initConc   = 0;
    allglycans.get(i).initAmount = 0;
end

end

function [allspecies,allrxns] = setGlySpeciesRxninCompt(allspecies,allrxns,...
    listofstructs,thecompt,enzdb)
mgat1  = enzdb('mgat1');
mgat2  = enzdb('mgat2');
mgat3  = enzdb('mgat3');
mgat4  = enzdb('mgat4');
mgat5  = enzdb('mgat5');
mani   = enzdb('mani');
galt   = enzdb('galt');
manii  = enzdb('manii');

%set up rxn
m9species =  GlycanSpecies(listofstructs('M9'),thecompt);
m8species =  GlycanSpecies(listofstructs('M8'),thecompt);
m7species =  GlycanSpecies(listofstructs('M7'),thecompt);
m6species =  GlycanSpecies(listofstructs('M6'),thecompt);
m5species =  GlycanSpecies(listofstructs('M5'),thecompt);
m5gnspecies =  GlycanSpecies(listofstructs('M5gn'),thecompt);
m4gnspecies =  GlycanSpecies(listofstructs('M4gn'),thecompt);
m3gnspecies =  GlycanSpecies(listofstructs('M3gn'),thecompt);
m3gn2species =  GlycanSpecies(listofstructs('M3gn2'),thecompt);
m3gn3species =  GlycanSpecies(listofstructs('M3gn3'),thecompt);
m3gn3bspecies =  GlycanSpecies(listofstructs('M3gn3b'),thecompt);
m3gn4species =  GlycanSpecies(listofstructs('M3gn4'),thecompt);
m5gngspecies =  GlycanSpecies(listofstructs('M5gng'),thecompt);
m4gngspecies =  GlycanSpecies(listofstructs('M4gng'),thecompt);
m3gngspecies =  GlycanSpecies(listofstructs('M3gng'),thecompt);
m3gn2gspecies =  GlycanSpecies(listofstructs('M3gn2g'),thecompt);
m3gn3gspecies = GlycanSpecies(listofstructs('M3gn3g'),thecompt);
m3gn3bgspecies = GlycanSpecies(listofstructs('M3gn3bg'),thecompt);
m3gn4gspecies = GlycanSpecies(listofstructs('M3gn4g'),thecompt);
m5gngnbspecies = GlycanSpecies(listofstructs('M5gngnb'),thecompt);
m4gngnbspecies = GlycanSpecies(listofstructs('M4gngnb'),thecompt);
m3gngnbspecies = GlycanSpecies(listofstructs('M3gngnb'),thecompt);
m3gn2gnbspecies = GlycanSpecies(listofstructs('M3gn2gnb'),thecompt);
m3gn3gnbspecies = GlycanSpecies(listofstructs('M3gn3gnb'),thecompt);
m3gn3bgnbspecies = GlycanSpecies(listofstructs('M3gn3bgnb'),thecompt);
m3gn4gnbspecies = GlycanSpecies(listofstructs('M3gn4gnb'),thecompt);
m5gngnbgspecies = GlycanSpecies(listofstructs('M5gngnbg'),thecompt);
m4gngnbgspecies = GlycanSpecies(listofstructs('M4gngnbg'),thecompt);
m3gngnbgspecies = GlycanSpecies(listofstructs('M3gngnbg'),thecompt);
m3gn2gnbgspecies = GlycanSpecies(listofstructs('M3gn2gnbg'),thecompt);
m3gn3gnbgspecies = GlycanSpecies(listofstructs('M3gn3gnbg'),thecompt);
m3gn3bgnbgspecies = GlycanSpecies(listofstructs('M3gn3bgnbg'),thecompt);
m3gn4gnbgspecies = GlycanSpecies(listofstructs('M3gn4gnbg'),thecompt);

allspecies.add(m9species);
allspecies.add(m8species);
allspecies.add(m7species);
allspecies.add(m6species);
allspecies.add(m5species);
allspecies.add(m5gnspecies);
allspecies.add(m4gnspecies);
allspecies.add(m3gnspecies);
allspecies.add(m3gn2species);
allspecies.add(m3gn3species);
allspecies.add(m3gn3bspecies);
allspecies.add(m3gn4species);
allspecies.add(m5gngspecies);
allspecies.add(m4gngspecies);
allspecies.add(m3gngspecies);
allspecies.add(m3gn2gspecies);
allspecies.add(m3gn3gspecies);
allspecies.add(m3gn3bgspecies);
allspecies.add(m3gn4gspecies);
allspecies.add(m5gngnbspecies);
allspecies.add(m4gngnbspecies);
allspecies.add(m3gngnbspecies);
allspecies.add(m3gn2gnbspecies);
allspecies.add(m3gn3gnbspecies);
allspecies.add(m3gn3bgnbspecies);
allspecies.add(m3gn4gnbspecies);
allspecies.add(m5gngnbgspecies);
allspecies.add(m4gngnbgspecies);
allspecies.add(m3gngnbgspecies);
allspecies.add(m3gn3gnbgspecies);
allspecies.add(m3gn2gnbgspecies);
allspecies.add(m3gn3bgnbgspecies);
allspecies.add(m3gn4gnbgspecies);

m8rxn   = Rxn(m9species,m8species,mani);
allrxns.add(m8rxn);

m7rxn   = Rxn(m8species,m7species,mani);
allrxns.add(m7rxn);

m6rxn   = Rxn(m7species,m6species,mani);
allrxns.add(m6rxn);

m5rxn   = Rxn(m6species,m5species,mani);
allrxns.add(m5rxn);

m5gnrxn = Rxn(m5species,m5gnspecies,mgat1);
allrxns.add(m5gnrxn);

m4gnrxn = Rxn(m5gnspecies,m4gnspecies,manii);
allrxns.add(m4gnrxn);

m3gnrxn = Rxn(m4gnspecies,m3gnspecies,manii);
allrxns.add(m3gnrxn);

m3gn2rxn = Rxn(m3gnspecies,m3gn2species,mgat2);
allrxns.add(m3gn2rxn);

m3gn3brxn = Rxn(m3gn2species,m3gn3bspecies,mgat5);
allrxns.add(m3gn3brxn);

m3gn3rxn  = Rxn(m3gn2species,m3gn3species,mgat4);
allrxns.add(m3gn3rxn);

m3gn4arxn = Rxn(m3gn3bspecies,m3gn4species,mgat4);
allrxns.add(m3gn4arxn);

m3gn4brxn = Rxn(m3gn3species,m3gn4species,mgat5);
allrxns.add(m3gn4brxn);

m5gngrxn = Rxn(m5gnspecies ,m5gngspecies,galt);
allrxns.add(m5gngrxn);

m4gngrxn = Rxn(m4gnspecies,m4gngspecies,galt);
allrxns.add(m4gngrxn);

m3gngrxn = Rxn(m3gnspecies,m3gngspecies,galt);
allrxns.add(m3gngrxn);

m3gn2grxn = Rxn(m3gn2species,m3gn2gspecies,galt);
allrxns.add(m3gn2grxn);

m3gn3bgrxn = Rxn(m3gn3bspecies,m3gn3bgspecies,galt);
allrxns.add(m3gn3bgrxn);

m3gn3grxn = Rxn(m3gn3species,m3gn3gspecies,galt);
allrxns.add(m3gn3grxn);

m3gn4grxn = Rxn(m3gn4species,m3gn4gspecies,galt);
allrxns.add(m3gn4grxn);

m5gngnbrxn = Rxn(m5gnspecies,m5gngnbspecies,mgat3);
allrxns.add(m5gngnbrxn);

m4gngnbrxn = Rxn(m4gnspecies,m4gngnbspecies,mgat3);
allrxns.add(m4gngnbrxn);

m3gngnbrxn = Rxn(m3gnspecies,m3gngnbspecies,mgat3);
allrxns.add(m3gngnbrxn);

m3gn2gnbrxn = Rxn(m3gn2species,m3gn2gnbspecies,mgat3);
allrxns.add(m3gn2gnbrxn);

m3gn3bgnbrxn = Rxn(m3gn3bspecies,m3gn3bgnbspecies,mgat3);
allrxns.add(m3gn3bgnbrxn);

m3gn3gnbrxn = Rxn(m3gn3species,m3gn3gnbspecies,mgat3);
allrxns.add(m3gn3gnbrxn);

m3gn4gnbrxn = Rxn(m3gn4species,m3gn4gnbspecies,mgat3);
allrxns.add(m3gn4gnbrxn);

m5gngnbgrxn = Rxn(m5gngnbspecies,m5gngnbgspecies,galt);
allrxns.add(m5gngnbgrxn);

m4gngnbgrxn = Rxn(m4gngnbspecies,m4gngnbgspecies,galt);
allrxns.add(m4gngnbgrxn);

m3gngnbgrxn = Rxn(m3gngnbspecies,m3gngnbgspecies,galt);
allrxns.add(m3gngnbgrxn);

m3gn2gnbgrxn = Rxn(m3gn2gnbspecies,m3gn2gnbgspecies,galt);
allrxns.add(m3gn2gnbgrxn);

m3gn3bgnbgrxn = Rxn(m3gn3bgnbspecies,m3gn3bgnbgspecies,galt);
allrxns.add(m3gn3bgnbgrxn);

m3gn3gnbgrxn = Rxn(m3gn3gnbspecies,m3gn3gnbgspecies,galt);
allrxns.add(m3gn3gnbgrxn);

m3gn4gnbgrxn = Rxn(m3gn4gnbspecies,m3gn4gnbgspecies,galt);
allrxns.add(m3gn4gnbgrxn);
end

function listofstructs = setGlyStruct()
listofstructs = containers.Map;
listofstructs('M9')=glycanMLread('M9.glycoct_xml');
listofstructs('M8')=glycanMLread('M8.glycoct_xml');
listofstructs('M7')=glycanMLread('M7.glycoct_xml');
listofstructs('M6')=glycanMLread('M6.glycoct_xml');
listofstructs('M5')=glycanMLread('M5.glycoct_xml');
listofstructs('M5gn')=glycanMLread('M5gn.glycoct_xml');
listofstructs('M4gn')=glycanMLread('M4gn.glycoct_xml');
listofstructs('M3gn')=glycanMLread('M3gn.glycoct_xml');
listofstructs('M3gn2')=glycanMLread('M3gn2.glycoct_xml');
listofstructs('M3gn3')=glycanMLread('M3gn3.glycoct_xml');
listofstructs('M3gn3b')=glycanMLread('M3gn3b.glycoct_xml');
listofstructs('M3gn4')=glycanMLread('M3gn4.glycoct_xml');
listofstructs('M5gng')=glycanMLread('M5gng.glycoct_xml');
listofstructs('M4gng')=glycanMLread('M4gng.glycoct_xml');
listofstructs('M3gng')=glycanMLread('M3gng.glycoct_xml');
listofstructs('M3gn2g')=glycanMLread('M3gn2g.glycoct_xml');
listofstructs('M3gn3g')=glycanMLread('M3gn3g.glycoct_xml');
listofstructs('M3gn3bg')=glycanMLread('M3gn3bg.glycoct_xml');
listofstructs('M3gn4g')=glycanMLread('M3gn4g.glycoct_xml');
listofstructs('M5gngnb')=glycanMLread('M5gngnb.glycoct_xml');
listofstructs('M4gngnb')=glycanMLread('M4gngnb.glycoct_xml');
listofstructs('M3gngnb')=glycanMLread('M3gngnb.glycoct_xml');
listofstructs('M3gn2gnb')=glycanMLread('M3gn2gnb.glycoct_xml');
listofstructs('M3gn3gnb')=glycanMLread('M3gn3gnb.glycoct_xml');
listofstructs('M3gn3bgnb')=glycanMLread('M3gn3bgnb.glycoct_xml');
listofstructs('M3gn4gnb')=glycanMLread('M3gn4gnb.glycoct_xml');
listofstructs('M5gngnbg')=glycanMLread('M5gngnbg.glycoct_xml');
listofstructs('M4gngnbg')=glycanMLread('M4gngnbg.glycoct_xml');
listofstructs('M3gngnbg')=glycanMLread('M3gngnbg.glycoct_xml');
listofstructs('M3gn2gnbg')=glycanMLread('M3gn2gnbg.glycoct_xml');
listofstructs('M3gn3gnbg')=glycanMLread('M3gn3gnbg.glycoct_xml');
listofstructs('M3gn3bgnbg')=glycanMLread('M3gn3bgnbg.glycoct_xml');
listofstructs('M3gn4gnbg')=glycanMLread('M3gn4gnbg.glycoct_xml');
end


function allenzs=setupEnzKinetics(enzdb,glycanstructdb)
mgat1  = enzdb('mgat1');
mgat2  = enzdb('mgat2');
mgat3  = enzdb('mgat3');
mgat4  = enzdb('mgat4');
mgat5  = enzdb('mgat5');
mani   = enzdb('mani');
galt   = enzdb('galt');
manii  = enzdb('manii');
m5gnspecies = glycanstructdb('M5gn');
m4gnspecies = glycanstructdb('M4gn');
m3gnspecies = glycanstructdb('M3gn');
m3gn3species = glycanstructdb('M3gn3');
m3gn3bspecies = glycanstructdb('M3gn3b');
m3gn4species = glycanstructdb('M3gn4');
m5gngnbspecies = glycanstructdb('M5gngnb');
m4gngnbspecies = glycanstructdb('M4gngnb');
m3gngnbspecies = glycanstructdb('M3gngnb');
m3gn2gnbspecies = glycanstructdb('M3gn2gnb');
m3gn3bgnbspecies = glycanstructdb('M3gn3bgnb');
m3gn3gnbspecies = glycanstructdb('M3gn3gnb');
m3gn4gnbspecies = glycanstructdb('M3gn4gnb');
    
rt1     = struct('Vm',450,'Km',0.65);
enzkineticsobj = MMenKinetics(TypeStringConst.Comp,TypeStringConst.Ensm,rt1);
mgat1.enzkinetics = enzkineticsobj;

rt2     = struct('Vm',140,'Km',0.475);
enzkineticsobj= MMenKinetics(TypeStringConst.Comp,TypeStringConst.Ensm,rt2);
mgat2.enzkinetics = enzkineticsobj;

rt3     = struct('Vm',0,'Km',0.475);
enzkineticsobj = MMenKinetics(TypeStringConst.Comp,TypeStringConst.Ensm,rt3);
specistruct= CellArrayList;
specistruct.add(m5gnspecies);
specistruct.add(m4gnspecies);
specistruct.add(m3gnspecies);
specifiratio = CellArrayList;
ratio1 = struct('Vm',1,'Km',21.053);
ratio2 = struct('Vm',1,'Km',21.053);
ratio3 = struct('Vm',1,'Km',21.053);
specifiratio.add(ratio1);
specifiratio.add(ratio2);
specifiratio.add(ratio3);
setSubstrateSpecificityK(enzkineticsobj,specistruct,specifiratio,1);
mgat3.enzkinetics = enzkineticsobj;

rt4     = struct('Vm',10,'Km',8.5);
enzkineticsobj = MMenKinetics(TypeStringConst.Comp,TypeStringConst.Ensm,rt4);
mgat4.enzkinetics = enzkineticsobj;

rt5     = struct('Vm',10,'Km',0.325);
enzkineticsobj = MMenKinetics(TypeStringConst.Comp,TypeStringConst.Ensm,rt5);
specistruct= CellArrayList;
specistruct.add(m3gn3species);
specifiratio = CellArrayList;
ratio1 = struct('Vm',1,'Km',0.692);
specifiratio.add(ratio1);
setSubstrateSpecificityK(enzkineticsobj,specistruct,specifiratio,1);
mgat5.enzkinetics = enzkineticsobj;

rt6     = struct('Vm',300,'Km',0.25);
enzkineticsobj6 = MMenKinetics(TypeStringConst.Comp,TypeStringConst.Ensm,rt6);
mani.enzkinetics = enzkineticsobj6;

rt7     = struct('Vm',580,'Km',0.325);
enzkineticsobj = MMenKinetics(TypeStringConst.Comp,TypeStringConst.Ensm,rt7);
specistruct= CellArrayList;
specistruct.add(m5gnspecies);
specistruct.add(m4gnspecies);
specistruct.add(m3gnspecies);
specistruct.add(m3gn3bspecies);
specistruct.add(m3gn3species);
specistruct.add(m3gn4species);
specistruct.add(m5gngnbspecies);
specistruct.add(m4gngnbspecies);
specistruct.add(m3gngnbspecies);
specistruct.add(m3gn2gnbspecies);
specistruct.add(m3gn3bgnbspecies);
specistruct.add(m3gn3gnbspecies);
specistruct.add(m3gn4gnbspecies);
specifiratio = CellArrayList;
ratio1 = struct('Vm',1,'Km',30.769);
ratio2 = struct('Vm',1,'Km',30.769);
ratio3 = struct('Vm',1,'Km',30.769);
ratio4 = struct('Vm',1,'Km',0.538);
ratio5 = struct('Vm',1,'Km',0.385);
ratio6 = struct('Vm',1,'Km',0.308);
ratio7 = struct('Vm',1,'Km',30.769);
ratio8 = struct('Vm',1,'Km',30.769);
ratio9 = struct('Vm',1,'Km',30.769);
ratio10 = struct('Vm',1,'Km',3.846);
ratio11 = struct('Vm',1,'Km',1.692);
ratio12 = struct('Vm',1,'Km',1.538);
ratio13 = struct('Vm',1,'Km',1.077);
specifiratio.add(ratio1);
specifiratio.add(ratio2);
specifiratio.add(ratio3);
specifiratio.add(ratio4);
specifiratio.add(ratio5);
specifiratio.add(ratio6);
specifiratio.add(ratio7);
specifiratio.add(ratio8);
specifiratio.add(ratio9);
specifiratio.add(ratio10);
specifiratio.add(ratio11);
specifiratio.add(ratio12);
specifiratio.add(ratio13);
setSubstrateSpecificityK(enzkineticsobj,specistruct,specifiratio,1);
galt.enzkinetics = enzkineticsobj;

rt8            = struct('Vm',300,'Km',0.25);
enzkineticsobj = MMenKinetics(TypeStringConst.Comp,TypeStringConst.Ensm,rt8);
specistruct    = CellArrayList;
specistruct.add(m5gnspecies);
specifiratio   = CellArrayList;
ratio1         = struct('Vm',1,'Km',2);
specifiratio.add(ratio1);
enzkineticsobj.setSubstrateSpecificityK(specistruct,specifiratio,1);
manii.enzkinetics = enzkineticsobj;

% set up enzyme pool
allenzs = CellArrayList;
allenzs.add(mgat1);
allenzs.add(mgat2);
allenzs.add(mgat3);
allenzs.add(mgat4);
allenzs.add(mgat5);
allenzs.add(mani);
allenzs.add(manii);
allenzs.add(galt);
end

function theRxns = setupRxnKinetics(theRxns,enzdb)
    galt   = enzdb('galt');
  for i = 1 : length(theRxns)
    if(~strcmpi(theRxns.get(i).enz.name, galt.name))
        if(strcmpi(theRxns.get(i).reac.compartment.name,'cis_golgi'))
            theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) = theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) ...
                *0.15;
        elseif(strcmpi(theRxns.get(i).reac.compartment.name,'media_golgi'))
           theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) = theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) ...
                *0.45;
        elseif(strcmpi(theRxns.get(i).reac.compartment.name,'trans_golgi'))
            theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) = theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) ...
                *0.30;
        elseif(strcmpi(theRxns.get(i).reac.compartment.name,'trans_golgi_network'))
            theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) = theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) ...
                *0.10;        
        end
  elseif(theRxns.get(i).enz.name == galt.name)
        if(strcmpi(theRxns.get(i).reac.compartment.name,'cis_golgi'))
              theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) = theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) ...
                *0;
        elseif(strcmpi(theRxns.get(i).reac.compartment.name,'media_golgi'))
              theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) = theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) ...
                *0.05;            
        elseif(strcmpi(theRxns.get(i).reac.compartment.name,'trans_golgi'))
              theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) = theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) ...
                *0.20;
        elseif(strcmpi(theRxns.get(i).reac.compartment.name,'trans_golgi_network'))
              theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) = theRxns.get(i).rxnkinetics.kineticlaws.parametervalues(1) ...
                *0.75;
        end
    end
  end
end