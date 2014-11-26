residueMap=load('residueTypes.mat');

% GCNT1
GCNTI                    = GTEnz([2;4;1;102]);
GCNTI.resfuncgroup       = residueMap.allresidues('GlcNAc');
GCNTI.resAtt2FG          = residueMap.allresidues('GalNAc');
galnacBond               = GlycanBond('?','1');
GCNTI.linkresAtt2FG      = struct('bond',galnacBond,'anomer','a');
glcnacbond                  = GlycanBond('6','1');
GCNTI.linkFG             = struct('anomer','b','bond',glcnacbond);
GCNTI.substMaxStruct     = glycanMLread('518.26.glycoct_xml');

% B3GALT1
B3GALT                     = GTEnz([2;4;1;147]);
B3GALT.resfuncgroup        = residueMap.allresidues('Gal');
B3GALT.resAtt2FG           = residueMap.allresidues('GalNAc');
galnacBond                 = GlycanBond('?','1');
B3GALT. linkresAtt2FG      = struct('bond',galnacBond,'anomer','a');
galbond                    = GlycanBond('3','1');
B3GALT.linkFG              = struct('anomer','b','bond',galbond);
B3GALT.substMinStruct      = glycanMLread('314.16.glycoct_xml');
B3GALT.substMaxStruct      = glycanMLread('314.16.glycoct_xml');
B3GALT.substNABranch       = glycanMLread('559.28b13.glycoct_xml');

% ST3GalNAc I
ST3GalNAc                     = GTEnz([2;4;99;4]);
ST3GalNAc.resfuncgroup        = residueMap.allresidues('NeuAc');
ST3GalNAc.resAtt2FG           = residueMap.allresidues('Gal');
galBond                       = GlycanBond('?','1');  
ST3GalNAc. linkresAtt2FG      = struct('bond',galBond,'anomer','b');
siabond                       = GlycanBond('3','2');
ST3GalNAc.linkFG              = struct('anomer','a','bond',siabond);
ST3GalNAc.isTerminalTarget    = true;

enzArray=CellArrayList;
enzArray.add(B3GALT);
enzArray.add(ST3GalNAc);
enzArray.add(GCNTI);

%define the end prod
a1  = GlycanSpecies(glycanMLread('314.16.glycoct_xml'));
a2  = GlycanSpecies(glycanMLread('1124.56.glycoct_xml'));

glycanArray = CellArrayList;
glycanArray.add(a1);
glycanArray.add(a2);

%Perform reaction
displayOptions = displayset('showMass',true,'showLinkage',true,'showRedEnd',true);
fprintf(1,'Input of glycan product structure is \n');
glycanViewer(a1.glycanStruct,displayOptions);
glycanViewer(a2.glycanStruct,displayOptions);
[isPath,olinkedpath]=inferGlyConnPath(glycanArray,enzArray);
fprintf(1,'Inferred network is shown below:\n');
glycanPathViewer(olinkedpath);