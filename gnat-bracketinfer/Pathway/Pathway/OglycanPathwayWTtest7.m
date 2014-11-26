residueMap=load('residueTypes.mat');

B3GALTI                     = GTEnz([2;4;1;122]);
B3GALTI.resfuncgroup        = residueMap.allresidues('Gal');
B3GALTI.resAtt2FG           = residueMap.allresidues('GalNAc');
galnacBond                 = GlycanBond('?','1');
B3GALTI. linkresAtt2FG      = struct('bond',galnacBond,'anomer','a');
galbond                    = GlycanBond('3','1');
B3GALTI.linkFG              = struct('anomer','b','bond',galbond);
B3GALTI.substMaxStruct      = glycanMLread('314.16.glycoct_xml');

B3GALTIV                     = GTEnz([2;4;1;86]);
B3GALTIV.resfuncgroup        = residueMap.allresidues('Gal');
B3GALTIV.resAtt2FG           = residueMap.allresidues('GlcNAc');
glcnacBond                 = GlycanBond('?','1');
B3GALTIV. linkresAtt2FG      = struct('bond',glcnacBond,'anomer','b');
galbond                    = GlycanBond('4','1');
B3GALTIV.linkFG              = struct('anomer','b','bond',galbond);

FUT4                     = GTEnz([2;4;1;152]);
FUT4.resfuncgroup        = residueMap.allresidues('Fuc');
FUT4.resAtt2FG           = residueMap.allresidues('GlcNAc');
glcnacBond                  = GlycanBond('?','1');
FUT4.linkresAtt2FG       = struct('bond',glcnacBond,'anomer','b');
futbond                  = GlycanBond('3','1');
FUT4.linkFG              = struct('anomer','a','bond',futbond);
FUT4.targetBranch        = glycanMLread('559.28b16.glycoct_xml');
FUT4.substMinStruct      = glycanMLread('967.48.glycoct_xml');
FUT4.isTerminalTarget    = false;

FUT7                     = GTEnz([2;4;1;152]);
FUT7.resfuncgroup        = residueMap.allresidues('Fuc');
FUT7.resAtt2FG           = residueMap.allresidues('GlcNAc');
glcnacBond                  = GlycanBond('?','1');
FUT7.linkresAtt2FG       = struct('bond',glcnacBond,'anomer','b');
futbond                  = GlycanBond('3','1');
FUT7.linkFG              = struct('anomer','a','bond',futbond);
FUT7.targetBranch        = glycanMLread('559.28b16.glycoct_xml');
FUT7.isTerminalTarget    = false;

FUT9                     = GTEnz([2;4;1;152]);
FUT9.resfuncgroup        = residueMap.allresidues('Fuc');
FUT9.resAtt2FG           = residueMap.allresidues('GlcNAc');
glcnacBond                  = GlycanBond('?','1');
FUT9.linkresAtt2FG       = struct('bond',glcnacBond,'anomer','b');
futbond                  = GlycanBond('3','1');
FUT9.linkFG              = struct('anomer','a','bond',futbond);
FUT9.targetBranch = glycanMLread('559.28b16.glycoct_xml');
FUT9.substMinStruct      = glycanMLread('967.48.glycoct_xml');
FUT9.isTerminalTarget    = false;

GCNTI                    = GTEnz([2;4;1;102]);
GCNTI.resfuncgroup       = residueMap.allresidues('GlcNAc');
GCNTI.resAtt2FG          = residueMap.allresidues('GalNAc');
galnacBond               = GlycanBond('?','1');
GCNTI.linkresAtt2FG      = struct('bond',galnacBond,'anomer','a');
glcnacbond               = GlycanBond('6','1');
GCNTI.linkFG             = struct('anomer','b','bond',glcnacbond);
GCNTI.substMaxStruct     = glycanMLread('518.26.glycoct_xml');
GCNTI.substNABranch      = glycanMLread('879.43.glycoct_xml');

ST3GalI                     = GTEnz([2;4;99;4]);
ST3GalI.resfuncgroup        = residueMap.allresidues('NeuAc');
ST3GalI.resAtt2FG           = residueMap.allresidues('Gal');
galBond                     = GlycanBond('?','1');  
ST3GalI. linkresAtt2FG      = struct('bond',galBond,'anomer','b');
siabond                     = GlycanBond('3','2');
ST3GalI.linkFG              = struct('anomer','a','bond',siabond);
ST3GalI.targetBranch        = glycanMLread('518.26.glycoct_xml');
ST3GalI.isTerminalTarget    = true;

ST3GalIII                     = GTEnz([2;4;99;6]);
ST3GalIII.resfuncgroup        = residueMap.allresidues('NeuAc');
ST3GalIII.resAtt2FG           = residueMap.allresidues('Gal');
galBond                       = GlycanBond('?','1');  
ST3GalIII. linkresAtt2FG      = struct('bond',galBond,'anomer','b');
siabond                       = GlycanBond('3','2');
ST3GalIII.linkFG              = struct('anomer','a','bond',siabond);
ST3GalIII.targetBranch        = glycanMLread('763.38tb.glycoct_xml');
ST3GalIII.isTerminalTarget    = true;

ST3GalIV                     = GTEnz([2;4;99;4]);
ST3GalIV.resfuncgroup        = residueMap.allresidues('NeuAc');
ST3GalIV.resAtt2FG           = residueMap.allresidues('Gal');
galBond                     = GlycanBond('?','1');  
ST3GalIV. linkresAtt2FG      = struct('bond',galBond,'anomer','b');
siabond                     = GlycanBond('3','2');
ST3GalIV.linkFG              = struct('anomer','a','bond',siabond);
ST3GalIV.targetBranch        = glycanMLread('763.38tb.glycoct_xml');
ST3GalIV.isTerminalTarget    = true;

ST6GalNAc                    = GTEnz([2;4;99;3]);
ST6GalNAc.resfuncgroup       = residueMap.allresidues('NeuAc');
ST6GalNAc.resAtt2FG          = residueMap.allresidues('GalNAc');
glcnacBond                   = GlycanBond('?','1');
ST6GalNAc.linkresAtt2FG      = struct('bond',glcnacBond,'anomer','a');
siabond                      = GlycanBond('6','2');
ST6GalNAc.linkFG             = struct('anomer','a','bond',siabond);
ST6GalNAc.substMinStruct     = glycanMLread('314.16.glycoct_xml');
ST6GalNAc.substMaxStruct     = glycanMLread('879.43.glycoct_xml');

enzArray=CellArrayList;
enzArray.add(B3GALTI);
enzArray.add(B3GALTIV);
enzArray.add(FUT4);
enzArray.add(FUT7);
enzArray.add(FUT9);
enzArray.add(GCNTI);
enzArray.add(ST3GalI);
enzArray.add(ST3GalIII);
enzArray.add(ST3GalIV);
enzArray.add(ST6GalNAc);

%define the end prod
% a1  = GlycanSpecies(glycanMLread('314.16.glycoct_xml'));
a2  = GlycanSpecies(glycanMLread('1863.92.glycoct_xml'));
a3  = GlycanSpecies(glycanMLread('1502.75.glycoct_xml'));



glycanArray = CellArrayList;
% glycanArray.add(a1);
glycanArray.add(a2);
glycanArray.add(a3);

%Perform reaction
displayOptions = displayset('showMass',true,'showLinkage',true,'showRedEnd',true);
fprintf(1,'Input of glycan product structure is \n');
% glycanViewer(a1.glycanStruct,displayOptions);
glycanViewer(a2.glycanStruct,displayOptions);
glycanViewer(a3.glycanStruct,displayOptions);
[isPath,olinkedpath]=inferGlyConnPath(glycanArray,enzArray);
fprintf(1,'Inferred network is shown below:\n');
glycanPathViewer(olinkedpath);