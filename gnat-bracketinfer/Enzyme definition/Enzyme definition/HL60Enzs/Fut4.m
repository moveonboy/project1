residueMap               = load('residueTypes.mat');
FUT4                     = GTEnz([2;4;1;152]);
FUT4.resfuncgroup        = residueMap.allresidues('Fuc');
FUT4.resAtt2FG           = residueMap.allresidues('GlcNAc');
glcnacBond               = GlycanBond('?','1');
FUT4.linkresAtt2FG       = struct('bond',glcnacBond,'anomer','b');
futbond                  = GlycanBond('3','1');
FUT4.linkFG              = struct('anomer','a','bond',futbond);
% FUT4.targetbranchcontain = glycanMLread('518.26fut4.glycoct_xml');
FUT4.targetBranch        = glycanMLread('559.28b16.glycoct_xml');
FUT4.substMinStruct      = glycanMLread('967.48.glycoct_xml');
FUT4.isTerminalTarget    = false;
enzViewer(FUT4);

