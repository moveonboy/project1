residueMap                 = load('residueTypes.mat');
B3GALTIV                     = GTEnz([2;4;1;86]);
B3GALTIV.resfuncgroup        = residueMap.allresidues('Gal');
B3GALTIV.resAtt2FG           = residueMap.allresidues('GlcNAc');
glcnacBond                 = GlycanBond('?','1');
B3GALTIV. linkresAtt2FG      = struct('bond',glcnacBond,'anomer','b');
galbond                    = GlycanBond('4','1');
B3GALTIV.linkFG              = struct('anomer','b','bond',galbond);
enzViewer(B3GALTIV);