%% Output HL60WT Model
HL60WTModelSBMLfileName = 'test2.xml';
HL60WTModel = glycanNetSBMLread(HL60WTModelSBMLfileName);
HL60WTModel.glycanNet_sbmlmodel = convertToSingleCompt(HL60WTModel.glycanNet_sbmlmodel);
ubODEFileName = 'test2';
WriteODEFunction(HL60WTModel.glycanNet_sbmlmodel,ubODEFileName);