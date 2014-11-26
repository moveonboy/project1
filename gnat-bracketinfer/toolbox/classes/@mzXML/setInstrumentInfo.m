function setInstrumentInfo(obj)
% SETINSTRUMENTINFO read MS instrument information
%   
% Syntax: setInstrumentInfo(obj)
%         obj.setInstrumentInfo
%
%See Also mzXML;

% Author: Gang Liu
% Date Lastly Updated: 1/26/14

mzxmljava        = obj.mzxmljava;
instrinfojava    = mzxmljava.getHeaderInfo.getInstrumentInfo;
if(~isempty(instrinfojava))
obj.instrumentinfo.dectector    = char(instrinfojava.getDetector);
obj.instrumentinfo.ionization   = char(instrinfojava.getIonization);
obj.instrumentinfo.manufacturer = char(instrinfojava.getManufacturer);
obj.instrumentinfo.massAnalyzer = char(instrinfojava.getMassAnalyzer);
obj.instrumentinfo.model        = char(instrinfojava.getModel);
obj.instrumentinfo.operator     = char(instrinfojava.getOperator);
obj.instrumentinfo.softwarename = char(instrinfojava.getSoftwareInfo.name);
obj.instrumentinfo.softwaretype = char(instrinfojava.getSoftwareInfo.type);
obj.instrumentinfo.softwareversion = char(instrinfojava.getSoftwareInfo.version);
else
   obj.instrumentinfo=[]; 
end

end