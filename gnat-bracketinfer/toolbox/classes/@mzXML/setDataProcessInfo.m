function setDataProcessInfo(obj)
% SETDATAPROCESSINFO read MS data processing information
%
% Syntax: obj.setDataProcessInfo
%         setDataProcessInfo(obj)
%
%See Also mzXML.

% Author: Gang Liu
% Date Lastly Updated: 9/25/13

mzxmljava           = obj.mzxmljava;
dataprocessingjava  = mzxmljava.getHeaderInfo.getDataProcessing;
if(~isempty(dataprocessingjava))
    obj.dataprocessinginfo.centroided   = dataprocessingjava.getCentroided;
    obj.dataprocessinginfo.deconvoluted = dataprocessingjava.getChargeDeconvoluted;
    obj.dataprocessinginfo.deisotoped   = dataprocessingjava.getDeisotoped;
    obj.dataprocessinginfo.intensitycufoff = dataprocessingjava.getIntensityCutoff;
    obj.dataprocessinginfo.softwareused =  char(dataprocessingjava.getSoftwareUsed.toString);
    obj.dataprocessinginfo.spotintegration = dataprocessingjava.getSpotIntegration;
else
    obj.dataprocessinginfo=[];
end
end