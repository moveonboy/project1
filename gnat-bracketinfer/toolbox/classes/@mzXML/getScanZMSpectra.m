function [scan,z,Mz,spectra]= getScanZMSpectra(obj)
% GETMSSpectra read MS Scan Data from mzXML file 
%   
% Example 1: mzXMLobj = mzXML('test2.mzXML',1);
%            [scan,z,mz,spectra] = mzXMLobj.getScanZMSpectra;
%
%See Also mzXML.

% Author: Gang Liu
% Date Lastly Updated: 1/26/14

ms2scan = obj.getMSScan('mslevel',2); 
for i = 1 : length(ms2scan)
   scan(i)    = ms2scan{i}.scannum;
   z(i)       = ms2scan{i}.precursorCharge;
   spectra{i} = ms2scan{i}.msintensitylist;
   Mz(i)      = ms2scan{i}.precursorMz;
end        

end