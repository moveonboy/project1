function [zcharge]= retrievezCharge(obj)
% GETMSSpectra read MS Scan Data from mzXML file 
%   
% Example 1: mzXMLobj = mzXML('test2.mzXML',1);
%            [scan] = mzXMLobj.retrieveScanNum;
%
%See Also mzXML.

% Author: Gang Liu
% Date Lastly Updated: 1/26/14

zcharge =[];
if(nargin==1)
   numscan =  obj.mzxmljava.getScanCount;
   zcharge=zeros(numscan,1);
   for i = 1 :numscan
      zcharge(i)=double(obj.mzxmljava.rap(i).getPrecursorCharge);
   end
elseif(nargin==2)
   zcharge=double(obj.mzxmljava.rap(varargin{1}).getPrecursorCharge);
end

end