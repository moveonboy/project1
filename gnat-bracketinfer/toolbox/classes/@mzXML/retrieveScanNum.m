function [scan]= retrieveScanNum(obj,varargin)
% retrieveScanNum read scan number from mzXML file 
%   
% Example 1: mzXMLobj = mzXML('test2.mzXML',0,'memsave');
%            [scan] = mzXMLobj.retrieveScanNum;
%
%See Also mzXML.

% Author: Gang Liu
% Date Lastly Updated: 1/26/14

scan =[];
if(nargin==1)
   numscan =  obj.mzxmljava.getScanCount;
   scan=zeros(numscan,1);
   for i = 1 :numscan
      scan(i)=double(obj.mzxmljava.rap(i).getNum);
   end
elseif(nargin==2)
  scan=double(obj.mzxmljava.rap(varargin{1}).getNum);
end

end