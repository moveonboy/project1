function [mz]= retrieveMz(obj,varargin)
% GETMSSpectra read MS Scan Data from mzXML file 
%   
% Example 1: mzXMLobj = mzXML('test2.mzXML',1);
%            [mz] = mzXMLobj.retrieveMz(23);
%
%See Also mzXML.

% Author: Gang Liu
% Date Lastly Updated: 1/26/14

mz =[];
if(nargin==1)
   numscan =  obj.mzxmljava.getScanCount;
   mz=zeros(numscan,1);
   for i = 1 :numscan
      mz(i)=double(obj.mzxmljava.rap(i).getPrecursorMz);
   end
elseif(nargin==2)
  mz=double(obj.mzxmljava.rap(varargin{1}).getPrecursorMz);
elseif(nargin==3)
  if(strcmpi(varargin{1},'scannum'))
    numscan =  obj.mzxmljava.getScanCount;
    for i = 1 : numscan
      if(obj.mzxmljava.rap(i).getNum==varargin{2})
        mz = double(obj.mzxmljava.rap(i).getPrecursorMz);
        return
      end
    end
  else
    error('MATLAB:GPAT:ERRORINPUT','WRONG INPUT');
  end
  error('MATLAB:GPAT:ERRORINPUT','WRONG INPUT');
end