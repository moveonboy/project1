function mslevel = retrieveMSlevel(obj,varargin)
% retrieveMSlevel read MS level from mzXML file 
%   
% Example 1: mzXMLobj  = mzXML('test1.mzXML',0,'memsave');
%            mslevel   = mzXMLobj.retrieveMSlevel(1);
%
% Example 2: mzXMLobj  = mzXML('test1.mzxml',0,'memsave');
%            mslevel = mzXMLobj.retrieveMSlevel('scannum',34);
%
%See Also mzXML.

% Author: Gang Liu
% Date Lastly Updated: 3/6/14

mslevel = -1;

if(nargin==2)
    % read all scan data in mzXML format
    mslevel = double(obj.extmzxmljava.rap(varargin{1}).getMsLevel);
elseif(nargin==3)
    numscan   = obj.mzxmljava.getScanCount;
    if(strcmpi(varargin{1},'scannum')...
            && isnumeric(varargin{2}))   
       scannum = varargin{2};
       for i = 1 : numscan
         if(obj.mzxmljava.rap(i).getNum==scannum)
            mslevel =double(obj.extmzxmljava.rap(i).getMsLevel); 
            return
         end         
       end
       error('MATLAB:GPAT:INCORRESCANNUM','SCAN NUMBER IS NOT FOUND IN MZXML');
    else
       error('MATLAB:GPAT:INCORRECTINPUT','INCORRECT INPUT TYPE');
    end
end

end