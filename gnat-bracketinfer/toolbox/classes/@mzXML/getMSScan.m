function msdata = getMSScan(obj,varargin)
% GETMSSpectra read MS Scan Data from mzXML file 
%   
% Example 1: mzXMLobj = mzXML('test1.mzXML');
%            scandata = mzXMLobj.getMSScan(1);
%
% Example 2: mzXMLobj = mzXML('test1.mzxml');
%            scandata = mzXMLobj.getMSScan('mslevel',1);
%
% Example 3: mzXMLobj = mzXML('test1.mzxml');
%            scandata = mzXMLobj.getMSScan('scannum',40);
%
%See Also mzXML.

% Author: Gang Liu
% Date Lastly Updated: 1/26/14

narginchk(2,3);

if(nargin==2)
    % read all scan data in mzXML format
    msdata    = obj.scandata{varargin{1}};
elseif(nargin==3)
    numscan   = obj.mzxmljava.getScanCount;
    if(strcmpi(varargin{1},'mslevel')...
            && isnumeric(varargin{2}))
       mslevel = varargin{2};
       count=0;
       for i = 1 : numscan
         if(obj.scandata{i}.msLevel==mslevel)
            count = count +1;
            msdata{count} = obj.scandata{i}; 
         end
       end
    elseif(strcmp(varargin{1},'scannum')...
            && isnumeric(varargin{2}))   
       scannum = varargin{2};
       count=0;
       for i = 1 : numscan
         if(obj.scandata{i}.scannum==scannum)
            count = count +1;
            msdata{count} = obj.scandata{i}; 
         end
       end
    else
       msdata = [];
       return;
    end
end

% if(length(msdata)==1)
%     msdata=msdata{:};
% end
end