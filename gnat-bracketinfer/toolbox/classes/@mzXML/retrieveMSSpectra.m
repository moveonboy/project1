function msspectra = retrieveMSSpectra(obj,varargin)
% RETRIEVEMSSPECTRA: Read MS Spectra from an ML file 
%  
% Syntax: 
%    msspectra = retrieveMSSpectra(mzXMLobj)
%    msspectra = retrieveMSSpectra(mzXMLobj,scannum)
%    msspectra = retrieveMSSpectra(mzXMLobj,'mslevel',MSLEVEL) 
%    msspectra = retrieveMSSpectra(mzXMLobj,'scannum',MSLEVEL,'charge',MSCHARGE) 
%  
% Input: 
%   mzXMLobj: an Object of mzXML class
%   scannum: the scan number
%   MSLEVEL: 1 or 2
%   MSCHARGE: z value for MS2 spectra
% 
% Examples:   
%  Example 1:mzXMLobj  = mzXMLread('test1.mzXML');
%            msspectra = retrieveMSSpectra(mzXMLobj,34);
%
%  Example 2: mzXMLobj  = mzXMLread('test1.mzXML');
%             ms1spectras = mzXMLobj.retrieveMSSpectra('mslevel',1);
%
%  Example 3: mzXMLobj  = mzXMLread('test1.mzXML');
%             msspectra = mzXMLobj.retrieveMSSpectra;
%
%See Also mzXML.

% Author: Gang Liu
% Date Lastly Updated: 8/8/14

msspectra = [];

if(nargin==1)
   numscan =  obj.mzxmljava.getScanCount;
   msspectra= cell(numscan,1);
   for i = 1 :numscan
      msspectra{i}=double(obj.mzxmljava.rap(i).getMassIntensityList);
   end
    
elseif(nargin==2)
    % read all scan data in mzXML format    
    msspectra = (obj.mzxmljava.rap(varargin{1}).getMassIntensityList)';
elseif(nargin==3)
    numscan   = obj.mzxmljava.getScanCount;
    if(strcmp(varargin{1},'mslevel')...
            && isnumeric(varargin{2}))
       mslevel = varargin{2};
       count=0;
       for i = 1 : numscan
         if(obj.mzxmljava.rap(i).getMsLevel==mslevel)
            count = count +1;
            msspectra{count} = (obj.mzxmljava.rap(i).getMassIntensityList)'; 
         end
       end       
       
    elseif(strcmp(varargin{1},'scannum')...
            && isnumeric(varargin{2}))   
       scannum = varargin{2};
       for i = 1 : numscan
         if(obj.mzxmljava.rap(i).getNum==scannum)
            msspectra{1} =(obj.mzxmljava.rap(i).getMassIntensityList)';
            break
         end
       end
    else
       error('MATLAB:GPAT:INPUTERROR','WRONG INPUT TYPE');
    end
elseif(nargin==5)
    numscan   = obj.mzxmljava.getScanCount;
    if(strcmpi(varargin{1},'scannum') && strcmpi(varargin{3},'charge'))
        for i =1: numscan
            if (varargin{2}==obj.mzxmljava.rap(i).getNum) && ...
                (varargin{4}==obj.mzxmljava.rap(i).getPrecursorCharge)               
                msspectra{1} =(obj.mzxmljava.rap(i).getMassIntensityList)';
                break;
            end
        end
    else
        error('MATLAB:GPAT:INPUTERROR','WRONG INPUT TYPE');
    end
end

if(length(msspectra)==1)
    msspectra = msspectra{1};
    msspectra = double(msspectra);
end

end