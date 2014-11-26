function actmethod = retrieveActMethod(obj,varargin)
%RETRIEVEACTMETHOD: Read MS activation method from an mzXML object
%  
% Syntax:
%   actmethod = retrieveActMethod(mzXMLObj,ithscan)
%   actmethod = retrieveActMethod(mzXMLObj,'scannum',scannum) 
%  
% Example: 
%   
%   Example 1:
%     mzXMLobj  = mzXMLread('fetuin_test.mzXML'');
%     actmethod = mzXMLobj.retrieveActMethod(34);
%     disp(actmethod)
%     Answer: 
%       'ETD+SA'
%
%  Example 2:  
%     mzXMLobj  = mzXMLread('fetuin_test.mzXML'');
%     actmethod = mzXMLobj.retrieveActMethod('scannum',3500);
%   
%
%See Also mzXML.

% Author: Gang Liu
% Date Lastly Updated: 8/8/14

actmethod = [];

if(nargin==2)
    % read all scan data in mzXML format
    actmethod = char(obj.extmzxmljava.rap(varargin{1}).getActivationMethod);
elseif(nargin==3)
    numscan   = obj.mzxmljava.getScanCount;
    if(strcmpi(varargin{1},'scannum')...
            && isnumeric(varargin{2}))   
       scannum = varargin{2};
       for i = 1 : numscan
         if(obj.mzxmljava.rap(i).getNum==scannum)
            actmethod =char(obj.extmzxmljava.rap(i).getActivationMethod); 
            return
         end
       end
    else
       error('MATLAB:GLYCOPAT:INCORRECTINPUT','INCORRECT INPUT TYPE');
    end
end

end