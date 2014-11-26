function mzXMLobj = readmzXML(varargin)
%READMZXML: Read mass spectrometry data from an mzxML file 
%
% Syntax: 
%     mzxmlobj = readmzXML(mzxmlfilefullname);
%     mzxmlobj = readmzXML(mzxmlfiledir, mzxmlfilename);
%   
% Input: 
%    mzxmlfilefullname: the full name for an mzXML file including its path 
%    mzxmlfiledir: the name of the file directory
%    filename: the name of an mzXML file 
%   (Note: if the file "test1.mzxml" is stored in the folder
%      'c:/glycopat/toolbox/test/data'
%      the "mzxmlfilefullname" is 'test1.mzxml'
%      the  "mzxmlfiledir" is 'c:/glycopat/toolbox/test/data'
%      the  "filename " is 'test1.mzxml')
%    
% Output:
%   mzxmlobj: an object of mzXML class.
%   
% Examples:
%
%   Example 1. 
%     mzxmlobj = readmzXML('c:/glycopat/toolbox/test/data',test1.mzXML);
%
%   Example 2. 
%     mzxmlobj = readmzXML('c:/glycopat/toolbox/test/data/test1.mzXML');
%  
%   Example 3. 
%       mzxmlobj = readmzXML('test1.mzXML'). The program will search 
%    for "test1.mzXML" in the current path.If not found,it will look for
%    the file in MATLAB search path.  
%
%See also mzDTAread,peptideread,varptmread,fixedptmread. 

% Author: Gang Liu
% Date Lastly Updated: 08/08/14

narginchk(1,2);

if(nargin==1)
     mzxmlfilename  = varargin{1};  
elseif(nargin==2)
     filepath = varargin{1};
     filename = varargin{2}; 
     mzxmlfilename = fullfile(filepath,filename);
end

fid = fopen(mzxmlfilename);
if(fid==-1) % if file not found
    mzxmlfilename = which(mzxmlfilename);
    if(isempty(mzxmlfilename))
        error('MATLAB:GlycoPAT:FILENOTFOUND','file is not found in the search path');
    end
end
fclose(fid);

import org.systemsbiology.jrap.*;
mzXMLobj = mzXML(mzxmlfilename,0,'memsave');

end