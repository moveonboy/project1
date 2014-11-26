function mzInt = readMSD(varargin)
%readMS convert MS raw peak data to a two-column 
%  matrix comprising of mz and intensity value. 
%
% mzInt = readMSD(MSDFULLFILENAME) reads MS data in msd file
%  format and returns a n*2 matrix containing mz and 
%  intensity value.
%
%
% mzInt = readMSD(MSDFILEPATH,MSDFILENAME) reads MS data in msd file
%  format and returns a n*2 matrix containing mz and 
%  intensity value.
% 
% Example: 
%       mzInt = readMSD('testCHO.msd');
%       mzInt = readMSD('c:\gnat\toolbox\test','testCHO.msd');
%       
%See also MSPROCESS.

% Author: Gang Liu
% Last Date Updated: 10/31/14 2:10 PM 

narginchk(1,2);

if(nargin==1)
    msdfullfilename = varargin{1};
elseif(nargin==2)
    msdfullfilepath = varargin{1};
    msdfilename = varargin{2};
    msdfullfilename = fullfile(msdfullfilepath,msdfilename);
end

if(exist(msdfullfilename,'file')~=2)
    error('MATLAB:GNAT:FILENOTFOUND','FILE IS NOT FOUND');
end

[pathstr,name,ext] = fileparts(msdfullfilename);
if(isempty(pathstr))
    msdfullfilename = which(msdfullfilename);
    if(isempty(msdfullfilename))
       error('MATLAB:GNAT:FILENOTFOUND','FILE IS NOT FOUND'); 
    end
end

mzInt= dlmread(msdfullfilename, '\t', 2, 0); % Skip the first two lines of RAW Data file. 

