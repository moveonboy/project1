function gly3letstr=readMSA(varargin)
%READMSA: Read MSA file annoated using glycoworkbench and Cartoonist
% 
% Syntax: 
%    gly3letstr = readMSA(msafullfilename)
%    gly3letstr = readMSA(msafilepath,msafilename)
%
% Input:
%    msafullfilename: the full name of a MSA file including the path
%    msafilepath: the path of a MSA file
%
% Output:
%    gly3letstr: a cell array of glycan composition strings
%
% Example 1:
%         msafilename='testmsa.txt';
%         gly3letstr=readMSA(msafilename);
%
%See also readMSD.

% Author: Yushen Zhou & Gang Liu 
% Date Lastly Updated, 10/31/14 by Gang Liu

narginchk(1,2);
if(nargin==1)
    msafullfilename = varargin{1};
elseif(nargin==2)
    msafullfilepath = varargin{1};
    msafilename = varargin{2};
    msafullfilename = fullfile(msafullfilepath,msafilename);
end

if(exist(msafullfilename,'file')~=2)
    error('MATLAB:GNAT:FILENOTFOUND','FILE IS NOT FOUND');
end

[pathstr,name,ext] = fileparts(msafullfilename);
if(isempty(pathstr))
    msafullfilename = which(msafullfilename);
    if(isempty(msafullfilename))
       error('MATLAB:GNAT:FILENOTFOUND','FILE IS NOT FOUND'); 
    end
end

fid=fopen(msafullfilename);
datatable=textscan(fid,...
    '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %*[^\n]',...
    'headerlines',2);
gly1letstr=datatable{1,16};
gly3letstr = cell(length(gly1letstr),1);
for i=1:length(gly1letstr)
    glystring_msa=gly1letstr{i,1};
    numHex=length(regexp(glystring_msa,'M','match'))+ ...
        length(regexp(glystring_msa,'A','match'))+...
        length(regexp(glystring_msa,'Ga','match'));    
    numHexNac=length(regexp(glystring_msa,'GN','match'));
    numNeuAc=length(regexp(glystring_msa,'NN','match'));    
    numNeuGc=length(regexp(glystring_msa,'NJ','match'));
    numFuc=length(regexp(glystring_msa,'Fa','match'));
    gly3letstr{i,1} = composmake(numHex,numHexNac,numNeuAc,numNeuGc,numFuc);    
end
fclose(fid);

end

function glycomp = composmake(numHex,numHexNac,numNeuAc,numNeuGc,numFuc)    
glycomp = '';  
if(numHex>0)
   glycomp=[glycomp 'Hex' int2str(numHex)]; 
end

if(numHexNac>0)
   glycomp=[glycomp 'HexNAc' int2str(numHexNac)]; 
end

if(numNeuAc>0)
   glycomp=[glycomp 'NeuAc' int2str(numNeuAc)]; 
end

if(numNeuGc>0)
   glycomp=[glycomp 'NeuGc' int2str(numNeuGc)]; 
end

if(numFuc>0)
  glycomp=[glycomp 'Fuc' int2str(numFuc)];   
end

end