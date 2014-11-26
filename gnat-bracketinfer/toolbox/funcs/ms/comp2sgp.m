function sgpstring=comp2sgp(glycompstring)
%COMP2SGP: Convert glycan composition string to small glypep format
%   
% Syntax:
%   sgpstring = comp2sgp(glycompstring)
%  
% Input: 
%   glycompstring: a string of glycan composition
%
% Output:
%   sgpstring: small glypep format.
% 
% Examples:
%   Example 1:
%      glycan3letterstring='Hex5HexNAc2'; % Hex5HexNAc2
%      sgpstring=comp2sgp(glycan3letterstring)
%  
%  Example 2:
%     glycan3letterstring ='Hex25HexNAc25';
%     sgpstring=comp2sgp(glycan3letterstring)  
%     
% See also glycomp.

% Author: Gang Liu
% Date Lastly Updated: 10/31/2014 by Gang Liu

sgpstring =[];

% HexNAC -h
nHexNAC        =  regexp(glycompstring,'(HexNAc\d+)','match');
numHexNAc      =  str2double(nHexNAC{1,1}(7:end));
for i = 1 : numHexNAc
    sgpstring = strcat(sgpstring,'n');
end

% Hex
nHex       =  regexp(glycompstring,'(Hex\d+)','match');
numHex     =  str2double(nHex{1,1}(4:end));
for i = 1 : numHex
    sgpstring = strcat(sgpstring,'h');
end

% Sailic Acid NeuAc
% nNeuAc=findstr('s',gly);
nNeuAc             =  regexp(glycompstring,'(NeuAc\d+)','match');

if(~isempty(nNeuAc))
numNeuAc       =  str2double(nNeuAc{1,1}(6:end));
for i = 1 : numNeuAc
    sgpstring = strcat(sgpstring,'s');
end

end

% Sailic Acid NeuGc
% nNeuAc=findstr('g',gly);
nNeuGc              =  regexp(glycompstring,'(NeuGc\d+)','match') ;
if(~isempty(nNeuGc))
numNeuGc       =  str2double(nNeuGc{1,1}(6:end));
for i = 1 : numNeuGc
    sgpstring = strcat(sgpstring,'g');
end
end

% Fucose f
nFuc               =  regexp(glycompstring,'(Fuc\d+)','match') ;  
if(~isempty(nFuc))

numFuc  =  str2double(nFuc{1,1}(4:end));
for i = 1 : numFuc
    sgpstring = strcat(sgpstring,'f');
end

end

% pentose
nPen           = regexp(glycompstring,'(Pen\d+)','match') ;
if(~isempty(nPen))

numPen         =  str2double(nPen{1,1}(4:end)); 
for i = 1 : numPen
    sgpstring = strcat(sgpstring,'x');
end

end

end

