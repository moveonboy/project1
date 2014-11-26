function glycanlist = synGlycanDBBuild(glycantype,specialglycan)
%synGlycanDBBuild dynamically build a local glycan database withiout isomer 
%
%
% Example:
%  enzdbmatfilename = 'glyenzDB.mat';
%  enzdb = enzdbmatLoad(enzdbmatfilename);
%  fprintf(1,'All the enzymes stored in the database are : \n');
%  disp(enzdb.keys);
%  mgat1 = enzdb('mgat1');
%  enzViewer(mgat1);
%
% See also dbLocalConnect,database,close.

% Author: Gang Liu
% Date Lastly Updated: 8/22/13
[path,name,ext]=fileparts(enzdbmatfilename);
if(~strcmpi(ext,'.mat'))
    error('MATLAB:GNAT:WRONGFILETYPE','FILE TYPE SHOULD BE MAT FILE');
end

enzdbstruct = load(enzdbmatfilename);
p = fieldnames(enzdbstruct);
if(length(p)~=1)
    error('MATLAB:GNAT:WRONGINPUT','wrong variable stored');
end

enzdb = enzdbstruct.(p{1});
if(~isa(enzdb,'containers.Map'))
    error('MATLAB:GNAT:WRONGINPUT','wrong variable stored');
end

keynames = enzdb.keys;
for ii = 1 :  length(keynames)
    keyname = keynames{ii};
    ithenz = enzdb(keyname);
    enzprop = properties(ithenz);
    for i = 1:length(enzprop)
        if(isa(ithenz.(enzprop{i}),'CellArrayList'))
            for j = 1 : ithenz.(enzprop{i}).length
                if(isa(ithenz.(enzprop{i}).get(j),'GlycanStruct'))
                    ithenz.(enzprop{i}).get(j).resetjava;
                end
            end
        elseif(isa(ithenz.(enzprop{i}),'GlycanStruct'))
            ithenz.(enzprop{i}).resetjava;
        end
    end
end

end