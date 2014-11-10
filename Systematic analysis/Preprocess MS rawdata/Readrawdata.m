function data = Readrawdata(varargin)
%Readrawdata covert the format '.txt' into Matlab format using readMS. 
%
% data = Readrawdata(varargin) if input is empty using defalut rawdata adress
% and filename, otherwise, input the file adress and name to convert into Matlab
% format.
%
% Example:
%     data = Readrawdata('C:\Users\Yusen\Dropbox\HL60cells\N_glycan\Rawdata',...
%               'HL60ST3Gal4CRISPR50%N-Glycans.txt') 
%
%See also readMS.

% Author: Yusen Zhou & Gang Liu
% Date Lastly Updated: 9/22/14
narginck(0,2)

if(isempty(varargin))
    msrawdatapath     = 'C:\Users\Yusen\Dropbox\HL60cells\N_glycan\Rawdata';
    msrawdatafilename ='HL60ST3Gal4CRISPR50%N-Glycans.txt';
elseif(length(varargin)==1)
    msrawdatapath     = 'C:\Users\Yusen\Dropbox\HL60cells\N_glycan\Rawdata';
    msrawdatafilename = varargin{1};
elseif(length(varargin)==2)
    msrawdatapath     = varargin{1};
    msrawdatafilename = varargin{2};
end
    
msrawdatafullfilename=fullfile(msrawdatapath,msrawdatafilename);
data = readMS(msrawdatafullfilename);
save('rawdata.mat','data');
end