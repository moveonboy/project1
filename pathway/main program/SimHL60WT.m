function therofraction = SimHL60WT()
%load glycanNetReactor Model
glymodel = GlycanNetModel.loadmat('HL60Model.mat');
odemfile = 'HL60WT.m';

% apply glycanSSim
[exitflag, speciesconc, speciesnames,funceval] = ...
    glycanSSim(glymodel,odemfile);

speciesconcdb = containers.Map;
for i = 1 : length(speciesnames)
  speciesconcdb(speciesnames{i}) = speciesconc(i);
end

save('specconcdb.mat','speciesconcdb');

glycangroupexptHL60 = createGlycanHL60Input();
numgroups       = glycangroupexptHL60.length;
theroabundance  = zeros(numgroups,1);

lastcomptname = glymodel.glycanpathway.findLastCompt.name;


totalabundance = 0;
for i = 1 : numgroups;
    ithgroup = glycangroupexptHL60.get(i);
    numglycansinithgroup = length(ithgroup.glycanlist);
    glycanindexinithgroup = [];
    for j = 1 : numglycansinithgroup
       jthglycaninithgroup = ithgroup.glycanlist(j,1);
       [jthspecies,speciesindex] = glymodel.glycanpathway.findSpeciesByStruct(...
           jthglycaninithgroup,lastcomptname);
       glycanindexinithgroup = [glycanindexinithgroup;speciesindex];
    end
    
    theroabundance(i) = sum(speciesconc(glycanindexinithgroup));
    totalabundance    = totalabundance+theroabundance(i);
end
therofraction = theroabundance/totalabundance;
end