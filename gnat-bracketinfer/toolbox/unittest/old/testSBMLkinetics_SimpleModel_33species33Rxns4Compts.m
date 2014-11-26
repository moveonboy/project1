% test example 
function simBailey()
clc;clear

theCompt   = CellArrayList;
cisCompt   = Compt('cis_golgi');
mediaCompt = Compt('media_golgi');
transCompt = Compt('trans_golgi');
tgnCompt   = Compt('trans_golgi_network');
cisCompt.posteriorcompt   = mediaCompt;
mediaCompt.priorcompt     = cisCompt;
mediaCompt.posteriorcompt = transCompt;
transCompt.priorcompt     = mediaCompt;
transCompt.posteriorcompt = tgnCompt;
tgnCompt.priorcompt       = transCompt;
theCompt.add(cisCompt);
theCompt.add(mediaCompt);
theCompt.add(transCompt);
theCompt.add(tgnCompt);

m8species_cis         = GlycanSpecies(glycanMLread('M8.glycoct_xml'),cisCompt);
m7species_cis         = GlycanSpecies(glycanMLread('M7.glycoct_xml'),cisCompt);
m6species_cis         = GlycanSpecies(glycanMLread('M6.glycoct_xml'),cisCompt);
m5species_cis         = GlycanSpecies(glycanMLread('M5.glycoct_xml'),cisCompt);
m5gnspecies_cis       = GlycanSpecies(glycanMLread('M5gn.glycoct_xml'),cisCompt);
m4gnspecies_cis       = GlycanSpecies(glycanMLread('M4gn.glycoct_xml'),cisCompt);
m3gnspecies_cis       = GlycanSpecies(glycanMLread('M3gn.glycoct_xml'),cisCompt);
m3gn2species_cis      = GlycanSpecies(glycanMLread('M3gn2.glycoct_xml'),cisCompt);
m3gn3species_cis      = GlycanSpecies(glycanMLread('M3gn3.glycoct_xml'),cisCompt);
m3gn3bspecies_cis     = GlycanSpecies(glycanMLread('M3gn3b.glycoct_xml'),cisCompt);
m3gn4species_cis      = GlycanSpecies(glycanMLread('M3gn4.glycoct_xml'),cisCompt);
m5gngspecies_cis      = GlycanSpecies(glycanMLread('M5gng.glycoct_xml'),cisCompt);
m4gngspecies_cis      = GlycanSpecies(glycanMLread('M4gng.glycoct_xml'),cisCompt);
m3gngspecies_cis      = GlycanSpecies(glycanMLread('M3gng.glycoct_xml'),cisCompt);
m3gn2gspecies_cis     = GlycanSpecies(glycanMLread('M3gn2g.glycoct_xml'),cisCompt);
m3gn3gspecies_cis     = GlycanSpecies(glycanMLread('M3gn3g.glycoct_xml'),cisCompt);
m3gn3bgspecies_cis    = GlycanSpecies(glycanMLread('M3gn3bg.glycoct_xml'),cisCompt);
m3gn4gspecies_cis     = GlycanSpecies(glycanMLread('M3gn4g.glycoct_xml'),cisCompt);
m5gngnbspecies_cis    = GlycanSpecies(glycanMLread('M5gngnb.glycoct_xml'),cisCompt);
m4gngnbspecies_cis    = GlycanSpecies(glycanMLread('M4gngnb.glycoct_xml'),cisCompt);
m3gngnbspecies_cis    = GlycanSpecies(glycanMLread('M3gngnb.glycoct_xml'),cisCompt);
m3gn2gnbspecies_cis   = GlycanSpecies(glycanMLread('M3gn2gnb.glycoct_xml'),cisCompt);
m3gn3gnbspecies_cis   = GlycanSpecies(glycanMLread('M3gn3gnb.glycoct_xml'),cisCompt);
m3gn3bgnbspecies_cis  = GlycanSpecies(glycanMLread('M3gn3bgnb.glycoct_xml'),cisCompt);
m3gn4gnbspecies_cis   = GlycanSpecies(glycanMLread('M3gn4gnb.glycoct_xml'),cisCompt);
m5gngnbgspecies_cis   = GlycanSpecies(glycanMLread('M5gngnbg.glycoct_xml'),cisCompt);
m4gngnbgspecies_cis   = GlycanSpecies(glycanMLread('M4gngnbg.glycoct_xml'),cisCompt);
m3gngnbgspecies_cis   = GlycanSpecies(glycanMLread('M3gngnbg.glycoct_xml'),cisCompt);
m3gn2gnbgspecies_cis  = GlycanSpecies(glycanMLread('M3gn2gnbg.glycoct_xml'),cisCompt);
m3gn3gnbgspecies_cis  = GlycanSpecies(glycanMLread('M3gn3gnbg.glycoct_xml'),cisCompt);
m3gn3bgnbgspecies_cis = GlycanSpecies(glycanMLread('M3gn3bgnbg.glycoct_xml'),cisCompt);
m3gn4gnbgspecies_cis  = GlycanSpecies(glycanMLread('M3gn4gnbg.glycoct_xml'),cisCompt);

function glycanspecies = setGlySpeciesinCompt(listofstructs,thecompt)
   glycanspecies = CellArrayList;
   for i = 1 : listofstructs.length
       glycanspecies.add(GlycanSpecies(listofstructs.get(i),thecompt));
   end
end

%% create an array of species
m9species_cis         = GlycanSpecies(glycanMLread('M9.glycoct_xml'),cisCompt);
m8species_cis         = GlycanSpecies(glycanMLread('M8.glycoct_xml'),cisCompt);
m7species_cis         = GlycanSpecies(glycanMLread('M7.glycoct_xml'),cisCompt);
m6species_cis         = GlycanSpecies(glycanMLread('M6.glycoct_xml'),cisCompt);
m5species_cis         = GlycanSpecies(glycanMLread('M5.glycoct_xml'),cisCompt);
m5gnspecies_cis       = GlycanSpecies(glycanMLread('M5gn.glycoct_xml'),cisCompt);
m4gnspecies_cis       = GlycanSpecies(glycanMLread('M4gn.glycoct_xml'),cisCompt);
m3gnspecies_cis       = GlycanSpecies(glycanMLread('M3gn.glycoct_xml'),cisCompt);
m3gn2species_cis      = GlycanSpecies(glycanMLread('M3gn2.glycoct_xml'),cisCompt);
m3gn3species_cis      = GlycanSpecies(glycanMLread('M3gn3.glycoct_xml'),cisCompt);
m3gn3bspecies_cis     = GlycanSpecies(glycanMLread('M3gn3b.glycoct_xml'),cisCompt);
m3gn4species_cis      = GlycanSpecies(glycanMLread('M3gn4.glycoct_xml'),cisCompt);
m5gngspecies_cis      = GlycanSpecies(glycanMLread('M5gng.glycoct_xml'),cisCompt);
m4gngspecies_cis      = GlycanSpecies(glycanMLread('M4gng.glycoct_xml'),cisCompt);
m3gngspecies_cis      = GlycanSpecies(glycanMLread('M3gng.glycoct_xml'),cisCompt);
m3gn2gspecies_cis     = GlycanSpecies(glycanMLread('M3gn2g.glycoct_xml'),cisCompt);
m3gn3gspecies_cis     = GlycanSpecies(glycanMLread('M3gn3g.glycoct_xml'),cisCompt);
m3gn3bgspecies_cis    = GlycanSpecies(glycanMLread('M3gn3bg.glycoct_xml'),cisCompt);
m3gn4gspecies_cis     = GlycanSpecies(glycanMLread('M3gn4g.glycoct_xml'),cisCompt);
m5gngnbspecies_cis    = GlycanSpecies(glycanMLread('M5gngnb.glycoct_xml'),cisCompt);
m4gngnbspecies_cis    = GlycanSpecies(glycanMLread('M4gngnb.glycoct_xml'),cisCompt);
m3gngnbspecies_cis    = GlycanSpecies(glycanMLread('M3gngnb.glycoct_xml'),cisCompt);
m3gn2gnbspecies_cis   = GlycanSpecies(glycanMLread('M3gn2gnb.glycoct_xml'),cisCompt);
m3gn3gnbspecies_cis   = GlycanSpecies(glycanMLread('M3gn3gnb.glycoct_xml'),cisCompt);
m3gn3bgnbspecies_cis  = GlycanSpecies(glycanMLread('M3gn3bgnb.glycoct_xml'),cisCompt);
m3gn4gnbspecies_cis   = GlycanSpecies(glycanMLread('M3gn4gnb.glycoct_xml'),cisCompt);
m5gngnbgspecies_cis   = GlycanSpecies(glycanMLread('M5gngnbg.glycoct_xml'),cisCompt);
m4gngnbgspecies_cis   = GlycanSpecies(glycanMLread('M4gngnbg.glycoct_xml'),cisCompt);
m3gngnbgspecies_cis   = GlycanSpecies(glycanMLread('M3gngnbg.glycoct_xml'),cisCompt);
m3gn2gnbgspecies_cis  = GlycanSpecies(glycanMLread('M3gn2gnbg.glycoct_xml'),cisCompt);
m3gn3gnbgspecies_cis  = GlycanSpecies(glycanMLread('M3gn3gnbg.glycoct_xml'),cisCompt);
m3gn3bgnbgspecies_cis = GlycanSpecies(glycanMLread('M3gn3bgnbg.glycoct_xml'),cisCompt);
m3gn4gnbgspecies_cis  = GlycanSpecies(glycanMLread('M3gn4gnbg.glycoct_xml'),cisCompt);

m9species_media         = GlycanSpecies(glycanMLread('M9.glycoct_xml'),mediaCompt);
m8species_media         = GlycanSpecies(glycanMLread('M8.glycoct_xml'),mediaCompt);
m7species_media         = GlycanSpecies(glycanMLread('M7.glycoct_xml'),mediaCompt);
m6species_media         = GlycanSpecies(glycanMLread('M6.glycoct_xml'),mediaCompt);
m5species_media         = GlycanSpecies(glycanMLread('M5.glycoct_xml'),mediaCompt);
m5gnspecies_media       = GlycanSpecies(glycanMLread('M5gn.glycoct_xml'),mediaCompt);
m4gnspecies_media       = GlycanSpecies(glycanMLread('M4gn.glycoct_xml'),mediaCompt);
m3gnspecies_media       = GlycanSpecies(glycanMLread('M3gn.glycoct_xml'),mediaCompt);
m3gn2species_media      = GlycanSpecies(glycanMLread('M3gn2.glycoct_xml'),mediaCompt);
m3gn3species_media      = GlycanSpecies(glycanMLread('M3gn3.glycoct_xml'),mediaCompt);
m3gn3bspecies_media     = GlycanSpecies(glycanMLread('M3gn3b.glycoct_xml'),mediaCompt);
m3gn4species_media      = GlycanSpecies(glycanMLread('M3gn4.glycoct_xml'),mediaCompt);
m5gngspecies_media      = GlycanSpecies(glycanMLread('M5gng.glycoct_xml'),mediaCompt);
m4gngspecies_media      = GlycanSpecies(glycanMLread('M4gng.glycoct_xml'),mediaCompt);
m3gngspecies_media      = GlycanSpecies(glycanMLread('M3gng.glycoct_xml'),mediaCompt);
m3gn2gspecies_media     = GlycanSpecies(glycanMLread('M3gn2g.glycoct_xml'),mediaCompt);
m3gn3gspecies_media     = GlycanSpecies(glycanMLread('M3gn3g.glycoct_xml'),mediaCompt);
m3gn3bgspecies_media    = GlycanSpecies(glycanMLread('M3gn3bg.glycoct_xml'),mediaCompt);
m3gn4gspecies_media     = GlycanSpecies(glycanMLread('M3gn4g.glycoct_xml'),mediaCompt);
m5gngnbspecies_media    = GlycanSpecies(glycanMLread('M5gngnb.glycoct_xml'),mediaCompt);
m4gngnbspecies_media    = GlycanSpecies(glycanMLread('M4gngnb.glycoct_xml'),mediaCompt);
m3gngnbspecies_media    = GlycanSpecies(glycanMLread('M3gngnb.glycoct_xml'),mediaCompt);
m3gn2gnbspecies_media   = GlycanSpecies(glycanMLread('M3gn2gnb.glycoct_xml'),mediaCompt);
m3gn3gnbspecies_media   = GlycanSpecies(glycanMLread('M3gn3gnb.glycoct_xml'),mediaCompt);
m3gn3bgnbspecies_media  = GlycanSpecies(glycanMLread('M3gn3bgnb.glycoct_xml'),mediaCompt);
m3gn4gnbspecies_media   = GlycanSpecies(glycanMLread('M3gn4gnb.glycoct_xml'),mediaCompt);
m5gngnbgspecies_media   = GlycanSpecies(glycanMLread('M5gngnbg.glycoct_xml'),mediaCompt);
m4gngnbgspecies_media   = GlycanSpecies(glycanMLread('M4gngnbg.glycoct_xml'),mediaCompt);
m3gngnbgspecies_media   = GlycanSpecies(glycanMLread('M3gngnbg.glycoct_xml'),mediaCompt);
m3gn2gnbgspecies_media  = GlycanSpecies(glycanMLread('M3gn2gnbg.glycoct_xml'),mediaCompt);
m3gn3gnbgspecies_media  = GlycanSpecies(glycanMLread('M3gn3gnbg.glycoct_xml'),mediaCompt);
m3gn3bgnbgspecies_media = GlycanSpecies(glycanMLread('M3gn3bgnbg.glycoct_xml'),mediaCompt);
m3gn4gnbgspecies_media  = GlycanSpecies(glycanMLread('M3gn4gnbg.glycoct_xml'),mediaCompt);

m9species_trans         = GlycanSpecies(glycanMLread('M9.glycoct_xml'),transCompt);
m8species_trans         = GlycanSpecies(glycanMLread('M8.glycoct_xml'),transCompt);
m7species_trans         = GlycanSpecies(glycanMLread('M7.glycoct_xml'),transCompt);
m6species_trans         = GlycanSpecies(glycanMLread('M6.glycoct_xml'),transCompt);
m5species_trans         = GlycanSpecies(glycanMLread('M5.glycoct_xml'),transCompt);
m5gnspecies_trans       = GlycanSpecies(glycanMLread('M5gn.glycoct_xml'),transCompt);
m4gnspecies_trans       = GlycanSpecies(glycanMLread('M4gn.glycoct_xml'),transCompt);
m3gnspecies_trans       = GlycanSpecies(glycanMLread('M3gn.glycoct_xml'),transCompt);
m3gn2species_trans      = GlycanSpecies(glycanMLread('M3gn2.glycoct_xml'),transCompt);
m3gn3species_trans      = GlycanSpecies(glycanMLread('M3gn3.glycoct_xml'),transCompt);
m3gn3bspecies_trans     = GlycanSpecies(glycanMLread('M3gn3b.glycoct_xml'),transCompt);
m3gn4species_trans      = GlycanSpecies(glycanMLread('M3gn4.glycoct_xml'),transCompt);
m5gngspecies_trans      = GlycanSpecies(glycanMLread('M5gng.glycoct_xml'),transCompt);
m4gngspecies_trans      = GlycanSpecies(glycanMLread('M4gng.glycoct_xml'),transCompt);
m3gngspecies_trans      = GlycanSpecies(glycanMLread('M3gng.glycoct_xml'),transCompt);
m3gn2gspecies_trans     = GlycanSpecies(glycanMLread('M3gn2g.glycoct_xml'),transCompt);
m3gn3gspecies_trans     = GlycanSpecies(glycanMLread('M3gn3g.glycoct_xml'),transCompt);
m3gn3bgspecies_trans    = GlycanSpecies(glycanMLread('M3gn3bg.glycoct_xml'),transCompt);
m3gn4gspecies_trans     = GlycanSpecies(glycanMLread('M3gn4g.glycoct_xml'),transCompt);
m5gngnbspecies_trans    = GlycanSpecies(glycanMLread('M5gngnb.glycoct_xml'),transCompt);
m4gngnbspecies_trans    = GlycanSpecies(glycanMLread('M4gngnb.glycoct_xml'),transCompt);
m3gngnbspecies_trans    = GlycanSpecies(glycanMLread('M3gngnb.glycoct_xml'),transCompt);
m3gn2gnbspecies_trans   = GlycanSpecies(glycanMLread('M3gn2gnb.glycoct_xml'),transCompt);
m3gn3gnbspecies_trans   = GlycanSpecies(glycanMLread('M3gn3gnb.glycoct_xml'),transCompt);
m3gn3bgnbspecies_trans  = GlycanSpecies(glycanMLread('M3gn3bgnb.glycoct_xml'),transCompt);
m3gn4gnbspecies_trans   = GlycanSpecies(glycanMLread('M3gn4gnb.glycoct_xml'),transCompt);
m5gngnbgspecies_trans   = GlycanSpecies(glycanMLread('M5gngnbg.glycoct_xml'),transCompt);
m4gngnbgspecies_trans   = GlycanSpecies(glycanMLread('M4gngnbg.glycoct_xml'),transCompt);
m3gngnbgspecies_trans   = GlycanSpecies(glycanMLread('M3gngnbg.glycoct_xml'),transCompt);
m3gn2gnbgspecies_trans  = GlycanSpecies(glycanMLread('M3gn2gnbg.glycoct_xml'),transCompt);
m3gn3gnbgspecies_trans  = GlycanSpecies(glycanMLread('M3gn3gnbg.glycoct_xml'),transCompt);
m3gn3bgnbgspecies_trans = GlycanSpecies(glycanMLread('M3gn3bgnbg.glycoct_xml'),transCompt);
m3gn4gnbgspecies_trans  = GlycanSpecies(glycanMLread('M3gn4gnbg.glycoct_xml'),transCompt);

m9species_tgn         = GlycanSpecies(glycanMLread('M9.glycoct_xml'),tgnCompt);
m8species_tgn         = GlycanSpecies(glycanMLread('M8.glycoct_xml'),tgnCompt);
m7species_tgn         = GlycanSpecies(glycanMLread('M7.glycoct_xml'),tgnCompt);
m6species_tgn         = GlycanSpecies(glycanMLread('M6.glycoct_xml'),tgnCompt);
m5species_tgn         = GlycanSpecies(glycanMLread('M5.glycoct_xml'),tgnCompt);
m5gnspecies_tgn       = GlycanSpecies(glycanMLread('M5gn.glycoct_xml'),tgnCompt);
m4gnspecies_tgn       = GlycanSpecies(glycanMLread('M4gn.glycoct_xml'),tgnCompt);
m3gnspecies_tgn       = GlycanSpecies(glycanMLread('M3gn.glycoct_xml'),tgnCompt);
m3gn2species_tgn      = GlycanSpecies(glycanMLread('M3gn2.glycoct_xml'),tgnCompt);
m3gn3species_tgn      = GlycanSpecies(glycanMLread('M3gn3.glycoct_xml'),tgnCompt);
m3gn3bspecies_tgn     = GlycanSpecies(glycanMLread('M3gn3b.glycoct_xml'),tgnCompt);
m3gn4species_tgn      = GlycanSpecies(glycanMLread('M3gn4.glycoct_xml'),tgnCompt);
m5gngspecies_tgn      = GlycanSpecies(glycanMLread('M5gng.glycoct_xml'),tgnCompt);
m4gngspecies_tgn      = GlycanSpecies(glycanMLread('M4gng.glycoct_xml'),tgnCompt);
m3gngspecies_tgn      = GlycanSpecies(glycanMLread('M3gng.glycoct_xml'),tgnCompt);
m3gn2gspecies_tgn     = GlycanSpecies(glycanMLread('M3gn2g.glycoct_xml'),tgnCompt);
m3gn3gspecies_tgn     = GlycanSpecies(glycanMLread('M3gn3g.glycoct_xml'),tgnCompt);
m3gn3bgspecies_tgn    = GlycanSpecies(glycanMLread('M3gn3bg.glycoct_xml'),tgnCompt);
m3gn4gspecies_tgn     = GlycanSpecies(glycanMLread('M3gn4g.glycoct_xml'),tgnCompt);
m5gngnbspecies_tgn    = GlycanSpecies(glycanMLread('M5gngnb.glycoct_xml'),tgnCompt);
m4gngnbspecies_tgn    = GlycanSpecies(glycanMLread('M4gngnb.glycoct_xml'),tgnCompt);
m3gngnbspecies_tgn    = GlycanSpecies(glycanMLread('M3gngnb.glycoct_xml'),tgnCompt);
m3gn2gnbspecies_tgn   = GlycanSpecies(glycanMLread('M3gn2gnb.glycoct_xml'),tgnCompt);
m3gn3gnbspecies_tgn   = GlycanSpecies(glycanMLread('M3gn3gnb.glycoct_xml'),tgnCompt);
m3gn3bgnbspecies_tgn  = GlycanSpecies(glycanMLread('M3gn3bgnb.glycoct_xml'),tgnCompt);
m3gn4gnbspecies_tgn   = GlycanSpecies(glycanMLread('M3gn4gnb.glycoct_xml'),tgnCompt);
m5gngnbgspecies_tgn   = GlycanSpecies(glycanMLread('M5gngnbg.glycoct_xml'),tgnCompt);
m4gngnbgspecies_tgn   = GlycanSpecies(glycanMLread('M4gngnbg.glycoct_xml'),tgnCompt);
m3gngnbgspecies_tgn   = GlycanSpecies(glycanMLread('M3gngnbg.glycoct_xml'),tgnCompt);
m3gn2gnbgspecies_tgn  = GlycanSpecies(glycanMLread('M3gn2gnbg.glycoct_xml'),tgnCompt);
m3gn3gnbgspecies_tgn  = GlycanSpecies(glycanMLread('M3gn3gnbg.glycoct_xml'),tgnCompt);
m3gn3bgnbgspecies_tgn = GlycanSpecies(glycanMLread('M3gn3bgnbg.glycoct_xml'),tgnCompt);
m3gn4gnbgspecies_tgn  = GlycanSpecies(glycanMLread('M3gn4gnbg.glycoct_xml'),tgnCompt);

allglycans  = CellArrayList;
allglycans.add(m9species_cis);
allglycans.add(m8species_cis);
allglycans.add(m7species_cis);
allglycans.add(m6species_cis);
allglycans.add(m5species_cis);
allglycans.add(m5gnspecies_cis);
allglycans.add(m4gnspecies_cis);
allglycans.add(m3gnspecies_cis);
allglycans.add(m3gn2species_cis);
allglycans.add(m3gn3species_cis);
allglycans.add(m3gn3bspecies_cis);
allglycans.add(m3gn4species_cis);
allglycans.add(m5gngspecies_cis);
allglycans.add(m4gngspecies_cis);
allglycans.add(m3gngspecies_cis);
allglycans.add(m3gn2gspecies_cis);
allglycans.add(m3gn3gspecies_cis);
allglycans.add(m3gn3bgspecies_cis);
allglycans.add(m3gn4gspecies_cis);
allglycans.add(m5gngnbspecies_cis);
allglycans.add(m4gngnbspecies_cis);
allglycans.add(m3gngnbspecies_cis);
allglycans.add(m3gn2gnbspecies_cis);
allglycans.add(m3gn3gnbspecies_cis);
allglycans.add(m3gn3bgnbspecies_cis);
allglycans.add(m3gn4gnbspecies_cis);
allglycans.add(m5gngnbgspecies_cis);
allglycans.add(m4gngnbgspecies_cis);
allglycans.add(m3gngnbgspecies_cis);
allglycans.add(m3gn2gnbgspecies_cis);
allglycans.add(m3gn3gnbgspecies_cis);
allglycans.add(m3gn3bgnbgspecies_cis);
allglycans.add(m3gn4gnbgspecies_cis);
allglycans.add(m9species_media);
allglycans.add(m8species_media);
allglycans.add(m7species_media);
allglycans.add(m6species_media);
allglycans.add(m5species_media);
allglycans.add(m5gnspecies_media);
allglycans.add(m4gnspecies_media);
allglycans.add(m3gnspecies_media);
allglycans.add(m3gn2species_media);
allglycans.add(m3gn3species_media);
allglycans.add(m3gn3bspecies_media);
allglycans.add(m3gn4species_media);
allglycans.add(m5gngspecies_media);
allglycans.add(m4gngspecies_media);
allglycans.add(m3gngspecies_media);
allglycans.add(m3gn2gspecies_media);
allglycans.add(m3gn3gspecies_media);
allglycans.add(m3gn3bgspecies_media);
allglycans.add(m3gn4gspecies_media);
allglycans.add(m5gngnbspecies_media);
allglycans.add(m4gngnbspecies_media);
allglycans.add(m3gngnbspecies_media);
allglycans.add(m3gn2gnbspecies_media);
allglycans.add(m3gn3gnbspecies_media);
allglycans.add(m3gn3bgnbspecies_media);
allglycans.add(m3gn4gnbspecies_media);
allglycans.add(m5gngnbgspecies_media);
allglycans.add(m4gngnbgspecies_media);
allglycans.add(m3gngnbgspecies_media);
allglycans.add(m3gn2gnbgspecies_media);
allglycans.add(m3gn3gnbgspecies_media);
allglycans.add(m3gn3bgnbgspecies_media);
allglycans.add(m3gn4gnbgspecies_media);
allglycans.add(m9species_trans);
allglycans.add(m8species_trans);
allglycans.add(m7species_trans);
allglycans.add(m6species_trans);
allglycans.add(m5species_trans);
allglycans.add(m5gnspecies_trans);
allglycans.add(m4gnspecies_trans);
allglycans.add(m3gnspecies_trans);
allglycans.add(m3gn2species_trans);
allglycans.add(m3gn3species_trans);
allglycans.add(m3gn3bspecies_trans);
allglycans.add(m3gn4species_trans);
allglycans.add(m5gngspecies_trans);
allglycans.add(m4gngspecies_trans);
allglycans.add(m3gngspecies_trans);
allglycans.add(m3gn2gspecies_trans);
allglycans.add(m3gn3gspecies_trans);
allglycans.add(m3gn3bgspecies_trans);
allglycans.add(m3gn4gspecies_trans);
allglycans.add(m5gngnbspecies_trans);
allglycans.add(m4gngnbspecies_trans);
allglycans.add(m3gngnbspecies_trans);
allglycans.add(m3gn2gnbspecies_trans);
allglycans.add(m3gn3gnbspecies_trans);
allglycans.add(m3gn3bgnbspecies_trans);
allglycans.add(m3gn4gnbspecies_trans);
allglycans.add(m5gngnbgspecies_trans);
allglycans.add(m4gngnbgspecies_trans);
allglycans.add(m3gngnbgspecies_trans);
allglycans.add(m3gn2gnbgspecies_trans);
allglycans.add(m3gn3gnbgspecies_trans);
allglycans.add(m3gn3bgnbgspecies_trans);
allglycans.add(m3gn4gnbgspecies_trans);
allglycans.add(m9species_tgn);
allglycans.add(m8species_tgn);
allglycans.add(m7species_tgn);
allglycans.add(m6species_tgn);
allglycans.add(m5species_tgn);
allglycans.add(m5gnspecies_tgn);
allglycans.add(m4gnspecies_tgn);
allglycans.add(m3gnspecies_tgn);
allglycans.add(m3gn2species_tgn);
allglycans.add(m3gn3species_tgn);
allglycans.add(m3gn3bspecies_tgn);
allglycans.add(m3gn4species_tgn);
allglycans.add(m5gngspecies_tgn);
allglycans.add(m4gngspecies_tgn);
allglycans.add(m3gngspecies_tgn);
allglycans.add(m3gn2gspecies_tgn);
allglycans.add(m3gn3gspecies_tgn);
allglycans.add(m3gn3bgspecies_tgn);
allglycans.add(m3gn4gspecies_tgn);
allglycans.add(m5gngnbspecies_tgn);
allglycans.add(m4gngnbspecies_tgn);
allglycans.add(m3gngnbspecies_tgn);
allglycans.add(m3gn2gnbspecies_tgn);
allglycans.add(m3gn3gnbspecies_tgn);
allglycans.add(m3gn3bgnbspecies_tgn);
allglycans.add(m3gn4gnbspecies_tgn);
allglycans.add(m5gngnbgspecies_tgn);
allglycans.add(m4gngnbgspecies_tgn);
allglycans.add(m3gngnbgspecies_tgn);
allglycans.add(m3gn2gnbgspecies_tgn);
allglycans.add(m3gn3gnbgspecies_tgn);
allglycans.add(m3gn3bgnbgspecies_tgn);
allglycans.add(m3gn4gnbgspecies_tgn);

% setSpeciesid
speciesname ={'m9species_cis','m8species_cis','m7species_cis','m6species_cis',...
              'm5species_cis','m5gnspecies_cis','m4gnspecies_cis','m3gnspecies_cis',...
              'm3gn2species_cis','m3gn3species_cis','m3gn3bspecies_cis','m3gn4species_cis',...
              'm5gngspecies_cis','m4gngspecies_cis','m3gngspecies_cis','m3gn2gspecies_cis',...
              'm3gn3gspecies_cis','m3gn3bgspecies_cis','m3gn4bgspecies_cis','m5gngnbspecies_cis',...
              'm4gngnbspecies_cis','m3gngnbspecies_cis','m3gn2gnbspecies_cis','m3gn3gnbspecies_cis',...
              'm3gn3bgnbspecies_cis','m3gn4gnbspecies_cis','m5gngnbgspecies_cis','m4gngnbgspecies_cis',...
              'm3gngnbgspecies_cis','m3gn2gnbgspecies_cis','m3gn3gnbgspecies_cis','m3gn3bgnbgspecies_cis',...
              'm3gn4gnbgspecies_cis',...
              'm9species_media','m8species_media','m7species_media','m6species_media'...
              'm5species_media','m5gnspecies_media','m4gnspecies_media','m3gnspecies_media',...
              'm3gn2species_media','m3gn3species_media','m3gn3bspecies_media','m3gn4species_media',...
              'm5gngspecies_media','m4gngspecies_media','m3gngspecies_media','m3gn2gspecies_media',...
              'm3gn3gspecies_media','m3gn3bgspecies_media','m3gn4bgspecies_media','m5gngnbspecies_media',...
              'm4gngnbspecies_media','m3gngnbspecies_media','m3gn2gnbspecies_media','m3gn3gnbspecies_media',...
              'm3gn3bgnbspecies_media','m3gn4gnbspecies_media','m5gngnbgspecies_media','m4gngnbgspecies_media',...
              'm3gngnbgspecies_media','m3gn2gnbgspecies_media','m3gn3gnbgspecies_media','m3gn3bgnbgspecies_media',...
              'm3gn4gnbgspecies_media',...
              'm9species_trans','m8species_trans','m7species_trans','m6species_trans',...
              'm5species_trans','m5gnspecies_trans','m4gnspecies_trans','m3gnspecies_trans',...
              'm3gn2species_trans','m3gn3species_trans','m3gn3bspecies_trans','m3gn4species_trans',...
              'm5gngspecies_trans','m4gngspecies_trans','m3gngspecies_trans','m3gn2gspecies_trans',...
              'm3gn3gspecies_trans','m3gn3bgspecies_trans','m3gn4bgspecies_trans','m5gngnbspecies_trans',...
              'm4gngnbspecies_trans','m3gngnbspecies_trans','m3gn2gnbspecies_trans','m3gn3gnbspecies_trans',...
              'm3gn3bgnbspecies_trans','m3gn4gnbspecies_trans','m5gngnbgspecies_trans','m4gngnbgspecies_trans',...
              'm3gngnbgspecies_trans','m3gn2gnbgspecies_trans','m3gn3gnbgspecies_trans','m3gn3bgnbgspecies_trans',...
              'm3gn4gnbgspecies_trans',...
              'm9species_tgn','m8species_tgn','m7species_tgn','m6species_tgn',...
              'm5species_tgn','m5gnspecies_tgn','m4gnspecies_tgn','m3gnspecies_tgn',...
              'm3gn2species_tgn','m3gn3species_tgn','m3gn3bspecies_tgn','m3gn4species_tgn',...
              'm5gngspecies_tgn','m4gngspecies_tgn','m3gngspecies_tgn','m3gn2gspecies_tgn',...
              'm3gn3gspecies_tgn','m3gn3bgspecies_tgn','m3gn4bgspecies_tgn','m5gngnbspecies_tgn',...
              'm4gngnbspecies_tgn','m3gngnbspecies_tgn','m3gn2gnbspecies_tgn','m3gn3gnbspecies_tgn',...
              'm3gn3bgnbspecies_tgn','m3gn4gnbspecies_tgn','m5gngnbgspecies_tgn','m4gngnbgspecies_tgn',...
              'm3gngnbgspecies_tgn','m3gn2gnbgspecies_tgn','m3gn3gnbgspecies_tgn','m3gn3bgnbgspecies_tgn',...
              'm3gn4gnbgspecies_tgn'};

for i=1:length(allglycans)
    allglycans.get(i).id = speciesname{i};
end

for i=1:length(allglycans)
    allglycans.get(i).initConc   = 0;
    allglycans.get(i).initAmount = 0;
end

allglycans.get(1).initConc   = 1;
allglycans.get(1).initAmount = 1;
allglycans.get(34).initConc   = 1;
allglycans.get(34).initAmount = 1;
allglycans.get(67).initConc   = 1;
allglycans.get(67).initAmount = 1;
allglycans.get(100).initConc   = 1;
allglycans.get(100).initAmount = 1;


%% create an array of reactions
allrxns  = CellArrayList;

% load enzyme database
enzdbmatfilename   = 'glyenzDB.mat';
enzdb              = enzdbmatLoad(enzdbmatfilename);
mgat1  = enzdb('mgat1');
mgat2  = enzdb('mgat2');
mgat3  = enzdb('mgat3');
mgat4  = enzdb('mgat4');
mgat5  = enzdb('mgat5');
mani   = enzdb('mani');
galt   = enzdb('galt');
manii  = enzdb('manii');

rt1     = struct('Vm',485.96,'Km',0.70);
enzkineticsobj = MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt1);
mgat1.enzkinetics = enzkineticsobj;

rt2     = struct('Vm',151.19,'Km',0.51);
enzkineticsobj= MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt2);
mgat2.enzkinetics = enzkineticsobj;

rt3     = struct('Vm',4319.65,'Km',0.51);
enzkineticsobj = MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt3);
specistruct= CellArrayList;
specistruct.add(m5gnspecies_cis.glycanStruct);
specistruct.add(m4gnspecies_cis.glycanStruct);
specistruct.add(m3gnspecies_cis.glycanStruct);
specifiratio = CellArrayList;
ratio1 = struct('Vm',1,'Km',21.053);
ratio2 = struct('Vm',1,'Km',21.053);
ratio3 = struct('Vm',1,'Km',21.053);
specifiratio.add(ratio1);
specifiratio.add(ratio2);
specifiratio.add(ratio3);
setSubstrateSpecificityK(enzkineticsobj,specistruct,specifiratio,1);
mgat3.enzkinetics = enzkineticsobj;

rt4     = struct('Vm',10.80,'Km',9.18);
enzkineticsobj = MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt4);
mgat4.enzkinetics = enzkineticsobj;

rt5     = struct('Vm',10.80,'Km',0.24);
enzkineticsobj = MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt5);
specistruct= CellArrayList;
specistruct.add(m3gn3species_cis.glycanStruct);
specifiratio = CellArrayList;
ratio1 = struct('Vm',1,'Km',0.692);
specifiratio.add(ratio1);
setSubstrateSpecificityK(enzkineticsobj,specistruct,specifiratio,1);
mgat5.enzkinetics = enzkineticsobj;

rt6     = struct('Vm',323.97,'Km',0.27);
enzkineticsobj6 = MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt6);
mani.enzkinetics = enzkineticsobj6;

rt7     = struct('Vm',626.35,'Km',0.35);
enzkineticsobj = MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt7);
specistruct= CellArrayList;
specistruct.add(m5gnspecies_cis.glycanStruct);
specistruct.add(m4gnspecies_cis.glycanStruct);
specistruct.add(m3gnspecies_cis.glycanStruct);
specistruct.add(m3gn3bspecies_cis.glycanStruct);
specistruct.add(m3gn3species_cis.glycanStruct);
specistruct.add(m3gn4species_cis.glycanStruct);
specistruct.add(m5gngnbspecies_cis.glycanStruct);
specistruct.add(m4gngnbspecies_cis.glycanStruct);
specistruct.add(m3gngnbspecies_cis.glycanStruct);
specistruct.add(m3gn2gnbspecies_cis.glycanStruct);
specistruct.add(m3gn3bgnbspecies_cis.glycanStruct);
specistruct.add(m3gn3gnbspecies_cis.glycanStruct);
specistruct.add(m3gn4gnbspecies_cis.glycanStruct);
specifiratio = CellArrayList;
ratio1 = struct('Vm',1,'Km',30.769);
ratio2 = struct('Vm',1,'Km',30.769);
ratio3 = struct('Vm',1,'Km',30.769);
ratio4 = struct('Vm',1,'Km',0.538);
ratio5 = struct('Vm',1,'Km',0.385);
ratio6 = struct('Vm',1,'Km',0.308);
ratio7 = struct('Vm',1,'Km',30.769);
ratio8 = struct('Vm',1,'Km',30.769);
ratio9 = struct('Vm',1,'Km',30.769);
ratio10 = struct('Vm',1,'Km',3.846);
ratio11 = struct('Vm',1,'Km',1.692);
ratio12 = struct('Vm',1,'Km',1.538);
ratio13 = struct('Vm',1,'Km',1.077);
specifiratio.add(ratio1);
specifiratio.add(ratio2);
specifiratio.add(ratio3);
specifiratio.add(ratio4);
specifiratio.add(ratio5);
specifiratio.add(ratio6);
specifiratio.add(ratio7);
specifiratio.add(ratio8);
specifiratio.add(ratio9);
specifiratio.add(ratio10);
specifiratio.add(ratio11);
specifiratio.add(ratio12);
specifiratio.add(ratio13);
setSubstrateSpecificityK(enzkineticsobj,specistruct,specifiratio,1);
galt.enzkinetics = enzkineticsobj;

rt8            = struct('Vm',323.97,'Km',0.27);
enzkineticsobj = MMenKinetics(MMenKinetics.Comp,MMenKinetics.Ensm,rt8);
specistruct    = CellArrayList;
specistruct.add(m5gnspecies_cis.glycanStruct);
specifiratio   = CellArrayList;
ratio1         = struct('Vm',1,'Km',2);
specifiratio.add(ratio1);
enzkineticsobj.setSubstrateSpecificityK(specistruct,specifiratio,1);
manii.enzkinetics = enzkineticsobj;

%set up rxn
m8rxn_cis   = Rxn(m9species_cis,m8species_cis,mani);
m8rxn_media   = Rxn(m9species_media,m8species_media,mani);
m8rxn_trans   = Rxn(m9species_trans,m8species_trans,mani);
m8rxn_tgn   = Rxn(m9species_tgn,m8species_tgn,mani);
%m8rxn.setrxnkinetics;
allrxns.add(m8rxn_cis);
allrxns.add(m8rxn_media);
allrxns.add(m8rxn_trans);
allrxns.add(m8rxn_tgn);

m7rxn_cis   = Rxn(m8species_cis,m7species_cis,mani);
m7rxn_media   = Rxn(m8species_media,m7species_media,mani);
m7rxn_trans   = Rxn(m8species_trans,m7species_trans,mani);
m7rxn_tgn   = Rxn(m8species_tgn,m7species_tgn,mani);
%m7rxn.setrxnkinetics;
allrxns.add(m7rxn_cis);
allrxns.add(m7rxn_media );
allrxns.add(m7rxn_trans);
allrxns.add(m7rxn_tgn);

m6rxn_cis   = Rxn(m7species_cis,m6species_cis,mani);
m6rxn_media   = Rxn(m7species_media,m6species_media,mani);
m6rxn_trans   = Rxn(m7species_trans,m6species_trans,mani);
m6rxn_tgn   = Rxn(m7species_tgn,m6species_tgn,mani);
%m6rxn.setrxnkinetics;
allrxns.add(m6rxn_cis);
allrxns.add(m6rxn_media);
allrxns.add(m6rxn_trans);
allrxns.add(m6rxn_tgn);

m5rxn_cis   = Rxn(m6species_cis,m5species_cis,mani);
m5rxn_media   = Rxn(m6species_media,m5species_media,mani);
m5rxn_trans   = Rxn(m6species_trans,m5species_trans,mani);
m5rxn_tgn   = Rxn(m6species_tgn,m5species_tgn,mani);
%m5rxn.setrxnkinetics;
allrxns.add(m5rxn_cis);
allrxns.add(m5rxn_media );
allrxns.add(m5rxn_trans);
allrxns.add(m5rxn_tgn);

m5gnrxn_cis = Rxn(m5species_cis,m5gnspecies_cis,mgat1);
m5gnrxn_media = Rxn(m5species_media,m5gnspecies_media,mgat1);
m5gnrxn_trans = Rxn(m5species_trans,m5gnspecies_trans,mgat1);
m5gnrxn_tgn = Rxn(m5species_tgn,m5gnspecies_tgn,mgat1);
%m5gnrxn.setrxnkinetics;
allrxns.add(m5gnrxn_cis);
allrxns.add(m5gnrxn_media);
allrxns.add(m5gnrxn_trans);
allrxns.add(m5gnrxn_tgn);

m4gnrxn_cis = Rxn(m5gnspecies_cis,m4gnspecies_cis,manii);
m4gnrxn_media = Rxn(m5gnspecies_media,m4gnspecies_media,manii);
m4gnrxn_trans = Rxn(m5gnspecies_trans,m4gnspecies_trans,manii);
m4gnrxn_tgn = Rxn(m5gnspecies_tgn,m4gnspecies_tgn,manii);
%m4gnrxn.setrxnkinetics;
allrxns.add(m4gnrxn_cis);
allrxns.add(m4gnrxn_media);
allrxns.add(m4gnrxn_trans);
allrxns.add(m4gnrxn_tgn);

m3gnrxn_cis = Rxn(m4gnspecies_cis,m3gnspecies_cis,manii);
m3gnrxn_media = Rxn(m4gnspecies_media,m3gnspecies_media,manii); 
m3gnrxn_trans = Rxn(m4gnspecies_trans,m3gnspecies_trans,manii); 
m3gnrxn_tgn = Rxn(m4gnspecies_tgn,m3gnspecies_tgn,manii); 
%m3gnrxn.setrxnkinetics;
allrxns.add(m3gnrxn_cis);
allrxns.add(m3gnrxn_media);
allrxns.add(m3gnrxn_trans);
allrxns.add(m3gnrxn_tgn);

m3gn2rxn_cis = Rxn(m3gnspecies_cis,m3gn2species_cis,mgat2);
m3gn2rxn_media = Rxn(m3gnspecies_media,m3gn2species_media,mgat2);
m3gn2rxn_trans = Rxn(m3gnspecies_trans,m3gn2species_trans,mgat2);
m3gn2rxn_tgn = Rxn(m3gnspecies_tgn,m3gn2species_tgn,mgat2);
%m3gn2rxn.setrxnkinetics;
allrxns.add(m3gn2rxn_cis);
allrxns.add(m3gn2rxn_media);
allrxns.add(m3gn2rxn_trans);
allrxns.add(m3gn2rxn_tgn);

m3gn3brxn_cis = Rxn(m3gn2species_cis,m3gn3bspecies_cis,mgat5);
m3gn3brxn_media = Rxn(m3gn2species_media,m3gn3bspecies_media,mgat5);
m3gn3brxn_trans = Rxn(m3gn2species_trans,m3gn3bspecies_trans,mgat5);
m3gn3brxn_tgn = Rxn(m3gn2species_tgn,m3gn3bspecies_tgn,mgat5);
%m3gn3brxn.setrxnkinetics;
allrxns.add(m3gn3brxn_cis);
allrxns.add(m3gn3brxn_media);
allrxns.add(m3gn3brxn_trans);
allrxns.add(m3gn3brxn_tgn);

m3gn3rxn_cis  = Rxn(m3gn2species_cis,m3gn3species_cis,mgat4);
m3gn3rxn_media  = Rxn(m3gn2species_media,m3gn3species_media,mgat4);
m3gn3rxn_trans  = Rxn(m3gn2species_trans,m3gn3species_trans,mgat4);
m3gn3rxn_tgn  = Rxn(m3gn2species_tgn,m3gn3species_tgn,mgat4);
%m3gn3rxn.setrxnkinetics;
allrxns.add(m3gn3rxn_cis);
allrxns.add(m3gn3rxn_media);
allrxns.add(m3gn3rxn_trans);
allrxns.add(m3gn3rxn_tgn);

m3gn4arxn_cis = Rxn(m3gn3bspecies_cis,m3gn4species_cis,mgat4);
m3gn4arxn_media = Rxn(m3gn3bspecies_media,m3gn4species_media,mgat4);
m3gn4arxn_trans = Rxn(m3gn3bspecies_trans,m3gn4species_trans,mgat4);
m3gn4arxn_tgn = Rxn(m3gn3bspecies_tgn,m3gn4species_tgn,mgat4);
%m3gn4arxn.setrxnkinetics;
allrxns.add(m3gn4arxn_cis);
allrxns.add(m3gn4arxn_media);
allrxns.add(m3gn4arxn_trans);
allrxns.add(m3gn4arxn_tgn);

m3gn4brxn_cis = Rxn(m3gn3species_cis,m3gn4species_cis,mgat5);
m3gn4brxn_media = Rxn(m3gn3species_media,m3gn4species_media,mgat5);
m3gn4brxn_trans = Rxn(m3gn3species_trans,m3gn4species_trans,mgat5);
m3gn4brxn_tgn = Rxn(m3gn3species_tgn,m3gn4species_tgn,mgat5);
%m3gn4brxn.setrxnkinetics;
allrxns.add(m3gn4brxn_cis);
allrxns.add(m3gn4brxn_media);
allrxns.add(m3gn4brxn_trans);
allrxns.add(m3gn4brxn_tgn);

m5gngrxn_cis = Rxn(m5gnspecies_cis ,m5gngspecies_cis,galt);
m5gngrxn_media = Rxn(m5gnspecies_media,m5gngspecies_media,galt);
m5gngrxn_trans = Rxn(m5gnspecies_trans,m5gngspecies_trans,galt);
m5gngrxn_tgn = Rxn(m5gnspecies_tgn,m5gngspecies_tgn,galt);
%m5gngrxn.setrxnkinetics;
allrxns.add(m5gngrxn_cis);
allrxns.add(m5gngrxn_media);
allrxns.add(m5gngrxn_trans);
allrxns.add(m5gngrxn_tgn);

m4gngrxn_cis = Rxn(m4gnspecies_cis,m4gngspecies_cis,galt);
m4gngrxn_media = Rxn(m4gnspecies_media,m4gngspecies_media,galt);
m4gngrxn_trans = Rxn(m4gnspecies_trans,m4gngspecies_trans,galt);
m4gngrxn_tgn = Rxn(m4gnspecies_tgn,m4gngspecies_tgn,galt);
%m4gngrxn.setrxnkinetics;
allrxns.add(m4gngrxn_cis);
allrxns.add(m4gngrxn_media);
allrxns.add(m4gngrxn_trans);
allrxns.add(m4gngrxn_tgn);

m3gngrxn_cis = Rxn(m3gnspecies_cis,m3gngspecies_cis,galt);
m3gngrxn_media = Rxn(m3gnspecies_media,m3gngspecies_media,galt);
m3gngrxn_trans = Rxn(m3gnspecies_trans,m3gngspecies_trans,galt);
m3gngrxn_tgn = Rxn(m3gnspecies_tgn,m3gngspecies_tgn,galt);
%m3gngrxn.setrxnkinetics;
allrxns.add(m3gngrxn_cis);
allrxns.add(m3gngrxn_media);
allrxns.add(m3gngrxn_trans);
allrxns.add(m3gngrxn_tgn);

m3gn2grxn_cis = Rxn(m3gn2species_cis,m3gn2gspecies_cis,galt);
m3gn2grxn_media = Rxn(m3gn2species_media,m3gn2gspecies_media,galt);
m3gn2grxn_trans = Rxn(m3gn2species_trans,m3gn2gspecies_trans,galt);
m3gn2grxn_tgn = Rxn(m3gn2species_tgn,m3gn2gspecies_tgn,galt);
%m3gn2grxn.setrxnkinetics;
allrxns.add(m3gn2grxn_cis);
allrxns.add(m3gn2grxn_media);
allrxns.add(m3gn2grxn_trans);
allrxns.add(m3gn2grxn_tgn);

m3gn3bgrxn_cis = Rxn(m3gn3bspecies_cis,m3gn3bgspecies_cis,galt);
m3gn3bgrxn_media = Rxn(m3gn3bspecies_media,m3gn3bgspecies_media,galt);
m3gn3bgrxn_trans = Rxn(m3gn3bspecies_trans,m3gn3bgspecies_trans,galt);
m3gn3bgrxn_tgn = Rxn(m3gn3bspecies_tgn,m3gn3bgspecies_tgn,galt);
%m3gn3bgrxn.setrxnkinetics;
allrxns.add(m3gn3bgrxn_cis);
allrxns.add(m3gn3bgrxn_media);
allrxns.add(m3gn3bgrxn_trans);
allrxns.add(m3gn3bgrxn_tgn);

m3gn3grxn_cis = Rxn(m3gn3species_cis,m3gn3gspecies_cis,galt);
m3gn3grxn_media = Rxn(m3gn3species_media,m3gn3gspecies_media,galt);
m3gn3grxn_trans = Rxn(m3gn3species_trans,m3gn3gspecies_trans,galt);
m3gn3grxn_tgn = Rxn(m3gn3species_tgn,m3gn3gspecies_tgn,galt);
%m3gn3grxn.setrxnkinetics;
allrxns.add(m3gn3grxn_cis);
allrxns.add(m3gn3grxn_media);
allrxns.add(m3gn3grxn_trans);
allrxns.add(m3gn3grxn_tgn);

m3gn4grxn_cis = Rxn(m3gn4species_cis,m3gn4gspecies_cis,galt);
m3gn4grxn_media = Rxn(m3gn4species_media,m3gn4gspecies_media,galt);
m3gn4grxn_trans = Rxn(m3gn4species_trans,m3gn4gspecies_trans,galt);
m3gn4grxn_tgn = Rxn(m3gn4species_tgn,m3gn4gspecies_tgn,galt);
%m3gn4grxn.setrxnkinetics;
allrxns.add(m3gn4grxn_cis);
allrxns.add(m3gn4grxn_media);
allrxns.add(m3gn4grxn_trans);
allrxns.add(m3gn4grxn_tgn);

m5gngnbrxn_cis = Rxn(m5gnspecies_cis,m5gngnbspecies_cis,mgat3);
m5gngnbrxn_media = Rxn(m5gnspecies_media,m5gngnbspecies_media,mgat3);
m5gngnbrxn_trans = Rxn(m5gnspecies_trans,m5gngnbspecies_trans,mgat3);
m5gngnbrxn_tgn = Rxn(m5gnspecies_tgn,m5gngnbspecies_tgn,mgat3);
%m5gngnbrxn.setrxnkinetics;
allrxns.add(m5gngnbrxn_cis);
allrxns.add(m5gngnbrxn_media);
allrxns.add(m5gngnbrxn_trans);
allrxns.add(m5gngnbrxn_tgn);

m4gngnbrxn_cis = Rxn(m4gnspecies_cis,m4gngnbspecies_cis,mgat3);
m4gngnbrxn_media = Rxn(m4gnspecies_media,m4gngnbspecies_media,mgat3);
m4gngnbrxn_trans = Rxn(m4gnspecies_trans,m4gngnbspecies_trans,mgat3);
m4gngnbrxn_tgn = Rxn(m4gnspecies_tgn,m4gngnbspecies_tgn,mgat3);
%m4gngnbrxn.setrxnkinetics;
allrxns.add(m4gngnbrxn_cis);
allrxns.add(m4gngnbrxn_media);
allrxns.add(m4gngnbrxn_trans);
allrxns.add(m4gngnbrxn_tgn);

m3gngnbrxn_cis = Rxn(m3gnspecies_cis,m3gngnbspecies_cis,mgat3);
m3gngnbrxn_media = Rxn(m3gnspecies_media,m3gngnbspecies_media,mgat3);
m3gngnbrxn_trans = Rxn(m3gnspecies_trans,m3gngnbspecies_trans,mgat3);
m3gngnbrxn_tgn = Rxn(m3gnspecies_tgn,m3gngnbspecies_tgn,mgat3);
%m3gngnbrxn.setrxnkinetics;
allrxns.add(m3gngnbrxn_cis);
allrxns.add(m3gngnbrxn_media);
allrxns.add(m3gngnbrxn_trans);
allrxns.add(m3gngnbrxn_tgn);

m3gn2gnbrxn_cis = Rxn(m3gn2species_cis,m3gn2gnbspecies_cis,mgat3);
m3gn2gnbrxn_media = Rxn(m3gn2species_media,m3gn2gnbspecies_media,mgat3);
m3gn2gnbrxn_trans = Rxn(m3gn2species_trans,m3gn2gnbspecies_trans,mgat3);
m3gn2gnbrxn_tgn = Rxn(m3gn2species_tgn,m3gn2gnbspecies_tgn,mgat3);
%m3gn2gnbrxn.setrxnkinetics;
allrxns.add(m3gn2gnbrxn_cis);
allrxns.add(m3gn2gnbrxn_media);
allrxns.add(m3gn2gnbrxn_trans);
allrxns.add(m3gn2gnbrxn_tgn);

m3gn3bgnbrxn_cis = Rxn(m3gn3bspecies_cis,m3gn3bgnbspecies_cis,mgat3);
m3gn3bgnbrxn_media = Rxn(m3gn3bspecies_media,m3gn3bgnbspecies_media,mgat3);
m3gn3bgnbrxn_trans = Rxn(m3gn3bspecies_trans,m3gn3bgnbspecies_trans,mgat3);
m3gn3bgnbrxn_tgn = Rxn(m3gn3bspecies_tgn,m3gn3bgnbspecies_tgn,mgat3);
%m3gn3bgnbrxn.setrxnkinetics;
allrxns.add(m3gn3bgnbrxn_cis);
allrxns.add(m3gn3bgnbrxn_media);
allrxns.add(m3gn3bgnbrxn_trans);
allrxns.add(m3gn3bgnbrxn_tgn);

m3gn3gnbrxn_cis = Rxn(m3gn3species_cis,m3gn3gnbspecies_cis,mgat3);
m3gn3gnbrxn_media = Rxn(m3gn3species_media,m3gn3gnbspecies_media,mgat3);
m3gn3gnbrxn_trans = Rxn(m3gn3species_trans,m3gn3gnbspecies_trans,mgat3);
m3gn3gnbrxn_tgn = Rxn(m3gn3species_tgn,m3gn3gnbspecies_tgn,mgat3);
%m3gn3gnbrxn.setrxnkinetics;
allrxns.add(m3gn3gnbrxn_cis);
allrxns.add(m3gn3gnbrxn_media);
allrxns.add(m3gn3gnbrxn_trans);
allrxns.add(m3gn3gnbrxn_tgn);

m3gn4gnbrxn_cis = Rxn(m3gn4species_cis,m3gn4gnbspecies_cis,mgat3);
m3gn4gnbrxn_media = Rxn(m3gn4species_media,m3gn4gnbspecies_media,mgat3);
m3gn4gnbrxn_trans = Rxn(m3gn4species_trans,m3gn4gnbspecies_trans,mgat3);
m3gn4gnbrxn_tgn = Rxn(m3gn4species_tgn,m3gn4gnbspecies_tgn,mgat3);
%m3gn4gnbrxn.setrxnkinetics;
allrxns.add(m3gn4gnbrxn_cis);
allrxns.add(m3gn4gnbrxn_media);
allrxns.add(m3gn4gnbrxn_trans);
allrxns.add(m3gn4gnbrxn_tgn);

m5gngnbgrxn_cis = Rxn(m5gngnbspecies_cis,m5gngnbgspecies_cis,galt);
m5gngnbgrxn_media = Rxn(m5gngnbspecies_media,m5gngnbgspecies_media,galt);
m5gngnbgrxn_trans = Rxn(m5gngnbspecies_trans,m5gngnbgspecies_trans,galt);
m5gngnbgrxn_tgn = Rxn(m5gngnbspecies_tgn,m5gngnbgspecies_tgn,galt);
%m5gngnbgrxn.setrxnkinetics;
allrxns.add(m5gngnbgrxn_cis);
allrxns.add(m5gngnbgrxn_media);
allrxns.add(m5gngnbgrxn_trans);
allrxns.add(m5gngnbgrxn_tgn);

m4gngnbgrxn_cis = Rxn(m4gngnbspecies_cis,m4gngnbgspecies_cis,galt);
m4gngnbgrxn_media = Rxn(m4gngnbspecies_media,m4gngnbgspecies_media,galt);
m4gngnbgrxn_trans = Rxn(m4gngnbspecies_trans,m4gngnbgspecies_trans,galt);
m4gngnbgrxn_tgn = Rxn(m4gngnbspecies_tgn,m4gngnbgspecies_tgn,galt);
%m4gngnbgrxn.setrxnkinetics;
allrxns.add(m4gngnbgrxn_cis);
allrxns.add(m4gngnbgrxn_media);
allrxns.add(m4gngnbgrxn_trans);
allrxns.add(m4gngnbgrxn_tgn);

m3gngnbgrxn_cis = Rxn(m3gngnbspecies_cis,m3gngnbgspecies_cis,galt);
m3gngnbgrxn_media = Rxn(m3gngnbspecies_media,m3gngnbgspecies_media,galt);
m3gngnbgrxn_trans = Rxn(m3gngnbspecies_trans,m3gngnbgspecies_trans,galt);
m3gngnbgrxn_tgn = Rxn(m3gngnbspecies_tgn,m3gngnbgspecies_tgn,galt);
%m3gngnbgrxn.setrxnkinetics;
allrxns.add(m3gngnbgrxn_cis);
allrxns.add(m3gngnbgrxn_media);
allrxns.add(m3gngnbgrxn_trans);
allrxns.add(m3gngnbgrxn_tgn);

m3gn2gnbgrxn_cis = Rxn(m3gn2gnbspecies_cis,m3gn2gnbgspecies_cis,galt);
m3gn2gnbgrxn_media = Rxn(m3gn2gnbspecies_media,m3gn2gnbgspecies_media,galt);
m3gn2gnbgrxn_trans = Rxn(m3gn2gnbspecies_trans,m3gn2gnbgspecies_trans,galt);
m3gn2gnbgrxn_tgn = Rxn(m3gn2gnbspecies_tgn,m3gn2gnbgspecies_tgn,galt);
%m3gn2gnbgrxn.setrxnkinetics;
allrxns.add(m3gn2gnbgrxn_cis);
allrxns.add(m3gn2gnbgrxn_media);
allrxns.add(m3gn2gnbgrxn_trans);
allrxns.add(m3gn2gnbgrxn_tgn);

m3gn3bgnbgrxn_cis = Rxn(m3gn3bgnbspecies_cis,m3gn3bgnbgspecies_cis,galt);
m3gn3bgnbgrxn_media = Rxn(m3gn3bgnbspecies_media,m3gn3bgnbgspecies_media,galt);
m3gn3bgnbgrxn_trans = Rxn(m3gn3bgnbspecies_trans,m3gn3bgnbgspecies_trans,galt);
m3gn3bgnbgrxn_tgn = Rxn(m3gn3bgnbspecies_tgn,m3gn3bgnbgspecies_tgn,galt);
%m3gn3bgnbgrxn.setrxnkinetics;
allrxns.add(m3gn3bgnbgrxn_cis);
allrxns.add(m3gn3bgnbgrxn_media);
allrxns.add(m3gn3bgnbgrxn_trans);
allrxns.add(m3gn3bgnbgrxn_tgn);

m3gn3gnbgrxn_cis = Rxn(m3gn3gnbspecies_cis,m3gn3gnbgspecies_cis,galt);
m3gn3gnbgrxn_media = Rxn(m3gn3gnbspecies_media,m3gn3gnbgspecies_media,galt);
m3gn3gnbgrxn_trans = Rxn(m3gn3gnbspecies_trans,m3gn3gnbgspecies_trans,galt);
m3gn3gnbgrxn_tgn = Rxn(m3gn3gnbspecies_tgn,m3gn3gnbgspecies_tgn,galt);
%m3gn3gnbgrxn.setrxnkinetics;
allrxns.add(m3gn3gnbgrxn_cis);
allrxns.add(m3gn3gnbgrxn_media);
allrxns.add(m3gn3gnbgrxn_trans);
allrxns.add(m3gn3gnbgrxn_tgn);

m3gn4gnbgrxn_cis = Rxn(m3gn4gnbspecies_cis,m3gn4gnbgspecies_cis,galt);
m3gn4gnbgrxn_media = Rxn(m3gn4gnbspecies_media,m3gn4gnbgspecies_media,galt);
m3gn4gnbgrxn_trans = Rxn(m3gn4gnbspecies_trans,m3gn4gnbgspecies_trans,galt);
m3gn4gnbgrxn_tgn = Rxn(m3gn4gnbspecies_tgn,m3gn4gnbgspecies_tgn,galt);
%m3gn4gnbgrxn.setrxnkinetics;
allrxns.add(m3gn4gnbgrxn_cis);
allrxns.add(m3gn4gnbgrxn_media);
allrxns.add(m3gn4gnbgrxn_trans);
allrxns.add(m3gn4gnbgrxn_tgn);

% set up enzyme pool
allenzs = CellArrayList;
allenzs.add(mgat1);
allenzs.add(mgat2);
allenzs.add(mgat3);
allenzs.add(mgat4);
allenzs.add(mgat5);
allenzs.add(mani);
allenzs.add(manii);
allenzs.add(galt);

testpathway = Pathway;
testpathway.addGlycans(allglycans);
testpathway.addRxns(allrxns)
testpathway.addEnzs(allenzs);

testpathway.setRxnsID;
testpathway.setEnzID;
testpathway.setEnzRxnsKinetics;

testpathway.setComptReactorType(Reactor.CSTR);  
% set qflow and k in each compts
qflow  = 10.8;
k      = 10.8;
testpathway.setComptsQflow(qflow);
testpathway.setComptskt(k);

%first compt
speciesInletConc = CellArrayList;
speciesInletConc.add(struct('species',m9species_cis ,'inletconc',1));
speciesInletConc.add(struct('species',m8species_cis ,'inletconc',0));
speciesInletConc.add(struct('species',m7species_cis ,'inletconc',0));
speciesInletConc.add(struct('species',m6species_cis ,'inletconc',0));
speciesInletConc.add(struct('species',m5species_cis,'inletconc',0));
speciesInletConc.add(struct('species',m5gnspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m4gnspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gnspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn2species_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn3species_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn3bspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn4species_cis,'inletconc',0));
speciesInletConc.add(struct('species',m5gngspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m4gngspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gngspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn2gspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn3gspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn3bgspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn4gspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m5gngnbspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m4gngnbspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gngnbspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn2gnbspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn3gnbspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn3bgnbspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn4gnbspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m5gngnbgspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m4gngnbgspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gngnbgspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn2gnbgspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn3gnbgspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn3bgnbgspecies_cis,'inletconc',0));
speciesInletConc.add(struct('species',m3gn4gnbgspecies_cis,'inletconc',0));

cisCompt.reactortype.setSpeciesInletDB(speciesInletConc)

testpathway.setTransporRxnKinetics;

modelname ='testSBML';
testmodel = GlycanNetModel(theCompt,testpathway,modelname);
testsbmlstruct= testmodel.toSBMLStruct;
glycanNetViewer(testmodel);
OutputSBML(testsbmlstruct);
end

function listofglycanstruct = setGlyStruct()
listofstructs = CellArrayList;
listofstructs.add((glycanMLread('M9.glycoct_xml')));
listofstructs.add((glycanMLread('M8.glycoct_xml')));
listofstructs.add((glycanMLread('M7.glycoct_xml')));
listofstructs.add((glycanMLread('M6.glycoct_xml')));
listofstructs.add((glycanMLread('M5.glycoct_xml')));
listofstructs.add((glycanMLread('M5gn.glycoct_xml')));
listofstructs.add((glycanMLread('M4gn.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn2.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn3.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn3b.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn4.glycoct_xml')));
listofstructs.add((glycanMLread('M5gng.glycoct_xml')));
listofstructs.add((glycanMLread('M4gng.glycoct_xml')));
listofstructs.add((glycanMLread('M3gng.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn2g.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn3g.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn3bg.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn4g.glycoct_xml')));
listofstructs.add((glycanMLread('M5gngnb.glycoct_xml')));
listofstructs.add((glycanMLread('M4gngnb.glycoct_xml')));
listofstructs.add((glycanMLread('M3gngnb.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn2gnb.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn3gnb.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn3bgnb.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn4gnb.glycoct_xml')));
listofstructs.add((glycanMLread('M5gngnbg.glycoct_xml')));
listofstructs.add((glycanMLread('M4gngnbg.glycoct_xml')));
listofstructs.add((glycanMLread('M3gngnbg.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn2gnbg.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn3gnbg.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn3bgnbg.glycoct_xml')));
listofstructs.add((glycanMLread('M3gn4nbg.glycoct_xml')));
end
