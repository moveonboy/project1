function charresidue = linkwithresidue(charresidue,residue)
charresidue = [charresidue '--'];
charresidue = [charresidue,residue.linkageParent.bonds.posParent];
charresidue = [charresidue,residue.anomer.symbol,residue.anomer.carbonPos];
charresidue = [charresidue,residue.stereoConfig.symbol];
charresidue = [charresidue,'-',residue.residueType.name];
charresidue = [charresidue,',p'];
end