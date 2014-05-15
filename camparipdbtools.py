import PDBParser

class CAMPARI_pdbException(Exception):
    pass

class CAMPARI_pdb:


    def __init__(self, filename):

        self.pdbfile = PDBParser.PDB_file(filename)

        
    def write_file(self, outputfilename):
        self.pdbfile.write_file(outputfilename)
        

    def convert_from_CAMPARI_to_GMX(self, capChange=True):

        # reorder ACE and NAC if needed


        chains = self.pdbfile.chains

        for chainID in chains:

            if capChange:

                if chains[chainID][0].res_name == "ACE":                            
                    self.pdbfile.define_residue_order(chainID, 1, ["CH3", "1H","2H","3H","C", "O"])

                    # rename ACE atoms...
                    self.pdbfile.rename_atom(chainID, 1, "1H", "1HH3")
                    self.pdbfile.rename_atom(chainID, 1, "2H", "2HH3")
                    self.pdbfile.rename_atom(chainID, 1, "3H", "3HH3")


                if chains[chainID][len(chains[chainID])-1].res_name == "NME":
                
                    finalRes = len(chains[chainID])
                
                    self.pdbfile.define_residue_order(chainID, finalRes, ["N", "HN", "CH3", "1H", "2H", "3H"])
                    self.pdbfile.rename_residue(chainID, finalRes, "NAC")
                
                    self.pdbfile.rename_atom(chainID, finalRes, "HN", "H")
                    self.pdbfile.rename_atom(chainID, finalRes, "1H", "1HH3")
                    self.pdbfile.rename_atom(chainID, finalRes, "2H", "2HH3")
                    self.pdbfile.rename_atom(chainID, finalRes, "3H", "3HH3")


            # no need to make the HIE/HID/HIP->HIS correction as GROMACS can typically deal with one of these


    def convert_from_GMX_to_CAMPARI(self, capChange=False):

        # reorder ACE and NAC if needed
        chains = self.pdbfile.chains
        
        for chainID in chains:
            

            ## So sometimes it seems like re-setting the CAP atom orders is important?
            ## In any case, if this is desired you can set capChange to true, else just leave it
            ## false...
            ##
            if capChange:
                if chains[chainID][0].res_name == "ACE":            
                    self.pdbfile.define_residue_order(chainID, 1, ["CH3", "C", "O", "1HH3", "2HH3", "3HH3"])

                    # rename ACE atoms...
                    self.pdbfile.rename_atom(chainID, 1, "1HH3", "1H")
                    self.pdbfile.rename_atom(chainID, 1, "2HH3", "2H")
                    self.pdbfile.rename_atom(chainID, 1, "3HH3", "3H")


                if chains[chainID][len(chains[chainID])-1].res_name == "NAC":

                    finalRes = len(chains[chainID])

                    self.pdbfile.define_residue_order(chainID, finalRes, ["N", "CH3", "H", "1HH3", "2HH3", "3HH3"])

                    # rename NAC to NME
                    self.pdbfile.rename_residue(chainID, finalRes, "NME")

               
                    self.pdbfile.rename_atom(chainID, finalRes, "H", "HN")
                    self.pdbfile.rename_atom(chainID, finalRes, "1HH3", "1H")
                    self.pdbfile.rename_atom(chainID, finalRes, "2HH3", "2H")
                    self.pdbfile.rename_atom(chainID, finalRes, "3HH3", "3H")

                
            chainlocal_residue_count=1
            for res in chains[chainID]:


                

                ##----------------------------------------------------------------------
                ## HISTADINE RENAMING
                ##
                # if you find a histadine CAMPARI expects HIE or HID
                # figure out where the hydrogen is and deal accordingly
                
                if res.res_name == "HIS":

                    # set initialization counts for histadine nitrogen 
                    # protons
                    HIE=0
                    HID=0

                    for atom in res:
                        if atom.atom_name == "HE1" or atom.atom_name == "HE2" :
                            HIE=HIE+1
                        if atom.atom_name == "HD1" or atom.atom_name == "HD2" :
                            HID=HID+1

                    if HIE == 2 and HID == 2:
                        self.pdbfile.rename_residue(chainID, chainlocal_residue_count, "HIP")
                    elif HIE == 2:
                        self.pdbfile.rename_residue(chainID, chainlocal_residue_count, "HIE")
                    elif HID == 2:
                        self.pdbfile.rename_residue(chainID, chainlocal_residue_count, "HID")
                    else:
                        raise CAMPARI_pdbException("ERROR: Histadine residue " + str(res.res_id) + " has an odd protonation state...")        
                ##----------------------------------------------------------------------


                # increment residue counter
                chainlocal_residue_count=chainlocal_residue_count+1
                            

