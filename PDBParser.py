class PDB_residueException(Exception):
        """ Generic exception from PDB_residue errors """        
        pass


class PDB_atom:
    """ Main class which holds an induvidual atom from a PBD file.
        Carries out the parsing of an atom line from a PDB file into a valid
        PDB_atom object.

        This is the ONLY place where parsing of the PDB atom lines should 
        occur. All further logic based on chain name, residue ID, residue
        name
    """


    class PDB_atomException(Exception):
        """ Generic exception from PDB_atom errors """        
        pass

    def __init__(self, line):
        self.parse(line)

    def parse(self, line):        
        """ This is the initialization parser. Converts a line
            from a PDB file and initializes the PDB_atom's variables
            to appropriate values. Uses implicit typecasting as a 
            failcheck for parsing the file correctly.

            Well formatted PDB files are 80 characters across and are
            well defined/easy to parse.
            
            Sadly not all PDB files are like this, so we use a heuristic
            to parse PDB files which are not 80 characters across.

            # INPUT
            line      :     String

            # OUTPUT
            -         :     None
        """



        # correctly formatted PDB files
        # TODO - assuming 80 chars means well formatted is
        #        perhaps risky. Need a more robust way to asses
        #        formatting validity
        if len(line) == 80:
            self.record_name    = line[0:6].strip()
            self.atom_id    = int(line[6:11].strip())
            self.atom_name      = line[12:16].strip()
            self.alt_location   = line[16]
            self.res_name       = line[17:20].strip()
            self.chain          = line[21]
            self.res_iod         = line[22:26].strip()
            self.res_ins_code   = line[26]
            self.coord_X        = float(line[30-38].strip())
            self.coord_Y        = float(line[38-46].strip())
            self.coord_Z        = float(line[46-54].strip())
            self.occupancy      = float(line[54-60].strip())
            self.beta           = float(line[60-66].strip())
            self.seg_ID         = line[72:76].strip()
            self.element        = line[76:78].strip()
            self.charge         = float(line[78:80].strip())
            self.chain_local_id = -1
            self.formatted_ok   = True

        # Heuristic section - split by space and then use
        # errors in casting as flags for things being issues
        # Note this may need to be expanded as malformed edge-cases
        # are identified...
        else:
            rawsplitline = filter(None, line.split(" "))


            splitline = []
            for i in rawsplitline:
                if i == "\n" or i == "\t":
                    pass
                else:
                    splitline.append(i)
                    
            num_cols = len(splitline)
            
            try:
                if num_cols == 10:
                    self.record_name    = splitline[0]   
                    self.atom_id        = int(splitline[1])
                    self.atom_name      = splitline[2]   
                    self.alt_location   = ""
                    self.res_name       = splitline[3]   
                    self.chain          = ""
                    self.res_id         = int(splitline[4])
                    self.res_ins_code   = ""
                    self.coord_X        = float(splitline[5])  
                    self.coord_Y        = float(splitline[6])  
                    self.coord_Z        = float(splitline[7])
                    self.occupancy       = float(splitline[8])
                    self.beta           = float(splitline[9])
                    self.seg_ID         = " "
                    self.element        = " "                
                    self.charge         = " "
                    self.chain_local_id = -1
                    self.formatted_ok   = False

                elif num_cols == 11:
                    self.record_name    = splitline[0]   
                    self.atom_id        = int(splitline[1])
                    self.atom_name      = splitline[2]   
                    self.alt_location   = " "
                    self.res_name       = splitline[3]   
                    self.chain          = splitline[4]
                    self.res_id         = int(splitline[5])
                    self.res_ins_code   = " "
                    self.coord_X        = float(splitline[6])  
                    self.coord_Y        = float(splitline[7])  
                    self.coord_Z        = float(splitline[8])  
                    self.occupancy       = float(splitline[9]) 
                    self.beta           = float(splitline[10])
                    self.seg_ID         = " "
                    self.element        = " "                
                    self.charge         = " "
                    self.chain_local_id = -1
                    self.formatted_ok   = False
                elif num_cols == 11:
                    self.record_name    = splitline[0]   
                    self.atom_id        = int(splitline[1])
                    self.atom_name      = splitline[2]   
                    self.alt_location   = " "
                    self.res_name       = splitline[3]   
                    self.chain          = splitline[4]
                    self.res_id         = int(splitline[5])
                    self.res_ins_code   = " "
                    self.coord_X        = float(splitline[6])  
                    self.coord_Y        = float(splitline[7])  
                    self.coord_Z        = float(splitline[8])  
                    self.occupancy       = float(splitline[9]) 
                    self.beta           = float(splitline[10])
                    self.seg_ID         = " "
                    self.element        = splitline[11]      
                    self.charge         = " "
                    self.chain_local_id = -1
                    self.formatted_ok   = False
                else:
                    raise PDB_atomException("Did not match number of columns")
            except ValueError,e:
                print "Error with columns (using " + str(num_cols) + ") columns"
                print "Tried to cast string to int/float"
                raise e
                                            
                                            
                    
    ## Overwritten default behaviours
    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "<PDB_Atom " + self.res_name + str(self.res_id) + " -> " + self.atom_name + "[atom " + str(self.atom_id) + "]>"
                    

class PDB_residue:
    """ Class for holding a single residue. """

    def __init__(self, residue_atoms_list):
        self.atoms          = residue_atoms_list
        self.res_name       = residue_atoms_list[0].res_name
        self.res_id         = residue_atoms_list[0].res_id
        self.chain          = residue_atoms_list[0].chain
        self.chain_local_id = residue_atoms_list[0].chain_local_id

    def get_alpha_carbon(self):
        for i in self.atoms:
            if i.atom_name == "CA":
                return i
            
    def set_residue_order(self, atomname_list):

        newOrder = []

        # collect the ordered atom_ids to reassign - cannot
        # assume they're incrmented by 1 each time (probably 
        # are but...)
        new_atom_ids = []
        for atom in self.atoms:
            new_atom_ids.append(atom.atom_id)

        id_index = 0

        for atomname in atomname_list:
            for atom in self.atoms:
                if atomname == atom.atom_name:
                    atom.atom_id = new_atom_ids[id_index]
                    
                    newOrder.append(atom)
                    id_index=id_index+1
        
        if not len(self.atoms) == len(newOrder):
            msg = "When reseting atom order in  residue must fully define the new order\n" + \
                  "Old = " + str(self.atoms) + "\n" + \
                  "New = " + str(newOrder) + "\n"
                  
            raise PDB_residueException(msg)

        self.atoms = newOrder

    def rename_residue(self, newName):
        for atom in self.atoms:
            atom.res_name = newName
        self.res_name = newName

    def rename_atom(self, oldName, newName):
        
        oldName = str(oldName)
        newName = str(newName)
        
        for atom in self.atoms:
            if atom.atom_name == str(oldName):
                atom.atom_name = str(newName)
                return
        msg = "ERROR: Unable to find " + oldName + " in residue " + str(self.res_id) + "(" + self.res_name  + ") in chain " + self.chain + " (to replace with " + newName + ")"
        raise PDB_residueException(msg)

    ## Overwritten default behaviours

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "<PDB_Residue " + str(self.atoms[0].res_id) + " " + self.atoms[0].res_name + ">"

    def __getitem__(self,indx):
        return self.atoms[indx]

    def __len__(self):
        return len(self.atoms)
    
    


class PDB_residue_organizer:

    def __init__(self):        
        pass

    def construct_chains(self, atomlines):

        chains = {}          
        parsed_atoms_list = self.__parse_atoms(atomlines)
        
        chainID_list = self.__get_chainlist(parsed_atoms_list)        

        for chainID in chainID_list:
            chains[chainID] = PDB_chain(self.__get_atoms_from_chain(parsed_atoms_list, chainID))  
            
        return chains


    def __get_chainlist(self, parsed_atoms_list):
        
        chainid_list = []

        for atom in parsed_atoms_list:
            if atom.chain in chainid_list:
                pass
            else:
                chainid_list.append(atom.chain)

        return chainid_list
            
        
    def __parse_atoms(self, atomlines):

        parsed_atoms_list = []
        
        for atom_line in atomlines:
            atom = PDB_atom(atom_line)
            parsed_atoms_list.append(atom)

        return parsed_atoms_list

    def __get_atoms_from_chain(self, parsed_atoms_list, chainid):
        
        chain_atoms =[]

        for atom in parsed_atoms_list:
            
            if atom.chain == chainid:
                chain_atoms.append(atom)

        return chain_atoms



class PDB_chain:

    def __init__(self, atomlist):                
        """ Initialization function which takes a list of atoms all from one chain
            and constructs a chain object. Importantly, the chain_local_id atom and residue
            level attribute is only initilized in the "context" of a chain - i.e. by this 
            initialization function, so this function actually completes atom and residue
            initialization.

            A crucial assumption is that all the atoms in atomlist come from the same chain.
            In fact - this list is how a chain object is defined, so if the list contains
            atoms/residues with different .chain values, things are going to go very wrong...
         """

        self.residues = []

        current_res = atomlist[0].res_id

        temp_res =[]

        for i in atomlist:
            if i.res_id == current_res:
                temp_res.append(i)
            else:
                self.residues.append(PDB_residue(temp_res))
                temp_res = []
                current_res = i.res_id
                temp_res.append(i)
                
                

        self.residues.append(PDB_residue(temp_res)) # get final residue

        if len(self.residues) > 0:            
            self.chain_name = self.residues[0].chain
        else:
            self.chain_name = None
        
        chain_local_id = 1        

        ## CHAIN LOCAL ID variables in atom and residue
        ## objects are set here!
        ##
        for res in self.residues:
            for atom in res.atoms:
                atom.chain_local_id = chain_local_id
            res.chain_local_id = chain_local_id
            chain_local_id=chain_local_id+1

            

    def get_residue(self, resid, chainLocal=False):
        """
            Function which returns the residue from a chain.

            If chainLocal = False (default) the global resID is used
            to search for a residue. If, on the other hand, chainLocal
            is set to true, then the resid local to that chain is used.
            
            chainLocal ID values ALWAYS start from 1 for the first residue
            in a chain and increment by one through the chain. They provide
            an easy internal way to get (for example) the first and last
            residue in the chain explicitly, or comparative residues when
            dealing with many identical species.
        """
        
        try:
            resid = int(resid)
        except ValueError, e:
            print "ERROR: get_resid() requires a numeric value"
            raise e
            

        if type(resid) == int:   
            
            if chainLocal:
                for res in self.residues:
                    if res.chain_local_id == resid:                        
                        return res
            else:                
                for res in self.residues:
                    if res.res_id == resid:
                        return res


    def get_chain_length(self):
        return len(self.residues)

    def __str__(self):
        if self.chain_name:
            return "<PDB_chain " + str(self.chain_name) + " - [" + str(len(self.residues)) + " residues]>"
        else:
            return "<PDB_chain of length " + str(len(self.residues)) + ">"
            
    def __repr__(self):
        return self.__str__()

    def __getitem__(self,indx):
        return self.residues[indx]

    def __len__(self):
        return len(self.residues)

        


class PDB_file:       
    class PDB_fileException(Exception):
        """ Generic exception from PDB_file errors """        
        pass

    def __init__(self, filename):
        content = self.__read_file(filename)
        self.chains = self.__parse_residues(content)
        self.header = self.__get_header(content)
        self.footer = self.__get_footer(content)

        self.residues = self.__get_residues_from_chains()
        self.atoms = self.__get_atoms_from_chain()
        
    def write_file(self,filename):
        f = open(filename,'w')
        for line in self.header:
            f.write(line)
        
        self.__write_chains(f)

        for line in self.footer:
            f.write(line)            

    def rename_residue(self, chainID, resID, newName):
        chain = self.chains[chainID]
        residue = chain.get_residue(resID, chainLocal=True)
        residue.rename_residue(newName)

    def define_residue_order(self, chainID, resID, atomOrder):
        """ Function where you can specify the atom order for a residue.
        Useful because different forcefields have annoying defaults
        for how residues are parsed which *occasionally* differ, so
        this provides an easy way to define the order and re-write.
        """

        chain = self.chains[chainID]
        residue = chain.get_residue(resID, chainLocal=True)
        residue.set_residue_order(atomOrder)

    def rename_atom(self, chainID, resID, oldName, newName):
        

        chain = self.chains[chainID]
        residue = chain.get_residue(resID, chainLocal=True)
        
        residue.rename_atom(oldName, newName)
        

    def __read_file(self, filename):        
        """ Reads a file into a content list line by line                                                                                                                                                                                    
        and returns that list
        """
        with open(filename) as f:
            content = f.readlines()
            
        return content

    def __get_header(self, content):
        header = [] 

        for line in content:
            if filter(None, line.split(" "))[0].upper() == "ATOM": 
                break
            header.append(line)
            
        return header

    def __get_footer(self, content):
        """ 
	"""
	footer = [] 

	passedAtoms = False
	for line in content:
            if filter(None, line.split(" "))[0].upper() == "ATOM" or filter(None, line.split(" "))[0].upper() == "TER": 
                passedAtoms = True
	    else:
                if passedAtoms:
                    if line == "TER\n": # TER lines are actually dealt with in chains to signify the end of a chain"
                        pass
		    else:
                        footer.append(line)        
        return footer
        
        
    def __parse_residues(self, content):
        """ Main parsing function build object based chains """
        
	atomlines =[]

        # create a list of atom lines
        for line in content:
            if filter(None, line.split(" "))[0].upper() == "ATOM": 
                atomlines.append(line)
               
        organizer = PDB_residue_organizer()

        return(organizer.construct_chains(atomlines))

    def __get_residues_from_chains(self):
        
        residues = []


        for i in self.chains:
            chain = self.chains[i]
            
            for res in chain.residues:
                residues.append(res)
                
        return residues
            
    def __get_atoms_from_chain(self):

        atoms =[]

        for res in self.residues:      
            for atom in res.atoms:
                atoms.append(atom)

        return atoms


            
        

    def __write_chains(self, handle):
        for chainname in self.chains:
            chain = self.chains[chainname]

            for residue in chain:
                for atom in residue:
                    self.__atom_write(handle,atom)

            handle.write("TER\n")

    def __atom_write(self, handle, atom):
        record_name = self.__string_padder(atom.record_name, 6, "L")
        atom_id     = self.__string_padder(atom.atom_id, 5, "R")
        space1       = " "
        

        if len(str(atom.atom_name)) < 4:
            space2 = " "
            atom_name   = self.__string_padder(atom.atom_name, 3, "L")
            
        elif len(str(atom.atom_name)) == 4:
            space2 = ""
            atom_name = atom.atom_name
        else:
            raise PDB_fileException("ERROR writing PDB formatting atom name")

        alt_location = self.__string_padder(atom.alt_location, 1, "C")
        res_name = self.__string_padder(atom.res_name,3, "C")
        space3 = " "
        chain = self.__string_padder(atom.chain, 1, "R")
        
        res_seq = self.__string_padder(atom.res_id, 4, "R")
        res_ins_code = self.__string_padder(atom.res_ins_code, 1, "R")
        space4 = "   "
        
        xcoord = self.__string_padder("%.3f" % atom.coord_X, 8, "R")
        ycoord = self.__string_padder("%.3f" % atom.coord_Y, 8, "R")
        zcoord = self.__string_padder("%.3f" % atom.coord_Z, 8, "R")
        occupancy = self.__string_padder(atom.occupancy, 6, "R")
        beta = self.__string_padder(atom.beta, 6, "R")
        space5 = "       "
        segID = self.__string_padder(atom.seg_ID, 4, "L")
        element = self.__string_padder(atom.element, 1, "R")
        charge = self.__string_padder(atom.charge, 1, "R")

        # now build the actual line from the various components
        outline = record_name + atom_id + space1 + space2 + atom_name + \
                  alt_location + res_name + space3 + chain + res_seq + \
                  res_ins_code + space4 + xcoord + ycoord + zcoord + \
                  occupancy + beta + space5 + segID + element + charge + "\n"

        handle.write(outline)

    def __string_padder(self, item,length,justification):

        if len(str(item)) > length:
            errormsg = "Error when writing PDB file - field is wider than allowed based on PDB formatting: " +\
                        str(item) + " - must be equal to or less than " + str(length) + " characters"
            raise PDB_fileException(errormsg)
                                    

        padding = length - len(str(item))
        
        if justification == "L":
            return str(item)+padding*" "
        if justification == "R":
            return padding*" " + str(item)
        if justification == "C":
            if padding % 2 == 0:
                RP = padding/2
                LP = padding/2
            else:
                RP = padding/2
                LP = (padding/2)+1
            return RP*" " + str(item) + LP*" "
                    
    def __len__(self):
        return len(self.chains)


    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "<PDB_file of " + str(len(self.chains)) + " chains, " + str(len(self.residues)) + " residues and " +  str(len(self.atoms)) + " atoms>"

    def __getitem__(self, key):
        if key in self.chains:
            return self.chains[key]
        else:
            print "No chain of name " + str(key) + " in file"
            return None

    
