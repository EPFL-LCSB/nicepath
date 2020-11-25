import re

class Reaction:
    def __init__(self):
        self.name = ''
        self.eligible = False
        self.old_entries = []
        self.entry = ''
        self.kegg = ''
        self.fileline = ''
        self.equation_name = ''
        self.rpairs = []
        self.cars = [] # list of pairs
        self.operators = ''
        self.thermo_feasible = {'forward': True, 'reverse': True}
        self.mechanism_feasible = {'forward': True, 'reverse': True}
        self.is_in_pathway = False
        self.equation = ''
        self.children = []  # duplicate reactions
        self.reactants = [] # [entry1, entry2, entry3, entry4]
        self.stoichiometry = [] # [-1, -1, 1, 2]
        self.reaction_hash = ''
        self.is_known = False



    def loadReactionFromLine(self, _line, _headerlist, _is_kegg_linking):
        linelist = _line.rstrip().split(';')
        self.fileline = _line
        if 'IS_ELIGIBLE' in _headerlist: self.eligible = bool(int(linelist[_headerlist.index('IS_ELIGIBLE')]))
        if 'KEGG' in _headerlist: self.kegg = linelist[_headerlist.index('KEGG')]
        if 'ENTRY' in _headerlist: self.entry = linelist[_headerlist.index('ENTRY')]
        if 'NAME' in _headerlist: self.name = linelist[_headerlist.index('NAME')]
        if 'OPERATORS' in _headerlist: self.operators = ','.join(list(filter(None,linelist[_headerlist.index('OPERATORS')].split('|'))))
        if 'EQUATION_NAME' in _headerlist: self.equation_name = linelist[_headerlist.index('EQUATION_NAME')]
        if 'EQUATION' in _headerlist: self.equation = linelist[_headerlist.index('EQUATION')]
        if not _is_kegg_linking:
            self.parseEquation()
            if self.kegg != '':
                self.is_known = True
        else:
            self.is_known = True
            self.kegg = self.entry

        if linelist[_headerlist.index('RPAIRS')] == "": # if no reactant pairs calculated, stop here
            print("Warning: Reaction number %s had not reactant pair"%self.entry)
            return

        pairslist = linelist[_headerlist.index('RPAIRS')].split('|')
        for rpair in pairslist:
            if _is_kegg_linking:  # correctly handle tautomers in KEGG coverage file
                tautomer_matches = re.findall(r'(C\d{5}.*)_(C\d{5}.*)', rpair)
                self.rpairs.append((tautomer_matches[0][0], tautomer_matches[0][1]))
            else:
                self.rpairs.append((rpair.split('_')[0], rpair.split('_')[1]))

        for x in linelist[_headerlist.index('CONSERVEDRATIO')].split('|'):
            if x != '':
                self.cars.append(float(x))

    def determineThermoFeasibility(self):
        '''determineThermoFeasibility relies on the name equation for the moment, should be changed in future'''
        left_side = self.equation_name.split('<==>')[0]
        right_side = self.equation_name.split('<==>')[1]

        # if molecular oxygen on the product side
        if 'Oxygen' in right_side:
            self.thermo_feasible['forward'] = False
        elif 'Oxygen' in left_side:
            self.thermo_feasible['reverse'] = False

        # if CO2 on the substrate side
        if 'CO2' in right_side:
            self.thermo_feasible['reverse'] = False
        elif 'CO2' in left_side:
            self.thermo_feasible['forward'] = False


    def determineMechanismFeasibility(self):
        '''determineMechanismFeasibility relies on the name equation for the moment, should be changed in future'''
        # do not allow demethylation via SAM/SAH
        left_side = self.equation_name.split('<==>')[0]
        right_side = self.equation_name.split('<==>')[1]

        # if SAH -> SAM
        if 'S-Adenosyl-L-methionine' in right_side and 'S-Adenosyl-L-homocysteine' in left_side:
            self.mechanism_feasible['forward'] = False
        elif 'S-Adenosyl-L-methionine' in left_side and 'S-Adenosyl-L-homocysteine' in right_side:
            self.mechanism_feasible['reverse'] = False



    def isThermoFeasible(self, _direction):
         return self.thermo_feasible[_direction]

    def isMechanismFeasible(self, _direction):
         return self.mechanism_feasible[_direction]

    def markAsUsed(self):
        self.is_in_pathway = True

    def parseEquation(self):
        assert (self.equation != '') or ('<==>' in self.equation), 'ERROR in splitting equation %s for reaction entry %s'%(self.equation, self.entry)
        left = self.equation.split('<==>')[0][7:].split(' + ')
        right = self.equation.split('<==>')[1].split(' + ')
        for x in left:
            stoich_match = re.search(r'\((\d+)\)\s(\d+)',x)
            if stoich_match:
                stoich = -int(stoich_match.group(1))
                cpd = stoich_match.group(2)
            else:
                match = re.search(r'\d+', x)
                cpd = match.group()
                stoich = -1
            self.stoichiometry.append(stoich)
            self.reactants.append(cpd)
        for x in right:
            stoich_match = re.search(r'\((\d+)\)\s(\d+)',x)
            if stoich_match:
                stoich = int(stoich_match.group(1))
                cpd = stoich_match.group(2)
            else:
                match = re.search(r'\d+', x)
                cpd = match.group()
                stoich = 1
            self.stoichiometry.append(stoich)
            self.reactants.append(cpd)


    def generateReactionHash(self, _rootmap):
        if self.reaction_hash == '':
            left = []
            right = []
            for i, entry in enumerate(self.reactants):
                if self.stoichiometry[i] < 0 :
                    if self.stoichiometry[i] < -1:
                        left.append('(%d)%s'%(-self.stoichiometry[i],_rootmap[entry]))
                    else:
                        left.append(_rootmap[entry])
                if self.stoichiometry[i] > 0 :
                    if self.stoichiometry[i] > 1:
                        right.append('(%d)%s'%(self.stoichiometry[i],_rootmap[entry]))
                    else:
                        right.append(_rootmap[entry])
            self.reaction_hash = '-'.join(sorted(['+'.join(sorted(left)),'+'.join(sorted(right))]))
        return self.reaction_hash

    def addChild(self, other):
        self.children.append(other)

    def getAllDependentCompounds(self):
        entry_list = set()
        entry_list.update(set(self.reactants))
        [entry_list.update(set(rxn.reactants)) for rxn in self.children]
        return list(entry_list)




