import re

class Compound:
    def __init__(self):
        self.name = ''
        self.eligible = False
        self.formula = ''
        self.old_entries = []
        self.entry = ''
        self.inchikey = ''
        self.database_links = ''
        self.kegg = ''
        self.fileline = ''
        self.generation = -1
        self.is_coa = False
        self.children = [] # duplicate compounds
        self.is_in_pathway = False
        self.is_known = False



    def loadCompoundFromLine(self, _line, _headerlist, _is_kegg_linking):
        linelist = _line.rstrip().split(';')
        self.fileline = _line
        if 'FORMULA' in _headerlist: self.formula = linelist[_headerlist.index('FORMULA')]
        if 'ELIGIBILITY' in _headerlist: self.eligible = bool(int(linelist[_headerlist.index('ELIGIBILITY')]))
        if 'INCHIKEY' in _headerlist: self.inchikey = linelist[_headerlist.index('INCHIKEY')].split('-')[0] # only first part of inchikey
        if 'DATABASE_LINKS' in _headerlist: self.database_links = linelist[_headerlist.index('DATABASE_LINKS')] # not parsed
        if 'KEGG' in _headerlist: self.kegg = linelist[_headerlist.index('KEGG')][:6] # tautomer information is removed to ensure correct lookup in model
        if 'ENTRY' in _headerlist: self.entry = linelist[_headerlist.index('ENTRY')]
        if 'NAME' in _headerlist: self.name = linelist[_headerlist.index('NAME')]
        if 'GENERATION' in _headerlist: self.generation = int(linelist[_headerlist.index('GENERATION')])
        self.checkForCoa()
        if _is_kegg_linking:
            self.is_known = True
            self.kegg = self.entry
        elif self.kegg != '':
            self.is_known = True

    # check this function, it will be used for extended tautomers output
    def mergeCompounds(self, other):
        self.old_entries.append(other.entry)
        # get name/ID info from KEGG compounds
        if other.kegg != '':
            self.other_names.append(self.name)
            self.kegg = other.kegg
            self.name = other.name
        if self.name == '':
            self.name = other.name
        if self.database_links != other.database_links:
            self.database_links += ('|'+other.database_links)
        if other.eligible: self.eligible = True

    def checkForCoa(self):
        if self.kegg == 'C00010':
            self.is_coa = True
        elif self.entry  == 'C00010':
            self.is_coa = True
        elif self.name == 'CoA' or self.name == 'coa' or self.name  == 'Coenzyme A':
            self.is_coa = True

    def addChild(self, other):
        self.children.append(other)

    def markAsUsed(self):
        self.is_in_pathway = True
        for child in self.children:
            child.is_in_pathway = True