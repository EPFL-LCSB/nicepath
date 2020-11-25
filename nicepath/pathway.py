class Pathway:
    def __init__(self):
        self.index = -1
        self.source = ""
        self.target = ""
        self.length = -1
        self.intermediates = []
        self.score = -1
        self.average_car = -1
        self.reactions = []
        self.is_partial_reactions = False
        self.num_novel_intermediates = -1
        self.KEGG_ids = []
        self.num_known_reactions = -1
        self.operators = []
        self.co_reactants = [] # co-substrates, or co-products in the case of reverse search, that are not in the bigg_model

    def writeOutputLine(self, header, is_KEGG_linking):
        outlist = []
        for category in header:
            if category == 'pathwayIndex':
                outlist.append(str(self.index))
            elif category == 'Source':
                outlist.append(self.source)
            elif category == 'Target':
                outlist.append(self.target)
            elif category == 'Pathway_Length':
                outlist.append(str(self.length))
            elif category == 'Intermediates':
                if not is_KEGG_linking:
                    c_intermediates = []
                    [c_intermediates.append('C'+i) for i in self.intermediates]
                else:
                    c_intermediates = self.intermediates
                outlist.append(';'.join(c_intermediates))
            elif category == 'Pathway_score':
                outlist.append(str(self.score))
            elif category == 'Average_CAR':
                outlist.append(str(self.average_car))
            elif category == 'Reaction_Entries':
                if not is_KEGG_linking:
                    reactions = []
                    [reactions.append('R'+i) for i in self.reactions]
                else:
                    reactions = self.reactions
                outlist.append(';'.join(reactions))
            elif category == 'Partial_reactions':
                outlist.append(str(self.is_partial_reactions))
            elif category == 'NumNovelIntermediates':
                outlist.append(str(self.num_novel_intermediates))
            elif category == 'Reaction_KEGGID':
                outlist.append(';'.join(self.KEGG_ids))
            elif category == 'Number_of_known_steps_(reactions)':
                outlist.append(str(self.num_known_reactions))
            elif category == 'Reaction_Operators':
                outlist.append(';'.join(self.operators))
            elif category == 'Co_reactants':
                outlist.append(';'.join(self.co_reactants))
            else:
                outlist.append('')

        return '|'.join(outlist)

    def checkConsistency(self):
        assert self.length == len(self.reactions), "Problem in reaction %s: reactions and pathway length are not consistent" %(self.index)
        assert self.length == len(self.intermediates)-1, "Problem in reaction %s: intermediates and pathway length are not consistent" % (self.index)
