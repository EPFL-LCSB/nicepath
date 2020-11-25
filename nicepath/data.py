from util import *
import pickle
import os
import re
import time
import multiprocessing
import shutil
import subprocess
from parameters import *
from compound import *
from reaction import *
from pathway import *

class Data:
	def __init__(self, _project_name):
		# name of project
		self.project_name = _project_name

		# containers
		self.Dc = {}        # Compounds: { c_entry : compound instance }
		self.Dr = {}        # Reactions: { c_entry : reaction instance }
		self.Dp = {}        # Pairs with weight: {c_entry1_c_entry2_weight: ['1', '2', ... ] (reaction list) }
		self.G = nx.DiGraph() # Graph structure, directed
		self.Pathways = []  # Pathway list
		self.Dinchikey = {} # Inchikey map: { inchikey: c_entry root }
		self.Drootmap = {}  # Root map : {c_entry: c_entry root }

		self.D_reaction_rootmap = {} # Reaction root map : {r_entry : r_entry root }
		self.D_reaction_hash = {}    # Reaction has map : {reaction hash : r_entry root }

		# files
		self.inputfile_path = ''
		self.outputfolder_path = ''
		self.root_folder = '../'

		# system input
		self.compound_header = ''
		self.reaction_header = ''

		# parameters
		self.parameters = Parameters()
		self.is_iter_object = False
		self.current_precursor = 0

		# bigg database
		self.bigg_db_file = self.root_folder + 'data/Database/BiGG_metabolites.txt'

	def go(self):
		pool = multiprocessing.Pool()
		return pool.map(self.analyzePrecursor, self.parameters.precursor_compounds)


	def readSystemFile(self, _inputfile_path):
		"""
		Reads an BNICE.ch system File to create a compound dictionary, a reaction dictionary, and a weighted Graph
		"""

		# Open file
		try:
			inputfile = open(_inputfile_path)
		except IOError:
			print("Could not open ", _inputfile_path)
			quit()

		# Read file
		reaction = False
		rheader = ''
		cheader = ''
		read_headerline = False
		for line in inputfile:
			new = line.rstrip()
			if new.rstrip() == 'OPERATORS': break

			# category reaction
			if new.rstrip() == 'REACTIONS':
				reaction = True
				self.G.add_nodes_from(self.Dc.keys())
				print('--> Loading Reactions')
				read_headerline = True
				continue

			# check reaction header
			if read_headerline == True and reaction == True:
				self.reaction_header = new.rstrip()
				rheader = self.reaction_header.split(';')
				for x in ['CONSERVEDRATIO', 'RPAIRS', 'ENTRY']:
					assert x in rheader, "Element missing in header: " + x + " - Header: " + new
				if self.parameters.only_consider_eligible :
					assert 'IS_ELIGIBLE' in rheader, "Element missing in header: " + x + " - Header: " + new
				read_headerline = False
				continue

			# check compound header
			if read_headerline == True and reaction == False:
				self.compound_header = new.rstrip()
				cheader = self.compound_header.split(';')
				for x in ['ENTRY']:
					assert x in cheader, "Element missing in header: " + x + " - Header: " + new
				if self.parameters.only_consider_eligible :
					assert 'ELIGIBILITY' in cheader, "Element missing in header: " + x + " - Header: " + new
				if 'INCHIKEY' not in cheader:
					print('WARNING: INCHIKEY not found in header, duplicate compounds will not be merged')
				read_headerline = False
				continue

			# category compounds
			if new.rstrip() == 'COMPOUNDS':
				print('--> Loading Compounds')
				read_headerline = True
				continue

			# load compounds
			if reaction == False:
				new_compound = Compound()
				new_compound.loadCompoundFromLine(new, cheader, self.parameters.use_KEGG_linking_file)
				self.addCompound(new_compound)

			# load reactions
			if reaction == True:
				new_reaction = Reaction()
				new_reaction.loadReactionFromLine(new, rheader, self.parameters.use_KEGG_linking_file)
				self.addReaction(new_reaction)




	def printStatistics(self, _stats_list):
		statfile = open(self.outputfolder_path + 'Statistics_' + self.parameters.target_compound + '.txt', 'w')
		statfile.write('SourceID|TargetID|PathwayDistribution|Total\n')
		overall_total = 0
		for compound in _stats_list:
			distribution = ''
			total = 0
			for length, val in compound['Stats'].items():
				distribution += (str(length)+':'+str(val)+',')
				total += val
			statfile.write('|'.join([compound['Source'], compound['Target'],distribution[:-1],str(total)]) + '\n')
			overall_total += total
		statfile.write('Total||' + str(overall_total) + '\n')
		statfile.close()



	def addPathway(self, pathway):
		# get reactionstrings:
		reactionpaths_list = []
		is_partial_result = False

		l = 1
		reaction_graph = nx.DiGraph()
		prev_nodes = []

		while l < len(pathway):
			pathweight = self.G[pathway[l - 1]][pathway[l]]['weight']
			compounds = [pathway[l - 1], pathway[l]]
			pair_key = '_'.join(compounds) + '_' + str(pathweight)

			if self.parameters.max_reaction_alternatives < len(self.Dp[pair_key]):
				is_partial_result = True

			# check if there is at least one reaction catalyzing the biotransformation in the direction of the pathway
			if self.Dp.get(pair_key):
				reactions_considered = self.Dp[pair_key][:self.parameters.max_reaction_alternatives]
				if prev_nodes == []:
					[reaction_graph.add_edge('source', r) for r in reactions_considered]
				else:
					for node in prev_nodes:
						[reaction_graph.add_edge(node, r) for r in reactions_considered]
				prev_nodes = reactions_considered
			else:
				# if there is a step without valid reaction, the pathway is considered invalid
				return 0
			l += 1

		for node in prev_nodes:
			reaction_graph.add_edge(node,'target')

		try: paths = nx.all_simple_paths(reaction_graph,'source','target', len(pathway)+1)
		except: return 0


		for path in paths:
			reactionpaths_list.append(path[1:-1])

		# create Pathway instances
		for reaction_list in reactionpaths_list:
			new_pathway = Pathway()
			new_pathway.source = pathway[0]
			new_pathway.target = pathway[-1]
			new_pathway.length = len(pathway) - 1

			# remove inconsistent pathways
			if len(reaction_list) != new_pathway.length: # invalid shortcuts found within reaction_graph
				continue
			if len(set(reaction_list)) != len(reaction_list): # duplicate reactions in the pathway
				continue

			# get intermediates
			for c in pathway:
				new_pathway.intermediates.append(c)
			# calculate pathway score
			pathweight = 0
			pathcar = 0
			l = 1
			while l < len(pathway):
				pathweight += self.G[pathway[l - 1]][pathway[l]]['weight']
				pathcar += self.G[pathway[l - 1]][pathway[l]]['car']
				l += 1
			new_pathway.score = round(pathweight,2)
			new_pathway.average_car = round(pathcar/float(new_pathway.length), 2)
			new_pathway.reactions = reaction_list
			new_pathway.is_partial_reactions = is_partial_result

			# get number of novel intermediates
			new_pathway.num_novel_intermediates = 0
			for c in pathway[1:-1]:
				if not self.Dc[c].is_known:
					new_pathway.num_novel_intermediates +=1
			known_count = 0

			for r in reaction_list:
				# get known (KEGG) reaction
				new_pathway.KEGG_ids.append(self.Dr[r].kegg)
				if self.Dr[r].is_known:
					known_count+=1

				# add operators
				new_pathway.operators.append(self.Dr[r].operators)

				# add co_reactants not in model
				if self.parameters.load_bigg_model:
					for entry in self.Dr[r].reactants:
						root = self.Drootmap[entry]
						if root not in self.parameters.bigg_compounds and root not in new_pathway.intermediates and root not in new_pathway.co_reactants:
							new_pathway.co_reactants.append(root)

			new_pathway.num_known_reactions = known_count

			# filters
			# remove pathways with co-reactants not from model
			if self.parameters.ignore_non_model_coreactants and new_pathway.co_reactants != []:
				continue

			# remove compounds with weight-length difference higher then set threshold
			if (new_pathway.score - new_pathway.length) > self.parameters.max_diff_score_length:
				continue

			new_pathway.index = len(self.Pathways)
			new_pathway.checkConsistency()
			self.Pathways.append(new_pathway)
		return

	def writePathways(self, out_file):
		# assemble output line
		if len(self.Pathways) == 0:
			return

		for pathway in self.Pathways:
			out_line = pathway.writeOutputLine(self.parameters.header_list, self.parameters.use_KEGG_linking_file)
			out_file.write(out_line + '\n')


	def printLog(self, _runtime):
		logfile = open(self.outputfolder_path + 'Log.txt', 'w')
		logfile.write('RUN STATISTICS\n--------------------\n')
		logfile.write('Date|' + time.strftime("%d/%m/%Y") + '\n')
		logfile.write('Time|' + time.strftime("%H:%M:%S") + '\n')
		logfile.write('Runtime|' + str(_runtime) + '\n')
		logfile.write('\nPARAMETERS USED\n--------------------\n')
		logfile.write(str(self.parameters))
		logfile.close()

	def prepareOutput(self, out_file_path):
		out_file = open(out_file_path, 'w')
		out_file.write('|'.join(self.parameters.header_list)+ '\n')
		return out_file

	def loadDictionaries(self, root_folder):

		Dr_file = open(root_folder+'Dr.pkl', 'rb')
		self.Dr = pickle.load(Dr_file)
		Dr_file.close()

		Dc_file = open(root_folder+'Dc.pkl', 'rb')
		self.Dc = pickle.load(Dc_file)
		Dc_file.close()

		Dp_file = open(root_folder+'Dp.pkl', 'rb')
		self.Dp = pickle.load(Dp_file)
		Dp_file.close()

		G_file = open(root_folder+'G.pkl', 'rb')
		self.G = pickle.load(G_file)
		G_file.close()

		Dinchikey_file = open(root_folder+'Dinchikey.pkl', 'rb')
		self.Dinchikey = pickle.load(Dinchikey_file)
		Dinchikey_file.close()

		Drootmap_file = open(root_folder+'Drootmap.pkl', 'rb')
		self.Drootmap = pickle.load(Drootmap_file)
		Drootmap_file.close()

		D_reaction_rootmap_file = open(root_folder+'D_reaction_rootmap.pkl', 'rb')
		self.D_reaction_rootmap = pickle.load(D_reaction_rootmap_file)
		D_reaction_rootmap_file.close()

		D_reaction_hash_file = open(root_folder+'D_reaction_hash.pkl', 'rb')
		self.D_reaction_hash = pickle.load(D_reaction_hash_file)
		D_reaction_hash_file.close()

		return True


	def dumpDictionaries(self, root_folder):
		Dr_file = open(root_folder + 'Dr.pkl', 'wb')
		pickle.dump(self.Dr, Dr_file)
		Dc_file = open(root_folder + 'Dc.pkl', 'wb')
		pickle.dump(self.Dc, Dc_file)
		Dp_file = open(root_folder + 'Dp.pkl', 'wb')
		pickle.dump(self.Dp, Dp_file)
		G_file = open(root_folder + 'G.pkl', 'wb')
		pickle.dump(self.G, G_file)
		Dinchikey_file = open(root_folder + 'Dinchikey.pkl', 'wb')
		pickle.dump(self.Dinchikey, Dinchikey_file)
		Drootmap_file = open(root_folder + 'Drootmap.pkl', 'wb')
		pickle.dump(self.Drootmap, Drootmap_file)
		D_reaction_rootmap_file = open(root_folder + 'D_reaction_rootmap.pkl', 'wb')
		pickle.dump(self.D_reaction_rootmap, D_reaction_rootmap_file)
		D_reaction_hash_file = open(root_folder + 'D_reaction_hash.pkl', 'wb')
		pickle.dump(self.D_reaction_hash, D_reaction_hash_file)

	def loadParameters(self, default = False):

		if default:
			parameter_file = open(self.root_folder + 'input/parameters.default')
		else:
			assert os.path.isdir(self.root_folder + 'input/' + self.project_name), 'Invalid project name: %s - Check for typos or create a valid input folder' %(self.root_folder + 'input/' + self.project_name)

			# Read parameter file
			parameter_filename = self.root_folder + 'input/' + self.project_name + '/parameters.txt'
			try: parameter_file = open(parameter_filename)
			except IOError:
				print('Invalid parameter file: ', self.project_name,'. \nPlease specify a valid input.')

		param_dict = {}
		for line in parameter_file:
			if line[0] == '#' or line == '\n':
				continue
			else:
				param_list = line.rstrip().split('|')
				param_dict[param_list[0]] = param_list[1]
		parameter_file.close()


		for param in param_dict: # make this part prettier

			if param == 'precursor_compounds':
				if default:
					continue
				elif param_dict[param][-4:] == '.txt': # load precursors from text file
						precursor_file = open(self.root_folder + 'input/' + self.project_name + '/' + param_dict[param])
						for precursor in precursor_file.readlines():
							self.parameters.precursor_compounds.append(precursor.rstrip())
				elif param_dict[param].isdigit() or (param_dict[param][0]=="C" and param_dict[param][1:].isdigit()): # load only one precursor, indicated in the parameter file
					self.parameters.precursor_compounds.append(param_dict[param])
				elif param_dict[param] == 'model':
					self.parameters.search_against_model = True
				else:
					print('ERROR: Precursor compound format not recognized. Please define "precursor_compounds" correctly')
					quit()

			elif param == 'target_compounds':
				if default:
					continue
				elif param_dict[param][-4:] == '.txt': # load targets from text file
						target_file = open(self.root_folder + 'input/' + self.project_name + '/' + param_dict[param])
						for target in target_file.readlines():
							self.parameters.target_compounds.append(target.rstrip())
				else : # load only one target, indicated in the parameter file
					assert int(param_dict[param]), "PARAMETER ERROR:" +  param_dict[param] + "is not a number"
					self.parameters.target_compounds.append(param_dict[param])

			elif param == 'number_of_pathways':
				self.parameters.number_of_pathways = int(param_dict[param])
			elif param == 'min_conserved_ratio':
				self.parameters.min_conserved_ratio = float(param_dict[param])
			elif param == 'header_list':
				self.parameters.header_list = param_dict[param].split(',')
			elif param == 'compound_sequence_only':
				self.parameters.compound_sequence_only = getBoolean(param_dict[param])
			elif param == 'existing_dictionaries':
				self.parameters.existing_dictionaries = getBoolean(param_dict[param])
			elif param == 'exclude_pattern':
				self.parameters.exclude_pattern = param_dict[param]  # FIX THIS to have regexp in file
			elif param == 'is_test':
				self.parameters.is_test = getBoolean(param_dict[param])
			elif param == 'max_reaction_alternatives':
				self.parameters.max_reaction_alternatives = int(param_dict[param])
			elif param == 'max_num_steps':
				self.parameters.max_num_steps = int(param_dict[param])
			elif param == 'stop_if_precursor_hit':
				self.parameters.stop_if_precursor_hit = getBoolean(param_dict[param])
			elif param == 'use_parallel':
				self.parameters.use_parallel = getBoolean(param_dict[param])
			elif param == 'folder_name':
				self.parameters.folder_name = param_dict[param]
			elif param == 'only_consider_eligible':
				self.parameters.only_consider_eligible = getBoolean(param_dict[param])
			elif param == 'limit_pathways_length_to_generation':
				self.parameters.limit_pathways_length_to_generation = getBoolean(param_dict[param])
			elif param == 'ignore_infeasible_thermo':
				self.parameters.ignore_infeasible_thermo = getBoolean(param_dict[param])
			elif param == 'ignore_infeasible_mechanism':
				self.parameters.ignore_infeasible_mechanism = getBoolean(param_dict[param])
			elif param == 'ROOT':
				self.parameters.ROOT = param_dict[param]
				if self.parameters.ROOT != '':
					assert self.parameters.ROOT[-1] == '/', 'Parameter issue: ROOT should end with "/"'
			elif param == 'PROJECTNAME':
				self.parameters.PROJECTNAME = param_dict[param]
			elif param == 'COMPOUND':
				self.parameters.COMPOUND = param_dict[param]
			elif param == 'find_branched_pathways':
				self.parameters.find_branched_pathways = getBoolean(param_dict[param])
			elif param == 'input_structure':
				assert param_dict[param] in ['1', '2'], 'Parameter issue: invalid parameter "input_structure" '
				self.parameters.input_structure = param_dict[param]
			elif param == 'input_file_path':
				self.inputfile_path = param_dict[param]
			elif param == 'use_KEGG_linking_file' :
				self.parameters.use_KEGG_linking_file = getBoolean(param_dict[param])
			elif param == 'ignore_novel_reactions' :
				self.parameters.ignore_novel_reactions = getBoolean(param_dict[param])
			elif param == 'transformation_operator' :
				self.parameters.transformation_operator = param_dict[param]
			elif param == 'print_compound_sequence' :
				self.parameters.print_compound_sequence = getBoolean(param_dict[param])
			elif param == 'print_network_statistics' :
				self.parameters.print_network_statistics = getBoolean(param_dict[param])
			elif param == 'molfile_input_path' :
				self.parameters.molfile_input_path = param_dict[param]
			elif param == 'ignore_CoA':
				self.parameters.ignore_CoA = getBoolean(param_dict[param])
			elif param == 'print_bridgIT_system':
				self.parameters.print_bridgit_system = getBoolean(param_dict[param])
			elif param == 'search_reverse':
				self.parameters.search_reverse = getBoolean(param_dict[param])
			elif param == 'bigg_model':
				self.parameters.bigg_model = param_dict[param]
				self.parameters.load_bigg_model = True
			elif param == 'ignore_non_model_coreactants':
				self.parameters.ignore_non_model_coreactants = getBoolean(param_dict[param])
			elif param == 'max_diff_score_length':
				self.parameters.max_diff_score_length = float(param_dict[param])
			elif param == 'output_folder_path':
				self.outputfolder_path = param_dict[param]
				if self.outputfolder_path != '':
					assert self.outputfolder_path[-1] == '/', 'Parameter issue: output_folder_path should end with "/"'
			else:
				print('WARNING: Parameter ' + param + ' is not known!!! Check parameter.txt for typos...')

		if self.parameters.load_bigg_model:
			assert self.parameters.bigg_model != '', "ERROR: no bigg model specified. Please define parameter 'bigg_model'"
			self.loadKeggIDsFromBigg(self.parameters.bigg_model)

		# For option 1 :
		if self.parameters.input_structure == '1' and not default:

			# input
			self.inputfile_path = self.parameters.ROOT  + self.parameters.PROJECTNAME + '/' + self.parameters.COMPOUND + '/systemFile_' + self.parameters.COMPOUND + '.txt'
			if self.parameters.is_test: self.inputfile_path = self.parameters.ROOT + self.parameters.PROJECTNAME + '/' + self.parameters.COMPOUND + '/systemFile_' + self.parameters.COMPOUND + '_test.txt'

			# appendix
			APPENDIX = '_' + str(self.parameters.number_of_pathways)
			if self.parameters.folder_name != '':
				APPENDIX += '_' + self.parameters.folder_name
			if self.parameters.is_test:
				APPENDIX += '_test'
			if self.parameters.exclude_pattern != '':
				APPENDIX += ('_exclude-' + self.parameters.exclude_pattern)

			# output
			self.outputfolder_path = self.parameters.ROOT  +  self.parameters.PROJECTNAME + '/' + self.parameters.COMPOUND + '/PathwaySearch' + APPENDIX + '/'

		# Make sure output folder path does not point to ../data
		if not default:
			assert '../data/' not in self.outputfolder_path, 'Parameter issue: invalid output folder. Please specify an output folder outside "../data"'
			assert '/nicepath/data/' not in self.outputfolder_path, 'Parameter issue: invalid output folder. Please specify an output folder outside "../data"'

			print('Parameter file loaded successfully for project '+self.project_name)
			print('------------------------------------------------------------------')
			# print('Parameters:')
			# print(self.parameters)
			print('Input file path: '+self.inputfile_path)
			print('Output folder path: '+self.outputfolder_path)
			print('------------------------------------------------------------------')

		return True



	def analyzePrecursor(self, _precursor_compound):
		'''In case of tautomers, the precursor compound should be the first in the system file'''
		if self.parameters.search_reverse:
			start_compound = self.parameters.target_compound
			end_compound = _precursor_compound
		else:
			start_compound = _precursor_compound
			end_compound = self.parameters.target_compound

		pathwayStatistics = {}

		# Do pathway search
		pathway_count_old = len(self.Pathways)
		max_number_of_steps = self.parameters.max_num_steps

		try:
			k_shortest = k_shortest_paths(self.G, start_compound, end_compound,
								 self.parameters.number_of_pathways, 'weight')
		except:
			k_shortest = []

		for pathway in k_shortest:
			# Ignore pathways longer than the maximum indicated
			if self.parameters.limit_pathways_length_to_generation:
				max_number_of_steps = self.Dc[end_compound].generation

			if (len(pathway) - 1) > max_number_of_steps:
				continue

			# Stop if we hit a precursor in the pathway
			if self.parameters.stop_if_precursor_hit == True:
				intermediate_found = False
				# check the pathway for intermediates, but remove the target
				if self.parameters.search_reverse:
					pathway_without_target = pathway[:-1]
				else:
					pathway_without_target = pathway[1:]
				i = 0
				while intermediate_found == False and i < len(pathway_without_target):
					if pathway_without_target[i] in self.parameters.bigg_compounds:
						intermediate_found = True
					i += 1
				if intermediate_found == True:
					continue

			self.addPathway(pathway)

		diff = len(self.Pathways) - pathway_count_old
		if diff != 0:
			pathwayStatistics = collectStatistics(pathwayStatistics, pathway, diff)

		# print output
		print('Found %d pathways for %s -> %s' % (diff, start_compound, end_compound))
		if len(self.Pathways) != 0  and self.parameters.use_parallel:
			outputfile_path = self.outputfolder_path + 'Pathways_' + start_compound + '_' + end_compound + '.txt'
			outputfile = self.prepareOutput(outputfile_path)
			# write pathways to file
			self.writePathways(outputfile)
			outputfile.close()


		# return stattistics
		return {'Source': start_compound, 'Target': end_compound, 'Stats': pathwayStatistics}



	def loadData(self):
		# Create output folder (if not existing)
		if not os.path.isdir(self.outputfolder_path):
			os.mkdir(self.outputfolder_path)

		# Read Data from systemFile or existing dictionary
		if self.parameters.existing_dictionaries == True:
			# load Dp, Dr, Dc
			print('loading Data...')
			self.loadDictionaries(self.outputfolder_path)
			nx.freeze(self.G)
		else:
			# read systemFile
			print('Reading System File...')
			self.readSystemFile(self.inputfile_path)

		# If search against model, we convert KEGG ids in precursor list to entries
		if self.parameters.bigg_model != '':
			entries = []
			for compound in self.Dc.keys():
				if self.Dc[compound].kegg in self.parameters.bigg_compounds:
					entries.append(compound)
			self.parameters.bigg_compounds = entries

			if self.parameters.search_against_model:
				self.parameters.precursor_compounds = entries

		# Statistics output
		print('\nSystem file successfully loaded')
		print('------------------------------------------------------------------')
		print('Compounds loaded: ' + str(len(self.Dc)))
		print('Unique compounds loaded: ', len(self.Dinchikey.keys()))
		print('Reactions loaded: ' + str(len(self.Dr)))
		print('Edge library loaded: ' + str(len(self.Dp)))
		print('Graph loaded: ' + str(len(self.G.nodes())) + ' nodes, ' + str(len(self.G.edges())) + ' edges')
		print('Precursors loaded: ' + str(len( self.parameters.precursor_compounds)))


	def saveDictionaries(self):
		# Save data to dictionaries Dr, Dp, Dc and print network file for pairs
		print('Saving data to dictionaries...')
		if self.parameters.existing_dictionaries == False:
			pairout = open(self.outputfolder_path + 'Network_pairs.csv', 'w')
			pairout.write("Source,Target,Distance,CAR,Number_of_reactions\n")
			for pair in self.Dp.keys():
				pairout.write(','.join(pair.split('_')) + ',' + str(self.G[pair.split('_')[0]][pair.split('_')[1]]['car']) + ',' + str(len(self.Dp[pair])) + '\n')

			self.dumpDictionaries(self.outputfolder_path)
			pairout.close()

	def loadKeggIDsFromBigg(self, _model_name):
		bigg_ids = open(self.bigg_db_file)
		headerlist  = bigg_ids.readline().rstrip().split('\t')
		model_list_index = headerlist.index('model_list')
		db_link_index = headerlist.index('database_links')
		for line in bigg_ids:
			linelist = line.rstrip().split('\t')
			if _model_name in linelist[model_list_index]:
				match = re.findall(r'C\d\d\d\d\d', linelist[db_link_index])
				for kegg_id in match:
					if kegg_id not in self.parameters.bigg_compounds:
						self.parameters.bigg_compounds.append(kegg_id)

	def getInputMolfileFolder(self):
		input_folder = ''
		if os.path.isdir(self.parameters.molfile_input_path):
			input_folder = self.parameters.molfile_input_path
		elif os.path.isdir(
				self.parameters.ROOT + self.parameters.PROJECTNAME + '/' + self.parameters.COMPOUND + '/molfiles/'):
			input_folder = self.parameters.ROOT + self.parameters.PROJECTNAME + '/' + self.parameters.COMPOUND + '/molfiles/'
		else:
			print('Please specifiy input molfile folder, either by using INPUT OPTION 1, or by directly setting molfile_input_path')
			quit()
		return input_folder

	def createOutputFolder(self, _folder_name):
		output_folder = self.outputfolder_path + _folder_name+ '/'
		if os.path.isdir(output_folder):
			shutil.rmtree(output_folder)
		os.mkdir(output_folder)
		output_folder = self.outputfolder_path + _folder_name + '/'
		return output_folder

	def writeCompoundSequences(self):
		if not self.parameters.print_compound_sequence:
			return
		print('Writing compound sequence molfiles')
		input_folder = self.getInputMolfileFolder()
		output_folder = self.createOutputFolder('compound_sequences')
		for pathway in self.Pathways:
			os.mkdir(output_folder + str(pathway.index)+'/')
			for i, compound in enumerate(pathway.intermediates):
				new_molfile = open(output_folder + str(pathway.index) + '/' + str(i) + '.mol', 'w')
				input_compound = open(input_folder + compound + '.mol')
				molfile = input_compound.read()
				new_molfile.write(molfile)
				new_molfile.close()


	def printNetworkStatistics(self):
		print('Calculating network statistics...')
		stats_file = open(self.outputfolder_path + 'Network_statistics.txt', 'w')
		stats_file.write('ORIGINAL (WEIGHTED) GRAPH\n--------------------\n')
		stats_file.write('Number of nodes,%d\n'%self.G.number_of_nodes())
		stats_file.write('Number of edges,%d\n' % self.G.number_of_edges())
		stats_file.write('Number of disconnected nodes,%d\n' % len(list(nx.isolates(self.G))))
		stats_file.write('Number of connected nodes,%d\n' % (self.G.number_of_nodes() - len(list(nx.isolates(self.G)))))
		undirected_G = self.G.to_undirected()
		undirected_main_G = max(nx.connected_component_subgraphs(undirected_G), key=len)
		stats_file.write('Network diameter,%d\n' % nx.diameter(undirected_main_G))
		#stats_file.write('Average shortest path length,%d\n' % nx.average_shortest_path_length(undirected_main_G))
		#stats_file.write('Average shortest path length weight,%f\n' % nx.average_shortest_path_length(undirected_main_G, weight = 'weight'))


		CAR_cutoff = 0.34 # Best predictor for KEGG RPAIR of type "main"
		stats_file.write('\nUNWEIGHTED GRAPH - CAR THRESHOLD = %s\n--------------------\n'%(str(CAR_cutoff)))
		#create new unweighted graph with CAR >= 0.34
		simple_G = nx.Graph()
		for (u, v, c) in self.G.edges.data('car'):
			if c > CAR_cutoff:
				simple_G.add_edge(u, v)

		stats_file.write('Number of nodes,%d\n' % simple_G.number_of_nodes())
		stats_file.write('Number of edges,%d\n' % simple_G.number_of_edges())
		stats_file.write('Number of disconnected nodes,%d\n' % len(list(nx.isolates(simple_G))))
		stats_file.write('Number of components,%d\n' % nx.number_connected_components(simple_G))

		# Statistics on main / giant component
		stats_file.write('\nUNWEIGHTED GRAPH - CAR THRESHOLD = %s - LARGEST COMPONENT\n--------------------\n' % (str(CAR_cutoff)))
		self.main_component_G = max(nx.connected_component_subgraphs(simple_G), key=len)
		stats_file.write('Number of nodes,%d\n' % self.main_component_G.number_of_nodes())
		stats_file.write('Number of edges,%d\n' % self.main_component_G.number_of_edges())
		stats_file.write('Network diameter,%d\n' % nx.diameter(self.main_component_G))
		stats_file.write('Average shortest path length,%d\n' % nx.average_shortest_path_length(self.main_component_G))
		stats_file.write('Percentage nodes,%f\n' % (self.main_component_G.number_of_nodes() / simple_G.number_of_nodes() * 100))
		stats_file.write('Percentage edges,%f\n' % (self.main_component_G.number_of_edges() / simple_G.number_of_edges() * 100))
		#stats_file.write('Average degree connectivity,%f\n' % nx.average_degree_connectivity(self.main_component_G))
		stats_file.close()



	def addCompound(self, new_compound):
		''' Adds compound to compound dictionary (Dc)
			If there's an inchikey, also adds inchikey to inchikey map
			Remove duplicates based on inchikey '''

		if new_compound.inchikey != '' and self.Dinchikey.get(new_compound.inchikey): # check if compound already exists
			entry = self.Dinchikey[new_compound.inchikey]
			self.Dc[entry].addChild(new_compound)
			self.Drootmap[new_compound.entry] = entry
			# merge KEGG identifiers
			if self.Dc[entry].kegg == '':
				self.Dc[entry].kegg = new_compound.kegg

		else: # no duplicate found, new compound entry in Dc
			self.Dc[new_compound.entry] = new_compound
			self.Drootmap[new_compound.entry] = new_compound.entry
			self.Dinchikey[new_compound.inchikey] = new_compound.entry


	def addReaction(self, new_reaction):

		# check reaction
		if self.parameters.exclude_pattern != '':
			if re.search(self.parameters.exclude_pattern, new_reaction.equation_name):
				return
		if self.parameters.only_consider_eligible:
			if not new_reaction.eligible:
				return
		if self.parameters.ignore_novel_reactions:
			if not new_reaction.is_known:
				return

		# check directionality
		forward = True
		reverse = True
		if self.parameters.ignore_infeasible_thermo:
			new_reaction.determineThermoFeasibility()
			if not new_reaction.isThermoFeasible('forward'):
				forward = False
			if not new_reaction.isThermoFeasible('reverse'):
				reverse = False
		if self.parameters.ignore_infeasible_mechanism:
			new_reaction.determineMechanismFeasibility()
			if not new_reaction.isMechanismFeasible('forward'):
				forward = False
			if not new_reaction.isMechanismFeasible('reverse'):
				reverse = False

		add_reaction = False
		if self.parameters.use_KEGG_linking_file:
			add_reaction = True
		else:
			# generate reaction hash
			rhash = new_reaction.generateReactionHash(self.Drootmap)

			# check if hash already there
			if rhash != '' and self.D_reaction_hash.get(rhash):
				entry = self.D_reaction_hash[rhash]
				self.Dr[entry].addChild(new_reaction)
				self.D_reaction_rootmap[new_reaction.entry] = entry

			else: # no duplicate found, new compound entry in Dr

				self.D_reaction_rootmap[new_reaction.entry] = new_reaction.entry
				self.D_reaction_hash[rhash] = new_reaction.entry
				add_reaction = True

		if add_reaction:
			self.Dr[new_reaction.entry] = new_reaction
			# add pairs to graph
			for i, p in enumerate(new_reaction.rpairs):
				left_root = self.Drootmap[p[0]]
				right_root = self.Drootmap[p[1]]

				if self.parameters.min_conserved_ratio >= new_reaction.cars[i]:
					continue

				if self.parameters.only_consider_eligible:
					if self.Dc[left_root].eligible is False or self.Dc[right_root].eligible is False:
						continue

				if self.parameters.ignore_CoA:
					if self.Dc[left_root].is_coa is True or self.Dc[right_root].is_coa is True:
						continue

				conservedRatio = round(new_reaction.cars[i], 2)
				dist = calculateDistance(conservedRatio, self.parameters.transformation_operator)
				dist = round(dist, 2)
				if forward:
					self.G.add_edge(left_root, right_root, weight=dist, car=conservedRatio)
				if reverse:
					self.G.add_edge(right_root, left_root, weight=dist, car=conservedRatio)

				pair_forward = left_root + '_' + right_root + '_' + str(dist)
				pair_reverse = right_root + '_' + left_root + '_' + str(dist)

				if forward:
					if not self.Dp.get(pair_forward):
						self.Dp[pair_forward] = [new_reaction.entry]
					else:
						if self.Dr[new_reaction.entry].kegg != '':  # make sure that KEGG reactions are in the beginning of the list
							self.Dp[pair_forward].insert(0, new_reaction.entry)
						else:
							self.Dp[pair_forward].append(new_reaction.entry)
				if reverse:
					if not self.Dp.get(pair_reverse):
						self.Dp[pair_reverse] = [new_reaction.entry]
					else:
						if self.Dr[new_reaction.entry].kegg != '':  # make sure that KEGG reactions are in the beginning of the list
							self.Dp[pair_reverse].insert(0, new_reaction.entry)
						else:
							self.Dp[pair_reverse].append(new_reaction.entry)

	def checkQueryCompounds(self):
		#Check targets
		for i, target in enumerate(self.parameters.target_compounds):
			if not self.Drootmap.get(target):
				print('ERROR: compound %s not found in input system'%target)
				exit(1)
			elif target != self.Drootmap[target]:
				self.parameters.target_compounds[i] = self.Drootmap[target]
				print('WARNING: target compound %s has been replaced by root compound %s'%(target, self.Drootmap[target]))

			# Check if targets in main_component
			if self.parameters.print_network_statistics and self.Drootmap[target] not in self.main_component_G.nodes():
				print('WARNING: target compound %s may be disconnected from the main component'%self.Drootmap[target])

		# Check precursors
		for i, prec in enumerate(self.parameters.precursor_compounds):
			if not self.Drootmap.get(prec):
				print('ERROR: compound %s not found in input system' %prec)
				exit(1)
			elif prec != self.Drootmap[prec]:
				self.parameters.precursor_compounds[i] = self.Drootmap[prec]
				print('WARNING: precursor compound %s has been replaced by root compound %s' % (prec, self.Drootmap[prec]))

			# Check if precursor in main_component
			if self.parameters.print_network_statistics and self.Drootmap[prec] not in self.main_component_G.nodes():
				print('WARNING: precursor compound %s may be disconnected from the main component' %self.Drootmap[prec])

	def printBridgITSystem(self, compound_entry):
		''' Printing a cropped BridgIT system, only containing compounds and reactions from the pathway
			It also creates a molfile folder.
			Only supported for input option 1 (but can be adapted for option 2) '''
		if not self.parameters.print_bridgit_system:
			return
		if len(self.Pathways) == 0:
			return
		print("Printing BridgIT input system")
		source_molfile_folder = self.getInputMolfileFolder()
		bridgit_folder_name = 'BridgIt_Input_' + compound_entry
		output_folder = self.createOutputFolder(bridgit_folder_name)
		bridgit_file = output_folder+ "systemFile_BridgIt.txt"
		target_molfile_folder = self.createOutputFolder(bridgit_folder_name + '/molfiles')

		# Find all reactions and associated compounds used in the pathways
		for p in self.Pathways:
			for r in p.reactions:
				# mark reaction as used
				self.Dr[r].markAsUsed()
				# mark compound as used
				[self.Dc[self.Drootmap[cpd]].markAsUsed() for cpd in self.Dr[r].getAllDependentCompounds()]
				# mark reaction children as used
				for rchild in self.Dr[r].children:
					rchild.markAsUsed()
					[self.Dc[self.Drootmap[cpd]].markAsUsed() for cpd in rchild.getAllDependentCompounds()]

				# mark reaction children compounds as used



		all_compounds = self.Dc.keys()
		all_reactions = self.Dr.keys()
		#sorted_keys = sorted(all_compounds) # doesn't work well with strings.. oops

		# move molfiles
		for c in all_compounds:
			if self.Dc[c].is_in_pathway:
				molfile = source_molfile_folder + str(c) + '.mol'
				molfile2 = target_molfile_folder + str(c) + '.mol'
				subprocess.call(['cp', molfile, molfile2])
				for child in self.Dc[c].children:
					molfile = source_molfile_folder + str(child.entry) + '.mol'
					molfile2 = target_molfile_folder + str(child.entry) + '.mol'
					subprocess.call(['cp', molfile, molfile2])

		# write compounds section
		bridgit = open(bridgit_file, 'w')
		bridgit.write('COMPOUNDS\n' + self.compound_header + ';ROOT_COMPOUND'+ '\n')
		for c in all_compounds:
			if self.Dc[c].is_in_pathway:
				bridgit.write(self.Dc[c].fileline + ';' + self.Dc[c].entry + '\n')
				for child in self.Dc[c].children:
					bridgit.write(child.fileline + ';' + self.Dc[c].entry + '\n')

		# write reactions section
		bridgit.write('REACTIONS\n' + self.reaction_header + ';ROOT_REACTION;REACTION_HASH'+ '\n')
		for r in all_reactions:
			if self.Dr[r].is_in_pathway:
				bridgit.write(self.Dr[r].fileline + ';' + self.Dr[r].entry + ';' + self.Dr[r].reaction_hash + '\n')
				for child in self.Dr[r].children:
					bridgit.write(child.fileline + ';' + self.Dr[r].entry  + ';' + self.Dr[r].reaction_hash + '\n')

		return 0
