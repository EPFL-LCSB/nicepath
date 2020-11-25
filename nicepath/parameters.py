class Parameters:
	def __init__(self):
		self.ROOT = ''
		self.PROJECTNAME = ''
		self.COMPOUND = ''

		self.input_file_path = ''
		self.output_folder_path = ''

		self.input_structure = -1

		self.target_compound = ''
		self.precursor_compounds = []
		self.number_of_pathways = -1
		self.min_conserved_ratio = 0
		self.header_list = []
		self.existing_dictionaries = False
		self.exclude_pattern = ''
		self.is_test = False
		self.max_reaction_alternatives = -1
		self.max_num_steps = -1
		self.stop_if_precursor_hit = False
		self.use_parallel = False
		self.folder_name = ''
		self.search_against_model = False
		self.find_branched_pathways = False
		self.only_consider_eligible = False
		self.limit_pathways_length_to_generation = False
		self.target_compounds = []
		self.ignore_infeasible_thermo = False
		self.ignore_infeasible_mechanism = False
		self.ignore_non_model_coreactants = False
		self.use_KEGG_linking_file = False
		self.ignore_novel_reactions = False
		self.transformation_operator = ""
		self.print_compound_sequence = False
		self.molfile_input_path = ""
		self.print_network_statistics = False
		self.ignore_CoA = True
		self.print_bridgit_system = False
		self.search_reverse = False
		self.bigg_model = ''
		self.bigg_compounds = []
		self.load_bigg_model = False
		self.max_diff_score_length = -1
		self.search_against_model = False


	def __str__(self):
		parameter_string = \
			"ROOT|" + self.ROOT + \
			"\nPROJECTNAME|" + self.PROJECTNAME + \
			"\nCOMPOUND|" + self.COMPOUND + \
			"\ninput_file_path|" + self.input_file_path + \
			"\noutput_folder_path|" + self.output_folder_path + \
			"\ninput_structure|" + self.input_structure + \
			"\ntarget_compound|" + self.target_compound + \
			"\nprecursor_compounds|" + str(self.precursor_compounds) + \
			"\nnumber_of_pathways|" + str(self.number_of_pathways) + \
			"\nmin_conserved_ratio|" + str(self.min_conserved_ratio) + \
			"\nheader_list|" + str(self.header_list) + \
			"\nexisting_dictionaries|" + str(self.existing_dictionaries) + \
			"\nexclude_pattern|" + self.exclude_pattern + \
			"\nis_test|" + str(self.is_test) + \
			"\nmax_reaction_alternatives|" + str(self.max_reaction_alternatives) + \
			"\nmax_num_steps|" + str(self.max_num_steps) + \
			"\nstop_if_precursor_hit|" + str(self.stop_if_precursor_hit) + \
			"\nuse_parallel|" + str(self.use_parallel) + \
			"\nfolder_name|" + self.folder_name + \
			"\nfind_branched_pathways|" +  str(self.find_branched_pathways) + \
			"\nonly_consider_eligible|" + str(self.only_consider_eligible) + \
			"\nlimit_pathways_length_to_generation|" + str(self.limit_pathways_length_to_generation) + \
			"\ntarget_compounds|" + str(self.target_compounds) + \
			"\nignore_infeasible_thermo|" + str(self.ignore_infeasible_thermo)  + \
			"\nignore_infeasible_mechanism|" + str(self.ignore_infeasible_mechanism) + \
			"\nuse_KEGG_linking_file|" + str(self.use_KEGG_linking_file) + \
			"\nignore_novel_reactions|" + str(self.ignore_novel_reactions) + \
			"\ntransformation_operator|" + self.transformation_operator + \
			"\nprint_compound_sequence|" + str(self.print_compound_sequence) + \
			"\nmolfile_input_path|" + self.molfile_input_path + \
			"\nprint_network_statistics|" + str(self.print_network_statistics) + \
			"\nignore_CoA|" + str(self.ignore_CoA) + \
			"\nprint_bridgit_system|" + str(self.print_bridgit_system) + \
			"\nsearch_reverse|" + str(self.search_reverse) + \
			"\nbigg_model|" + self.bigg_model + \
			"\nignore_non_model_coreactants|" + str(self.ignore_non_model_coreactants) + \
			"\nmax_diff_score_length|" + str(self.max_diff_score_length)

		return parameter_string



