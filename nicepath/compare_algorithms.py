# This script evaluates the computational efficacy between
# 1) searching a weighted network using Yen's k-shortest path search algorithm (as in NICEpath) and
# 2) searching a KEGG "main" RPAIR network using a breadth-first search algorithm (simplified PathPred concept)

import networkx as nx
from util import *
import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def __main__():
    # Define filename
    KEGG_file_path = "../data/Databases/KEGG/Main_RPAIRS_KEGG.tsv"
    # Load KEGG networks
    G_CAR, G_RPAIR, G_NOCOFAC = load_Graphs(KEGG_file_path)
    # Reference pathways
    evaluation_file_path = "../data/Databases/KEGG/Pathways_for_evaluation.txt"
    output_folder_path = "../output/KEGG/"
    data_file_path = evaluate_pathways(G_CAR, G_RPAIR, G_NOCOFAC, evaluation_file_path, output_folder_path)
    create_figures(output_folder_path, 'Pathway_evaluation_table.csv')


def load_Graphs(file_path):
    '''Load simple graph from KEGG RPAIR main, and weighted graph from CAR
        G_RPAIR: Unweighted RPAIR graph
        G_CAR: Weighted reactant pair graph
    '''

    # Create network containers
    G_RPAIR = nx.Graph()
    G_CAR = nx.Graph()
    G_NOCOFAC = nx.Graph()
    # Open file
    KEGG_file = open(file_path, 'r')
    # Parse file
    header = KEGG_file.readline().rstrip().split()
    for line in KEGG_file:
        source = line.split()[header.index("Reactant_pair")].split('_')[0]
        target = line.split()[header.index("Reactant_pair")].split('_')[1]
        CAR = float(line.split()[header.index("CAR")])
        is_main = getBoolean(line.split()[header.index("RPAIR_main")])

        if is_main: G_RPAIR.add_edge(source, target) # add main pairs to G_RPAIR graph

        G_CAR.add_edge(source, target,
                       dflt=calculateDistance(CAR, ''),
                       sqrt=calculateDistance(CAR, 'sqrt'),
                       exp=calculateDistance(CAR,'exp'),
                       car=CAR)
        G_NOCOFAC.add_edge(source, target)

    # remove coA from network
    G_CAR.remove_node('C00010')
    G_RPAIR.remove_node('C00010')

    # remove standard cofactors from network
    cofactor_list = ['C00003', 'C00004', 'C00080', 'C00006', 'C00005', 'C00011', 'C00001', 'C00009', 'C00010', 'C00007', 'C00399', 'C00390', 'C00027', 'C00016', 'C01352', 'C00014', 'C00708', 'C00272', 'C01847', 'C00061', 'C00343', 'C00342', 'C00019', 'C00021', 'C04488', 'C01217', 'C00035', 'C00015', 'C00013', 'C01137', 'C00002', 'C00008', 'C00020', 'C00044']
    G_NOCOFAC.remove_nodes_from(cofactor_list)


    return G_CAR, G_RPAIR, G_NOCOFAC

def get_unweighted_graph(G):
    CAR_cutoff = 0.34  # Best predictor for KEGG RPAIR of type "main"
    # create new unweighted graph with CAR >= 0.34
    simple_G = nx.Graph()
    for (u, v, c) in G.edges.data('car'):
        if c > CAR_cutoff:
            simple_G.add_edge(u, v)
    return simple_G

def evaluate_pathways(G_CAR, G_RPAIR, G_NOCOFAC, ref_file_path, out_folder_path):
    # create a graph from G_CAR only keeping edges above threshold
    # according to previous analysis, should be similar to G_RPAIR
    G_CAR_unweighted = get_unweighted_graph(G_CAR)
    ref_file = open(ref_file_path)
    all_results = []
    count = 0
    # evaluate pathways one by one
    for line in ref_file:
        count += 1
        print('Evaluating pathway number',count)
        output = evaluate_reference_pathway(G_CAR, G_RPAIR, G_NOCOFAC, G_CAR_unweighted, line.rstrip().split(), count)
        all_results.append(output)

    # Save output to data table
    results_df = pd.DataFrame(all_results)
    out_file_path = out_folder_path + 'Pathway_evaluation_table.csv'
    results_df.to_csv(out_file_path)
    return out_file_path


def evaluate_reference_pathway(G_CAR, G_RPAIR, G_NOCOFAC, G_CAR_unweighted, reference_pathway, pathway_number):
    # define search parameters
    source = reference_pathway[0]
    target = reference_pathway[-1]
    k = 100
    evaluation_cases = ['CAR_exp', 'CAR_dflt', 'CAR_sqrt', 'CAR_034', 'KEGG_RPAIR', 'KEGG_nocofac']
    results = {"Entry": pathway_number, 'Length': len(reference_pathway)-1}

    for case in evaluation_cases:
        start_time = time.time()
        try: #try to run pathway search
            # if no weight is given, networkx does BFS from both source and target and meets in the middle
            # BFS: Breadth-first search
            if case == 'KEGG_RPAIR':
                paths = k_shortest_paths(G_RPAIR, source, target, k, 'None')

            # run Yen's k-shortest loopless paths algorithm on differently weighted network
            elif case == 'CAR_dflt':
                paths = k_shortest_paths(G_CAR, source, target, k, 'dflt')
            elif case == 'CAR_exp':
                paths = k_shortest_paths(G_CAR, source, target, k, 'exp')
            elif case == 'CAR_sqrt':
                paths = k_shortest_paths(G_CAR, source, target, k, 'sqrt')
            elif case == 'CAR_034':
                paths = k_shortest_paths(G_CAR_unweighted, source, target, k, 'None')
            else:  # case == 'KEGG_nocofac'
                paths = k_shortest_paths(G_NOCOFAC, source, target, k, 'None')
        except:
            paths = []
        runtime = time.time() - start_time

        if reference_pathway in paths:
            rank = str(paths.index(reference_pathway) + 1)
        elif paths != []:
            rank = 'Not found'
            runtime = 'NA'
        else:
            rank = 'Not possible'
            runtime = 'NA'

        results[case+'_rank'] = rank
        if runtime != 'NA':
            results[case+'_runtime'] = round(runtime, 6)
        else: results[case+'_runtime'] = runtime

    # Check if all necessary edges are in graphs
    comment = ''
    for index, compound in enumerate(reference_pathway[:-1]):
        source = compound
        target = reference_pathway[index+1]
        if not G_CAR.has_edge(source, target):
            comment += ('%s_%s not in G_CAR, '%(source, target))
        if not G_RPAIR.has_edge(source, target):
            comment += ('%s_%s not in G_RPAIR, ' % (source, target))
        if not G_NOCOFAC.has_edge(source, target):
            comment += ('%s_%s not in G_NOCOFAC, ' % (source, target))
        if not G_CAR_unweighted.has_edge(source, target):
            comment += ('%s_%s not in G_CAR_unweighted, ' % (source, target))



    results['comment'] = comment

    return results

def create_figures(output_folder_path, results_path):
    data_file_path = output_folder_path + results_path
    # read data
    df = pd.read_csv(data_file_path)

    # clean data: not possible/not found replaced with NaN
    df = df.replace({'Not possible', 'Not found'}, np.NaN)
    df = df.drop('comment', axis=1)
    df = df.drop('Unnamed: 0', axis=1)
    df = df.apply(pd.to_numeric) # convert to numeric (integers and floats)
    df.describe()

    ### SUPPLEMENTARY FIGURE 1 ###
    #------------------------------
    plot_distance_transformations(output_folder_path)

    ### SUPPLEMENTARY FIGURE 2 ###
    #------------------------------
    # Compare runtimes and rank distribution between networks
    # define figure settings
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 9), )
    plt.subplots_adjust(hspace=0.4)
    ax11 = axs[0][0]
    ax12 = axs[0][1]
    ax21 = axs[1][0]
    ax22 = axs[1][1]

    # create boxplots (new dataframe for each plot)
    df_a = df[['CAR_exp_runtime', 'CAR_dflt_runtime', 'CAR_sqrt_runtime']]
    df_a.columns = ['exp', 'dflt', 'sqrt']
    df_a.boxplot(ax=ax11, grid=False)
    df_b = df[['CAR_034_runtime', 'KEGG_RPAIR_runtime', 'KEGG_nocofac_runtime']]
    df_b.columns = ['CAR > 0.34', '"main" RPAIRs', 'no cofactors']
    df_b.boxplot(ax=ax12, grid=False)
    df_c = df[['CAR_exp_rank', 'CAR_dflt_rank', 'CAR_sqrt_rank']]
    df_c.columns = ['exp', 'dflt', 'sqrt']
    df_c.boxplot(ax=ax21, grid=False)
    df_d = df[['CAR_034_rank', 'KEGG_RPAIR_rank', 'KEGG_nocofac_rank']]
    df_d.columns = ['CAR > 0.34', '"main" RPAIRs', 'no cofactors']
    df_d.boxplot(ax=ax22, grid=False)

    # Axis y labels
    ax11.set_ylabel('Algorithm runtime [s]')
    ax12.set_ylabel('Algorithm runtime [s]')
    ax21.set_ylabel('Rank of reference pathway')
    ax22.set_ylabel('Rank of reference pathway')

    # Axis x labels
    ax11.set_xlabel('Weighted networks')
    ax12.set_xlabel('Unweighted networks')
    ax21.set_xlabel('Weighted networks')
    ax22.set_xlabel('Unweighted networks')

    # add panel indication
    ax11.set_title('a')
    ax12.set_title('b')
    ax21.set_title('c')
    ax22.set_title('d')

    # save figure to pdf
    figure_2 = plt.figure(2)
    figure_2.savefig(output_folder_path+"figure_2.pdf", format="pdf")


    ### SUPPLEMENTARY FIGURE 3 ###
    #------------------------------
    # show distribution of length in reference pathway
    length_df = pd.DataFrame(df.Length.value_counts())
    length_df = length_df.reset_index()
    length_df.columns = ['length', 'count']
    length_df = length_df.sort_values(by='length')
    ax = length_df.plot.bar(x='length', y='count', ylabel='Number of pathways', xlabel = 'Pathway length (number of reaction steps)', legend=None)

    d = 1.5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
                  linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax.plot([0.66, 0.67], [0, 0], transform=ax.transAxes, **kwargs)
    ax.plot([0.77, 0.78], [0, 0], transform=ax.transAxes, **kwargs)
    ax.plot([0.88, 0.89], [0, 0], transform=ax.transAxes, **kwargs)

    # save figure to pdf
    figure_3 = plt.figure(3)
    figure_3.savefig(output_folder_path+"figure_3.pdf", format="pdf")


    ### SUPPLEMENTARY FIGURE 4 ###
    # ------------------------------
    # show overall performance in finding the reference pathways
    # define categories of ranks I want to plot
    index = ['top_hit', 'hit_2', 'hit_3', 'hit_4-5', 'hit_6-10', 'hit_11-100', 'None']
    # create new dataframe with ranks only
    df_rank = df[['CAR_exp_rank', 'CAR_dflt_rank', 'CAR_sqrt_rank', 'CAR_034_rank', 'KEGG_RPAIR_rank', 'KEGG_nocofac_rank']]
    # count occurences of categories in the new dataframe
    data_dict = count_occurrences(index, df_rank)
    # create new dataframe with results
    df_summary = pd.DataFrame(data_dict, index=index)
    # change columnn names for plotting
    df_summary.columns = ['exp', 'dflt', 'sqrt','CAR > 0.34', '"main" RPAIRs', 'no cofactors']

    # draw figure :
    df_summary = df_summary.T # tanspose table
    # fetch nice color palette from seaborn package
    color = sns.color_palette("icefire", 7)
    #rename columns
    df_summary.columns = ['Rank 1', 'Rank 2', 'Rank 3', 'Rank 4-6', 'Rank 6-10', 'Rank 11-100', 'Not found']
    # create barplot
    df_summary.plot.bar(stacked=True, rot=20,color=color,ylabel='Number of pathways')

    # save figure to pdf
    figure_4 = plt.figure(4)
    figure_4.savefig(output_folder_path+"figure_4.pdf", format="pdf")

def count_occurrences(index, df):
    D = {}
    header = list(df.columns.values)
    for column in header:
        data_list = []
        for i, cat in enumerate(index):
            if cat == 'top_hit':
                data_list.append(len(df[df[column] == 1]))
            elif cat == 'hit_2':
                data_list.append(len(df[df[column] == 2]))
            elif cat == 'hit_3':
                data_list.append(len(df[df[column] == 3]))
            elif cat == 'hit_4-5':
                data_list.append(count_range_of_values(df,column,4,5))
            elif cat == 'hit_6-10':
                data_list.append(count_range_of_values(df,column,6,10))
            elif cat == 'hit_11-100':
                data_list.append(count_range_of_values(df,column, 11, 100))
            elif cat == 'None':
                data_list.append(len(df[df[column].isna()]))

        D[column] = data_list
    return D

def count_range_of_values(df, column, range_start, range_end):
    count = 0
    for x in range(range_start,range_end+1):
        count += len(df[df[column] == x])
    return count

def plot_distance_transformations(output_folder_path):
    # create range of CARs
    cars = np.arange(0,1,0.01)
    exp = []
    dflt = []
    sqrt = []
    for car in cars:
        exp.append(calculateDistance(car,'exp'))
        dflt.append(calculateDistance(car,''))
        sqrt.append(calculateDistance(car,'sqrt'))
    df = pd.DataFrame([cars,exp,dflt,sqrt])
    df = df.T
    df.columns = ['CAR', 'Exponential', 'Default', 'Square root']
    df.set_index('CAR', inplace=True)

    fig, ax = plt.subplots()
    sns.lineplot(data=df, palette='rocket_r')
    plt.vlines(0.34, 0, 32, colors='black', label = 'Optimal cutoff (CAR = 0.34)')
    ax.set_ylim(0, 32)
    ax.set_ylabel('Distance')
    plt.legend(title = 'Transformation operator')

    figure_1 = plt.figure(1)
    figure_1.savefig(output_folder_path+"figure_1.pdf", format="pdf")

__main__()







