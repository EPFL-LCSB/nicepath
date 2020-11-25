############################################
# NICEpath - EPFL,LCSB - 2018-2020         #
############################################
# Author: Jasmin Hafner

# main
#######################################

import multiprocessing
import datetime
from util import *
from data import * #make a data object that holds everything
import sys


def main(_project_name):
    launch()
    start_time = time.time()

    # check for command line arguments
    if len(sys.argv) == 2:
        _project_name = sys.argv[1] # the first argument is considered as project name
    elif len(sys.argv) > 2:
        print('Input error: Too many arguments')
        return False

    # Initiate data object
    data = Data(_project_name)

    # Load parameters and input data
    data.loadParameters(default=True)
    data.loadParameters()
    data.loadData()

    # Print network statistics report
    if data.parameters.print_network_statistics:
        data.printNetworkStatistics()

    # Make sure the query is possible
    data.checkQueryCompounds()

    # Run pathway search for each target compound
    for target in data.parameters.target_compounds:
        print('\n################ TARGET: %s ################' % (target))
        # Run pathway search for each precursor compound
        data.parameters.target_compound = target

        if data.parameters.use_parallel:
            print('----------------\nSTATUS PARALLEL:\n----------------')
            output_stats = data.go()

        else:
            output_stats = []
            for precursor in data.parameters.precursor_compounds:
                output_stats.append(data.analyzePrecursor(precursor))

            if len(data.Pathways) != 0:
                outputfile_path = '%sPathways_%s.txt' % (data.outputfolder_path,target)
                outputfile = data.prepareOutput(outputfile_path)
                # write pathways to file
                data.writePathways(outputfile)
                outputfile.close()

        data.printStatistics(output_stats)
        data.printBridgITSystem(target)
        data.writeCompoundSequences()

    # Save dictionaries in pickle format for later reuse
    data.saveDictionaries()

    # collect statistics on run and pathways, save to text files
    runtime = time.time() - start_time
    data.printLog(runtime)

    print("RUN TIME:", str(datetime.timedelta(seconds=runtime)))
    return 0


def launch():
    print(" _   _  _____  _____  _____                _    _     ")
    print("| \ | ||_   _|/  __ \|  ___|              | |  | |    ")
    print("|  \| |  | |  | /  \/| |__   _ __    __ _ | |_ | |__  ")
    print("| . ` |  | |  | |    |  __| | '_ \  / _` || __|| '_ \ ")
    print("| |\  | _| |_ | \__/\| |___ | |_) || (_| || |_ | | | |")
    print("\_| \_/ \___/  \____/\____/ | .__/  \__,_| \__||_| |_|")
    print("                            | |                       ")
    print("                            |_|                     \n")


main('Test')
