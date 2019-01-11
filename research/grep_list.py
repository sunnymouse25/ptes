# makes script with N large strings - commands for grepping a list of patterns in file in background
# python -i LIST -n N[10] -d DATABASE

import os
import argparse

from ptes.lib.general import init_file, writeln_to_file, shell_call

# arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str,                    
                    help="Name of input list to grep")
parser.add_argument("-n","--n_parts", type=int,
                    default = 10,                    
                    help="How many times to grep")                    
parser.add_argument("-db","--database", type=str,                    
                    help="Where to grep")
parser.add_argument("-o","--output", type=str,
                    default='to_grep.sh',
                    help="Name of script with grep commands")
args = parser.parse_args()


# Functions


# Main
inp_list_file = args.input
n_parts = args.n_parts
database = args.database
script_name = args.output

path_to_file = os.path.dirname(os.path.realpath(inp_list_file.strip()))
init_file(script_name, folder = path_to_file)

with open(inp_list_file,'r') as inp_list:
    filelist = inp_list.readlines()
    n_patterns = len(filelist) // n_parts   # how many patterns in one command
    common_prefix = "cat %s | " % database
    common_suffix = " >> %s.grep &" % inp_list_file
    if n_patterns > 1:
        for i in range(n_parts-1):
            writeln_to_file(common_prefix +
                  "grep '" + "\|".join([x.strip('\n') for x in filelist[(i*n_patterns): (i*n_patterns + n_patterns)]]) +
                  "'" + common_suffix,
                            script_name,
                            folder=path_to_file)
        writeln_to_file(common_prefix + "grep '" +
              "\|".join([x.strip('\n') for x in filelist[(i*n_patterns + n_patterns):]]) + "'" + common_suffix,
                        script_name,
                        folder=path_to_file)
    else:
        writeln_to_file(common_prefix + "grep '" + "\|".join([x.strip('\n') for x in filelist]) + "'" + common_suffix,
                        script_name,
                        folder=path_to_file)

shell_call('chmod +x %s/%s' % (path_to_file, script_name))