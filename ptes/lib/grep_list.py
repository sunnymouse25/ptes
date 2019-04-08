# makes script with N large strings - commands for grepping a list of patterns in file in background
# python -i LIST -n N[10] -d DATABASE

import argparse

from ptes.lib.general import shell_call

# Functions


def main():
    # Arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        help="Name of input list to grep")
    parser.add_argument("-n", "--n_parts", type=int,
                        default=10,
                        help="How many times to grep")
    parser.add_argument("-db", "--database", type=str,
                        help="Where to grep")
    parser.add_argument("-o", "--output", type=str,
                        default='to_grep.sh',
                        help="Name of script with grep commands")
    parser.add_argument("-p", "--prefix", type=str,
                        default='',
                        help="Common prefix, command after cat $input")
    parser.add_argument("-s", "--suffix", type=str,
                        default='',
                        help="Common suffix, command before >> $output &")
    args = parser.parse_args()
    # Main

    with open(args.input,'r') as inp_list:
        filelist = inp_list.readlines()
        n_patterns = len(filelist) // args.n_parts   # how many patterns in one command
        common_prefix = "cat %s | " % args.database + args.prefix
        common_suffix = args.suffix + " >> %s.grep &" % args.input
        out_list = []
        if n_patterns > 1:
            for i in range(args.n_parts-1):
                out_list.append(common_prefix +
                                "grep '" +
                                "\|".join(
                                    [x.strip('\n') for x in filelist[(i*n_patterns): (i*n_patterns + n_patterns)]]) +
                                "'" +
                                common_suffix)
            out_list.append(common_prefix +
                            "grep '" +
                            "\|".join([x.strip('\n') for x in filelist[(i*n_patterns + n_patterns):]]) +
                            "'" +
                            common_suffix)
        else:
            out_list.append(common_prefix +
                            "grep '" +
                            "\|".join([x.strip('\n') for x in filelist]) +
                            "'" +
                            common_suffix)
    with open(args.output, 'w') as out_file:
        out_file.write('\n'.join(out_list))

    shell_call('chmod +x %s' % args.output)


if __name__ == "__main__":
    main()