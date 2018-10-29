# Imports only from standard and side packages
# No imports from ptes package
# General functions

import subprocess
import os
from collections import defaultdict
import logging
import logging.config

import yaml


def write_to_file(line, filename, folder='./'):
    if folder[-1] != '/':
        folder = folder + '/'
    output_file = open(folder + filename, 'a')
    output_file.write(line)
    output_file.close()


def writeln_to_file(line, filename, folder='./'):
    if folder[-1] != '/':
        folder = folder + '/'
    output_file = open(folder + filename, 'a')
    output_file.write(line+'\n')
    output_file.close()    


def run_command(command):   # to read from bash stdout
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,shell=True)
    return iter(p.stdout.readline, b'')


def shell_call(command):
    subprocess.call(command, shell = True)    


def init_file(filename, folder = './'):   # make new empty file
    if folder[-1] != '/':
        folder = folder + '/'
    with open(folder+filename, 'w') as output_file:
        pass 


def setup_logging(default_path,
                  default_level=logging.INFO,
                  env_key='LOG_CFG'):
    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
    if os.path.exists(path):
        with open(path, 'rt') as f:
            string = f.read()
            config = yaml.load(string)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)        


def print_time():
    now = datetime.datetime.now()
    print now.strftime("%Y-%m-%d %H:%M")        


def sign(x): return 1 if x >= 0 else -1


def sort_by_value(dct):
    print sorted(dct.items(), key=lambda x:x[1]) 


def read_cigar(row):   # reads CIGAR field from SAM file
    M = 0    
    Del = 0    
    Ins = 0
    S = 0
    N = 0
    num = ''
    for letter in row:
        if letter in '0123456789':
            num += letter
        elif letter == 'M':
            M += int(num)
            num = ''
        elif letter == 'D':
            Del += int(num)
            num = ''
        elif letter == 'I':
            Ins += int(num)
            num = ''
        elif letter == 'S':
            S += int(num)
            num = ''
        elif letter == 'N':
            N += int(num)
            num = ''    
    return {"M": M, "D": Del, "I": Ins, "S" : S, "N": N}        


def leave_unique(edge_dict, return_copies=True):
    '''
    Two nodes are connected by edge if by some rule they are similar;
    We'll have both node1-node2 and node2-node1 pairs in a dictionary;
    Takes a dictionary of all edges;
    Returns set of nodes that are all different (only one copy from set of all similar)
    Can also return a list of copies for each node in unique_set
    '''
    isSeen = defaultdict(lambda:False)
    isCopy = set()
    def step(node):
        isSeen[node] = True 
        for any_node in edge_dict[node]:   
            if not isSeen[any_node]:
                isCopy.add(any_node)   # the first node is not a copy
                return any_node
        return None

    for node in edge_dict.keys():
        new_node = node
        while new_node:
            new_node = step(new_node)
    
    unique_set = set()
    copies_dict = defaultdict(list)    
    for node in edge_dict.keys():
        if node not in isCopy:
            unique_set.add(node)
            for value in edge_dict[node]:
                copies_dict[node].append(value)
    if return_copies:
        return unique_set, copies_dict
    else:
        return unique_set   

def worker(kwargs):
    kwargs['function'](**kwargs)