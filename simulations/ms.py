#! /usr/bin/python

"""
ms.py 

This module is for reading the simulation results from ms.

Sam Vohr (svohr@ucsc.edu)
October 21, 2012
"""

import collections

class MSResult:
    """ Wrapper class for a stream containing the results of a set of 
        ms simulations. """
    cmdraw = ''
    nsamples = 0
    nreps = 0
    args = collections.defaultdict(list)
    rseed = ''
    ms_input = None
    curr_sim = None

    def _read_args(self):
        items = self.cmdraw.split('-')
        required = items[0].split()
        self.nsamples = int(required[1])
        self.nreps = int(required[2])
        for arg in items[1:]:
            flag, values = arg[0:arg.find(' ')], arg[arg.find(' ')+1:].rstrip()
            self.args[flag].append(values)


    def __init__(self, ms_input):
        """ creates a new MSResult based on the input stream. """
        self.ms_input = ms_input
        # First line is the command with arguments for this simulation.
        self.cmdraw = self.ms_input.readline().rstrip() 
        self._read_args()
        self.rseed = self.ms_input.readline().rstrip()
        # Advance the file to the next simulation.
        line = self.ms_input.readline()
        while not line.startswith('//'):
            line = self.ms_input.readline()

    def next_sim(self):
        """ Reads the next simulation from the file. """
        lines = []
        line = self.ms_input.readline()
        if line == '':
            # end of file
            return None
        else:
            while not (line.startswith('//') or line == ''):
                if line != '\n':
                    # if line is not blank, add it to the list for this sim.
                    lines.append(line.rstrip())
                line = self.ms_input.readline()
            # hit the next simulation or end of file
            ms_current_sim = MSSim(self,lines)
            return ms_current_sim
    
    def pop_sizes( self ):
        """ Reads the population sizes from the ms arguments """
        if 'I' in self.args:
            pop_str = self.args['I'][0]
            items = pop_str.split()
            npops = int(items[0])
            pops = [ int(p) for p in items[1:] ]
            assert npops == len(pops)
            assert self.nsamples == sum(pops)
        else:
            pops = list( [ self.nsamples ] )
        return pops




class MSSim:
    """ Stores the results of a single simulation """
    nsites = 0
    pos = []
    genotypes = []
    ms_parent = None
    def __init__(self, parent, ms_lines):
        self.ms_parent = parent
        self.genotypes = list()
        site_line = ms_lines[0]
        assert site_line.startswith('segsites:')
        self.nsites = int( site_line[site_line.find(' ')+1:] )
        pos_line = ms_lines[1]
        assert pos_line.startswith('positions:')
        self.pos = pos_line.split()[1:]
        self.genotypes = ms_lines[2:]


    def pop_genotypes(self):
        pop_size = self.ms_parent.pop_sizes()
        pop_geno = list()
        cur = 0
        for p in pop_size:
            pop_geno.append( self.genotypes[cur:cur+p] )
            cur += p
        return pop_geno
