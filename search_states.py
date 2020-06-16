#!/usr/bin/env python3
#
# Copyright (c) 2017 Michel TERESE (seeterlab)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


VERSION = '1.2'
CHROMOSOME_COUNT = 5
NB_BASES = 1000

import sys
from bisect import bisect

# Check arguments
if len(sys.argv) -1 != 2 :
    print( "syntax: ", sys.argv[0], 
"""chrom_states_file data_file

INPUT:

chrom_states_file describe chromatin states along the chromosome. It contains 4 columns as below:

Chrom	From	To	State
1	1	1200	8
1	1201	3150	4
1	3151	4650	2
1	4651	6450	6
...

data_file must contain tab separated columns : chromosome number, gene ID, start, stop, direction. It must be sorted by chromosome number and start with chromosome 1
example:

1       AT1G01010       3631    5899    +
1       AT1G01020       6788    9130    -
1       AT1G03987       11101   11372   +
1       AT1G01030       11649   13714   -
1       AT1G01040       23121   31227   +
1       AT1G03993       23312   24099   -
...

OUTPUT:

The ouput file contains the 5 columns from dztz_file plus 2 new columns :
- States list.
This first supplementary column contain 3 fields separated with a "|". The field contains the states chain for intergenic space before the gene | gene body | intergenic space after the gene. The length of intergenic spaces are restricted to NB_BASES.
- States length
The lengths of each states are coma separated

1       AT1G01010       3631    5899    +        9|26|7       10|50,10|123
1       AT1G01020       6788    9130    -        31|42137|78  6,25|12,456,10,23,58|45,47
...

"""
)
    sys.exit(1)


#================================================================================
def get_state_len(state_ndx, chrom, range_start, range_stop):
    "return the lenght of the state indexed state_ndx depending on the start and stop of the considered interval"

    state_start = max(a_start[chrom][state_ndx], range_start)
    state_stop  = min(a_stop[chrom][state_ndx], range_stop)
    return str(state_stop - state_start + 1)


#================================================================================
def generate_state_list( chrom, start, stop, direction ):
    "generate a list of states with their length for the considered interval"

    # Search the state for the starting position of the current gene
    i = bisect(a_start[chrom], start) -1
    # If there is no state in the range "start-stop", a tuple of empty list is returned
    if a_stop[chrom][i] < start:
        return ( '', '' )

    state_list = a_state[chrom][i]
    len_list = []
    len_list.append( get_state_len(i, chrom, start, stop) )

    # Search for next states in the interval
    for j in range(i +1, len(a_start[chrom])) :
        if a_start[chrom][j] > stop :
            break
        # sate is a string!
        state_list += a_state[chrom][j]
        len_list.append( get_state_len(j, chrom, start, stop) )

    # If gene direction is minus, the list of states and lengths are reversed
    if direction == "-":
        state_list = state_list[::-1]
        len_list = len_list[::-1]

    return (state_list, ','.join(len_list))

#================================================================================
def print_line(line, before_state_list, state_list, after_state_list, direction):
    "print resulting line"

    # If prev_direction is '-' before and after state list are switched
    if direction == '-' :
        before_state_list, after_state_list = after_state_list, before_state_list

    print( line + "\t" + before_state_list[0] + '|' + state_list[0] + '|' + after_state_list[0] +
                  "\t" + before_state_list[1] + '|' + state_list[1] + '|' + after_state_list[1] )


#================================================================================
# main
#================================================================================
chrom_states_file = sys.argv[1]
data_file = sys.argv[2]

# The 3 tables below are indexed by chromosome number (from 1 to n)
# Each element of the table is a liste of position or states for the considered chromosome
# Index 0 is initialized to None as it is not used
a_start = []
a_stop = []
a_state = []

# Initialize the 3 tables
for i in range(0, CHROMOSOME_COUNT +1) :
    a_start.append(i)
    a_start[i] = [] if i > 0 else None
    a_stop.append(i)
    a_stop[i] = [] if i > 0 else None
    a_state.append(i)
    a_state[i] = [] if i > 0 else None

# read the state file and initialize the tables a_start, a_stop et a_state
with open( chrom_states_file ) as f:
    line_nb = 0
    for line in f :
        line_nb += 1
        # saute la 1ère ligne (en-tête)
        if line_nb == 1 :
            continue
        (chrom, start, stop, state) = line.split()
        chrom = int(chrom)
        a_start[chrom].append( int(start) )
        a_stop[chrom].append( int(stop) )
        a_state[chrom].append( state )

# Previous stop position
prev_stop = 0
# Previous direction
prev_direction = ''
# Previous line
prev_line = ''
# Previous gene state list
state_list = ''
before_state_list = ''

# Previous chomosom
prev_chrom = 1

# read the data file
with open( data_file ) as f:
    for line in f :
        fields = line.split()
        chrom     = int(fields[0])
        start     = int(fields[2])
        stop      = int(fields[3])
        direction = fields[4]

        # If it is a new chromosome
        if chrom != prev_chrom :
            # Compute the after_state_list at the end of the previous chromosome
            inter_gene_start = prev_stop + 1
            inter_gene_stop  = inter_gene_start + NB_BASES
            after_state_list = generate_state_list( prev_chrom, inter_gene_start, inter_gene_stop, prev_direction )
            # Print the previous line
            print_line(prev_line, before_state_list, state_list, after_state_list, prev_direction)

            # Reinitiate prev_stop for the next chromosome
            prev_chrom = chrom
            prev_stop = 0
            prev_line = ''

        # If at least one line has been read, compute the after_state_list and print the previous line
        if prev_line != '' :
            # If there is an intergenic space (no overlapping)
            inter_gen_length = start - prev_stop -1
            if inter_gen_length > 0 :
                inter_gene_start = prev_stop + 1
                inter_gene_stop  = inter_gene_start + min(NB_BASES, inter_gen_length)
                # Compute after_state_list
                after_state_list = generate_state_list( chrom, inter_gene_start, inter_gene_stop, prev_direction )
            else :
                after_state_list = ('' , '')

            # Print the previous line
            print_line(prev_line, before_state_list, state_list, after_state_list, prev_direction)

        # If there is an intergenic space (no overlapping)
        inter_gen_length = start - prev_stop -1
        if inter_gen_length > 0 :
            # Compute before_state_list
            inter_gene_stop   = start - 1
            inter_gene_start = inter_gene_stop - min(NB_BASES, inter_gen_length)
            if inter_gene_start <= 0:
                inter_gene_start = 1
            before_state_list = generate_state_list( chrom, inter_gene_start, inter_gene_stop, direction )
        else :
            before_state_list = ('', '')

        # Compute the state list of the gene
        state_list = generate_state_list( chrom, start, stop, direction )
        prev_line = line.rstrip()

        prev_stop = stop
        prev_direction = direction

# Print the last line
# Compute after_state_list at the end of the previous chromosome
inter_gene_start = prev_stop + 1
inter_gene_stop  = inter_gene_start + NB_BASES
after_state_list = generate_state_list( prev_chrom, inter_gene_start, inter_gene_stop, prev_direction )
# Print the previous line
print_line(prev_line, before_state_list, state_list, after_state_list, prev_direction)

