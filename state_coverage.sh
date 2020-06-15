#!/bin/bash
# This script, "state_coverage.sh" computes the number of bases covered by each state of the gene body in a search_state.py output file.
#
# Copyright (c) 2017 Michel TERESE
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

# For the gene body, return the number of bases covered by each state
# from a file containing the states and their lengths in the following format:
#
# 5 AT5G67330 26861266  26863897  + 3|3132|2  145|134,1649,449,396|496
#
# The file has no header
# Columns are:
# chromosome gene start end direction chrom_state lengths
#
# intergenic space before|gene body|intergenic space after

[[ $# -eq 1 && "$1" == '-h' ]] && printf "syntax:\n\t$0 search_state.py_output_file\nor\n\tsearch_state.py | $0\n" && exit 1

perl -wane '
  BEGIN {
    our @total_len= (0, 0, 0, 0, 0, 0, 0, 0, 0);

    sub handle_state {
      my $state= shift;

      $regexp_state= qr/$state/;
      if ( $state_string =~ /$regexp_state/ ) {
        ($len_string)= $F[6] =~ /^[^|]*\|([^|]+)/;
        @state_list= split //, $state_string;
        map { $_ = ($_ eq $state) ? "(\\d+)" : "\\d+" } @state_list;
        $regexp_str= "^" . join(",", @state_list);
        $regexp= qr/$regexp_str/;
        @len_list= $len_string =~ /$regexp/;
        map { $total_len[$state -1] += $_ } @len_list;
      }
    }

  }

  ($state_string)= $F[5] =~ /^[^|]*\|([^|]+)/;
  for my $state (1..9) {
    handle_state($state);
  }

  END {
    for my $state (1..9) {
      print "state $state : $total_len[$state -1]\n";
    }
  }
  ' $1
