#!/bin/bash
# Affiche le nombre de bases du corps du gène couvert par chaque état
# à partir d'un fichier qui contient les états et les longueurs correspondante au format suivant:
#
# 5 AT5G67330 26861266  26863897  + 3|3132|2  145|134,1649,449,396|496
#
# Le fichier ne doit pas contenir d'entête
# Les colonnes correspondent à:
# chromosome gene start stop sens état_chromatinien longueurs
#
# espace intergénique avant|corps du gène|espace intergénique après

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
