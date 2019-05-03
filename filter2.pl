  use strict;
  use warnings;
  use Text::CSV;

  my $file = 'mart_export.txt';
  my $csv = Text::CSV->new();

  my %genes_printed = ();

  open (CSV, "<", $file) or die $!;

  # print csv header
  print "tgID,gg_orthID,dN,dS,tgIdent,ggIdent,geneName,start\n";

    while (<CSV>) {
        if ($csv->parse($_)) {
            my @columns = $csv->fields();
            if ($columns[5] eq "") {
	       $columns[6] = $columns[7] = 0;
	    }   
	    if (($columns[14] eq 'ortholog_one2one') &&
	        ($columns[13] eq 'Z') &&
	        ($columns[11] < 0.5) &&
		($columns[6] < 70) &&
		($columns[7] < 70) &&
		!($genes_printed{$columns[0]})) {
		
			print "$columns[0],$columns[8],$columns[12],$columns[11],$columns[9],$columns[10],$columns[4],$columns[2]\n";
			$genes_printed{$columns[0]} = 1;
	    }
        } else {
            my $err = $csv->error_input;
            print "Failed to parse line: $err";
        }
    }
    

close CSV;
