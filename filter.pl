  use strict;
  use warnings;
  use Text::CSV;

  my $file = 'mart_export.txt';
  my $csv = Text::CSV->new();

  my %genes_printed = ();


  # print csv header
  print "tgID,gg_orthID,dN,dS,tgIdent,ggIdent,geneName,start\n";

  open (CSV, "<", $file) or die $!;
  while (<CSV>) {
        if ($csv->parse($_)) {
            my @columns = $csv->fields();
            if ($columns[5] eq "") {
	       $columns[6] = $columns[7] = 0;
	    }   
	    if (($columns[14] eq 'ortholog_one2one') &&
	        ($columns[13] eq 'Z') &&
	        ($columns[9] > 90) &&
		($columns[10] > 90) &&
		($columns[6] < 70) &&
		($columns[7] < 70) &&
		!($genes_printed{$columns[1]})) {
		
			print "$columns[1],$columns[8],$columns[12],$columns[11],$columns[9],$columns[10],$columns[4],$columns[2]\n";
			$genes_printed{$columns[1]} = 1;
	    }
        } else {
            my $err = $csv->error_input;
            print "Failed to parse line: $err";
        }
   }

   close CSV;
  
   my $pocet_genu = keys %genes_printed;
   print STDERR "number of selected genes: $pocet_genu\n";


