use strict;

BEGIN { $ENV{CLUSTALDIR} = 'c:/BioEdit/apps' }
use Bio::Tools::Run::Alignment::Clustalw;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Registry;

# use primer3 wrapper from primer3plus webapi
BEGIN { push @INC,'.'; }
use primer3plusFunctions;

use Text::CSV;

# marks conserved parts in given exonseq by aligning it to contig (Bio::Seq) 
# using $aligner 
sub mark_conserved {
    my ($aligner, $contig, $exonseq) = @_;
    
    my @alignseq = ($exonseq, $contig);
    
    my $seq_ref = \@alignseq;
    my $aln = $aligner->align($seq_ref);

    # returns ? for non-aligned bases    
    my $cons = $aln->consensus_string(100);
    
    my $rawseq = $exonseq->seq();
    
    my $found_one = 0;
    
    # for all fragments consisting only of ACGT (no ?)
    while($cons =~ /([ACGT]{20,})/g) {
    	my $consfrag = $1;
    	# print "conserved: $consfrag \n";
    	
    	# mark conserved fragment for primer3
      # suppose unique match for each conserved fragment  
    	$rawseq =~ s/$consfrag/>$consfrag</;
    	
    	# flag a hit
    	$found_one = 1;
    }
    
    # invert the marks (primer3 requires <exclude> syntax)
    $rawseq = '<' . $rawseq . '>';
    $rawseq =~ s/^<>//;
    $rawseq =~ s/<>$//;
    
    return $found_one ? $rawseq : "";
}

# input is two Bio::EnsEMBL::Gene refs
# $v is verbose
# returns mumber of suggested primers
sub find_conserved_primers {
    my ($tg_gene, $gg_gene, $v) = @_;    
    
    #selection, sorting and printing of exons/introns of zebrafinch sequences
    my @ex = @{ $tg_gene->get_all_Exons };
    my @srt = sort {$a->start() <=> $b->start()} @ex;
    
    
    # select introns with length between 500 and 1500 bp 
    
    my $prevex = shift(@srt);
     
    
    #vytvori seznam do ktereho priradi start a end intronu s delkou 500 a 1500 bp
    my @introns = ();
    my @flankexons = ();
    
    foreach (@srt) {
    	my $intronlength = ($_->start() - $prevex->end());
    	if($v) {print "intron start: " . $prevex->end() . " end: " . $_->start() . " length: " . $intronlength;} 
    
    	# pokud intron odpovida delkou, pridej ho do seznamu 	
    	if (($intronlength < 1500) && ($intronlength > 500)){
    		push (@introns, ($prevex->end(), $_->start()));
    		push (@flankexons, {'prev' => $prevex, 'next' => $_});
    		if($v) {print " +";}
    			
    	}
    	if($v) {print "\n";}
    	$prevex = $_;
    }
    #    
    # ziskani kurecich exonu a slepeni do jednoho kontigu
    #
    my @gg_ex = @{ $gg_gene->get_all_Exons };
    my $gg_contig = "";
    
    foreach (@gg_ex) {
    	# zkraceny prikaz $gg_contig .= $_->seq->seq();
    	$gg_contig = $gg_contig . $_->seq->seq();
    }

    #
    # for each selected intron find conserved exon sequences and build primer3 input
    #
        
    #  Build a clustalw alignment factory
    my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM', 'quiet' => '1');
    my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
    
    # create object holding the pseudocontig created for easier aligment
    my $contig_seq = Bio::Seq->new(-seq => $gg_contig, -id => 'contig');        
    
    # primer counter 
    my $primers = 0;
    
    foreach(@flankexons) {
      # skip 'pseudointrons' with exons on different strands
      next if ($_->{'prev'}->strand() != $_->{'next'}->strand());
      
      # set the strand to the one of current exons
      my $intron = $tg_gene->slice->subseq($_->{'prev'}->end()+1, $_->{'next'}->start()-1, $_->{'prev'}->strand());
    
      my $exon1 = mark_conserved($factory, $contig_seq, $_->{'prev'}->seq);
      my $exon2 = mark_conserved($factory, $contig_seq, $_->{'next'}->seq);
      
      next if ($exon1 eq "" || $exon2 eq "");
      my $primer3 = $exon1 . '[' . $intron . ']' . $exon2;
      
      print "primer3>\n" . $primer3 . "\n" if ($v); 
    
      # run primer3
      my (%paramsHash, %resHash);
      
      $paramsHash{"SCRIPT_DETECTION_PICK_LEFT"} = "1";
      $paramsHash{"SCRIPT_DETECTION_PICK_RIGHT"} = "1";
      $paramsHash{"SEQUENCE"} = $primer3;
      $paramsHash{"PRIMER_PRODUCT_OPT_SIZE"} = "800";
      $paramsHash{"PRIMER_PRODUCT_SIZE_RANGE"} = "500-1500";
      $paramsHash{"SCRIPT_DETECTION_USE_PRODUCT_SIZE"} = "0";
      $paramsHash{"PRIMER_NUM_RETURN"} = "1";
    
      # checkParameters parses the sequence 'markup' created just above
      checkParameters(%paramsHash);
    
      # call the external program (primer3.exe)
      $paramsHash{"SCRIPT_TASK"} = "Detection";
      findAllPrimers(\%paramsHash, \%resHash);
    
      next if (defined $resHash{"PRIMER_ERROR"} and $resHash{"PRIMER_ERROR"} ne "");
      
      # no primers found
      next if (!defined $resHash{"PRIMER_PRODUCT_SIZE_0"});
            
      # output the results - using only first one..
      # foreach my $i (0..4) {
    
      # use only first one
      my $i = 0;
    
    	# set output field separator to tab and record sep to newline
    	# storing previous values
    	my ($ofs, $ors) = ($, , $\);
    	($, , $\) = ("\t", "\n");
    	
    	# left primer
    	my ($start, $len) = split(/,/, $resHash{"PRIMER_LEFT\_$i"});
    	print $tg_gene->stable_id, $i, "left", $start, $len, $resHash{"PRIMER_LEFT\_$i\_TM"}, $resHash{"PRIMER_LEFT\_$i\_SEQUENCE"}, $resHash{"PRIMER_PRODUCT_SIZE\_$i"};
    
    	# right primer
    	($start, $len) = split(/,/, $resHash{"PRIMER_RIGHT_$i"});
    	print $tg_gene->stable_id, $i, "right", $start, $len, $resHash{"PRIMER_RIGHT\_$i\_TM"}, $resHash{"PRIMER_RIGHT\_$i\_SEQUENCE"}, $resHash{"PRIMER_PRODUCT_SIZE\_$i"};

      # restore old separators    	
    	($, , $\) = ($ofs, $ors);
      #}
      
      # increment primer counter
      $primers++;
    }
    
    return $primers;
}

#
# check arguments (single filename required)
# 
if($#ARGV != 0) {
  print STDERR "use: $0 input_csv_file\n";
  exit(-1);
}

#
# initialize EnsEMBL connection
#
my $registry = 'Bio::EnsEMBL::Registry';

print STDERR "connecting to registry\n";
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

# create csv parser
my $csv = Text::CSV->new();

# output csv header
print "tgID,primerID,side,start,len,tm,seq,prod_size\n";


open (CSV, "<", $ARGV[0]) or die $!;

# read a single line to skip the headers
<CSV>;

while (<CSV>) {
  if ($csv->parse($_)) {
    my @columns = $csv->fields();
    my ($tg_gene_id, $gg_gene_id) = @columns;
  
    print STDERR "$tg_gene_id: ";
    
    #download sequence of zebrafinch
    my $tg_adaptor = $registry->get_adaptor('zebrafinch', 'Core', 'Gene' );
    my $tg_gene = $tg_adaptor->fetch_by_stable_id($tg_gene_id);
    
    #download sequence of chicken
    my $gg_adaptor = $registry->get_adaptor('chicken', 'Core', 'Gene' );
    my $gg_gene = $gg_adaptor->fetch_by_stable_id($gg_gene_id);
  
    print STDERR find_conserved_primers($tg_gene, $gg_gene);
    print STDERR "\n";
  }
  else {
      my $err = $csv->error_input;
      print "Failed to parse line: $err\n";
  }

}
close CSV;
