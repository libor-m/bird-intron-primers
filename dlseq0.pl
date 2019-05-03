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

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $tg_gene_id = 'ENSTGUG00000000006';

#download sequence of zebrafinch
my $tg_adaptor = $registry->get_adaptor('zebrafinch', 'Core', 'Gene' );
my $tg_gene = $tg_adaptor->fetch_by_stable_id($tg_gene_id);

#download sequence of chicken
my $gg_adaptor = $registry->get_adaptor('chicken', 'Core', 'Gene' );
my $gg_gene = $gg_adaptor->fetch_by_stable_id('ENSGALG00000001843');


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
	print "intron start: " . $prevex->end() . " end: " . $_->start() . " length: " . $intronlength; 

	# pokud intron odpovida delkou, pridej ho do seznamu 	
	if (($intronlength < 1500) && ($intronlength > 500)){
		push (@introns, ($prevex->end(), $_->start()));
		push (@flankexons, {'prev' => $prevex, 'next' => $_});
		print " +";
			
	}
	print "\n";
	$prevex = $_;
}



# ziskani kurecich exonu a slepeni do jednoho kontigu
my @gg_ex = @{ $gg_gene->get_all_Exons };
my $gg_contig = "";

foreach (@gg_ex) {
	# zkraceny prikaz $gg_contig .= $_->seq->seq();
	$gg_contig = $gg_contig . $_->seq->seq();
}

#
# use clustalw to find conserved reginos between gg and tg
#

#  Build a clustalw alignment factory
my @params = ('ktuple' => 2, 'matrix' => 'BLOSUM');
my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);

# create object holding the pseudocontig created for easier aligment
my $contig_seq = Bio::Seq->new(-seq => $gg_contig, -id => 'contig');

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
    	print "conserved: $consfrag \n";
    	
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

#
# for each selected intron find conserved exon sequences and build primer3 input
#

foreach(@flankexons) {
  # skip 'pseudointrons' with exons on different strands
  next if ($_->{'prev'}->strand() != $_->{'next'}->strand());
  
  # set the strand to the one of current exons
  my $intron = $tg_gene->slice->subseq($_->{'prev'}->end()+1, $_->{'next'}->start()-1, $_->{'prev'}->strand());

  my $exon1 = mark_conserved($factory, $contig_seq, $_->{'prev'}->seq);
  my $exon2 = mark_conserved($factory, $contig_seq, $_->{'next'}->seq);
  
  next if ($exon1 eq "" || $exon2 eq "");
  my $primer3 = $exon1 . '[' . $intron . ']' . $exon2;
  print "primer3>\n" . $primer3 . "\n"; 

  # run primer3
  my (%paramsHash, %resHash);
  
  $paramsHash{"SCRIPT_DETECTION_PICK_LEFT"} = "1";
  $paramsHash{"SCRIPT_DETECTION_PICK_RIGHT"} = "1";
  $paramsHash{"SEQUENCE"} = $primer3;
  $paramsHash{"PRIMER_PRODUCT_OPT_SIZE"} = "800";
  $paramsHash{"PRIMER_PRODUCT_SIZE_RANGE"} = "500-1500";
  $paramsHash{"SCRIPT_DETECTION_USE_PRODUCT_SIZE"} = "0";

  # checkParameters parses the sequence 'markup' created just above
  checkParameters(%paramsHash);

  # call the external program
  $paramsHash{"SCRIPT_TASK"} = "Detection";
  findAllPrimers(\%paramsHash, \%resHash);

  # output the results
  # foreach my $i (0..4) {

  # use only first one
  my $i = 0;

	# set output field separator to tab and record sep to newline
	$, = "\t";
	$\ = "\n";
	
	# left primer
	my ($start, $len) = split(/,/, $resHash{"PRIMER_LEFT\_$i"});
	print $tg_gene_id, $i, "left", $start, $len, $resHash{"PRIMER_LEFT\_$i\_TM"}, $resHash{"PRIMER_LEFT\_$i\_SEQUENCE"}, $resHash{"PRIMER_PRODUCT_SIZE\_$i"}, "\n";

	# right primer
	($start, $len) = split(/,/, $resHash{"PRIMER_RIGHT_$i"});
	print $tg_gene_id, $i, "right", $start, $len, $resHash{"PRIMER_RIGHT\_$i\_TM"}, $resHash{"PRIMER_RIGHT\_$i\_SEQUENCE"}, $resHash{"PRIMER_PRODUCT_SIZE\_$i"}, "\n";
  #}
}
