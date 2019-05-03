use Bio::Root::IO;
use Bio::SeqIO;
use Bio::Seq;
use File::Copy;

$io = new Bio::Root::IO(-verbose => 1);
$td = $io->tempdir(CLEANUP => 0);
($tfh,$infilename) = $io->tempfile(-dir => $td);

$temp =  Bio::SeqIO->new('-fh'=>$tfh,
			 '-format' =>'Fasta');

my $contig_seq = Bio::Seq->new(-seq => 'tadayasdadasjdlasj', -id => 'contig');

$temp->write_seq($contig_seq);

$temp->close();
close($tfh);
undef $tfh;
copy($infilename, $td . '/xx');

print "$infilename\n";
