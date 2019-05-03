#!/cygdrive/c/Perl/bin/perl.exe

use warnings;
use strict;

use Bio::Tools::Run::Primer3;
use Bio::SeqIO;

my $file_in='piwil2.fa';
my $file_out='temp.out';
      
      # read a file and return a sequence object 
      my $seqio=Bio::SeqIO->new(-file=>$file_in);# read a file and return a sequence object
     
      # What are the returns from $seqio:
      #print " These are the returns from \$seqio\:\n";
      #while (my ($key, $value)=each %$seqio){print "\t$key=>\t$value\n"}
      
      #These are the results from the screen
      
      # These are the returns from $seqio:
      #  _seqio_seqfactory=>     Bio::Seq::SeqFastaSpeedFactory=HASH(0x105a3110)
      #  _file=> piwil2.fa
      #  _root_cleanup_methods=> ARRAY(0x105a3080)
      #  _flush_on_write=>       1
      #  _filehandle=>   GLOB(0x105a3044)
      #  _root_verbose=> 0
      #  _object_builder=>       Bio::Seq::SeqBuilder=HASH(0x105a2e10)
      #  _seqio_locfactory=>     Bio::Factory::FTLocationFactory=HASH(0x105a2d14)
      
        
      
     # use a loop to deal with multiple sequnces  
     
       # What are the returns from $seq: 
       #print " These are the returns from \$seq\:\n";                     

while (my $seq=$seqio->next_seq){ 
      
            #What are the returns from $seq
      		#while (my ($key, $value)=each %$seq){print "\t$key=>\t$value\n"}
      		
     #These are the results from the screen 		
     #These are the returns from $seq:
     #primary_id=>    gi|10946609|ref|NM_021308.1|
     #primary_seq=>   Bio::PrimarySeq=HASH(0x105b753c)
     #primary_id=>    gi|6671508|ref|NM_007393.1|
     #primary_seq=>   Bio::PrimarySeq=HASH(0x105b7bcc)
     # I get duplicate results due to 2 sequences in this file
     
     # base on the above results I figure out two subrountines/keys within $seq object:
     # 1)method "primary_id" returns a string, which is what I want
     # 2)method "primary_seq" returns a hash reference and I need to put"%" to deference it
    
     # access key/subroutine primary_id to get value primary_id
     my $primary_id=$seq->primary_id; 
     
     # parse the header line 
     my @primary_id=split( /\|/,$primary_id);
     print "@primary_id", "\n\n";
     #print "\n","\primary_id\t$primary_id\n\n";
     
     # access key/subroutine primary_seqto get value primary_seq
     # since it is a reference to a hash I need to derefernce it by adding % 
      my $primary_seq=$seq->seq;
     
      #print ref $primary_seq,"\n\n";       
      print "\n",$primary_seq, "\n\n" ;
      
	      
	  # read fasta file into Primer3
      my $primer3=Bio::Tools::Run::Primer3->new(
                                            -seq=>$seq,
                                            # I disable line -outfile 
                                            #because it only returns the information for the last sequences
                                            #-outfile=>$file_out, 
                                            -path=>"c:/Perl/local/primer3_1.0.0/src/primer3.exe"  # the path where primer3 executable file is located
                                            );
  
    #set your own parameters                       
    $primer3->add_targets(		
		'PRIMER_OPT_GC_PERCENT'=>  '50',		
		'PRIMER_OPT_SIZE'=>  '24',		
		'PRIMER_OPT_TM'=>  ' 60 ',		
		'PRIMER_PRODUCT_OPT_SIZE'=>  '500',		
		'PRIMER_PRODUCT_SIZE_RANGE'=>  '490-510',		
		'PRIMER_WT_GC_PERCENT_GT'=>  ' 50 ',		
		'PRIMER_WT_GC_PERCENT_LT'=>  ' 50 ',		
	);  
 
  	# run the primer3.exe program and return the results
   my $result=$primer3->run; 
   
   #What results are expected
   #print ref $result;
   
   # I get the results from screen
   # Bio::Tools::Primer
   # Based on this information I know this line asks Bio::Tools::Run::Primer to invoke Bio::Tools::Run::Primer3
   # And I think it explains why I can parse the results using the syntax from  Bio::Tools::Primer3   
   
  
  
   my $result1=$result->primer_results(1);#return a hash reference from  Bio::Tools::Run::Primer
   
	   
	  # This loop tells what the keys/values are expected to print out
	  # by default it retruns all the values for all primers it designs
	  # Since I just need part of the information  for the first pair of primer I provide 
	  #the individual key and ask perl to return its corresponding value
	  
	  print " These are the keys available","\n\n";
	         while (my($key,$value)=each %$result1) {print "$key\n"}
	   
	   my $key_PRIMER_LEFT_SEQUENCE='PRIMER_LEFT_SEQUENCE';   	   
	   print "\n$key_PRIMER_LEFT_SEQUENCE\t${$result1}{$key_PRIMER_LEFT_SEQUENCE}\n";
	   
	   my $key_PRIMER_RIGHT_SEQUENCE='PRIMER_RIGHT_SEQUENCE'; 
	   print "$key_PRIMER_RIGHT_SEQUENCE\t${$result1}{$key_PRIMER_RIGHT_SEQUENCE}\n";
	   
	   
       my $key_PRIMER_LEFT_TM='PRIMER_LEFT_TM';
       print "$key_PRIMER_LEFT_TM\t${$result1}{$key_PRIMER_LEFT_TM}\n";
       
       my $key_PRIMER_RIGHT_TM='PRIMER_RIGHT_TM';
       print "$key_PRIMER_RIGHT_TM\t${$result1}{$key_PRIMER_RIGHT_TM}\n"; 
       
       my $key_PRIMER_PRODUCT_SIZE='PRIMER_PRODUCT_SIZE';
        print "$key_PRIMER_PRODUCT_SIZE\t${$result1}{$key_PRIMER_PRODUCT_SIZE}\n"; 
       
       print "\n";   
  
}
   exit;
                        	
    
