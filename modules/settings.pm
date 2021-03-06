#!/usr/bin/perl -w

#  Copyright (c) 2006, 2007
#  by Andreas Untergasser and Harm Nijveen
#  All rights reserved.
# 
#  This file is part of Primer3Plus. Primer3Plus is a webinterface to primer3.
# 
#  The Primer3Plus is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  Primer3Plus is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Primer3Plus (file gpl.txt in the source distribution);
#  if not, write to the Free Software Foundation, Inc., 51 Franklin St,
#  Fifth Floor, Boston, MA  02110-1301  USA

package settings;
use strict;
use CGI;
use Carp;
use CGI::Carp qw(fatalsToBrowser);
use Exporter;
our (@ISA, @EXPORT, @EXPORT_OK, $VERSION);
@ISA = qw(Exporter);
@EXPORT = qw(&getDefaultSettings &getMachineSetting &setDoNotPick &getDoNotPick &getMisLibrary
             &getMachineSettings &getScriptTask &getLibraryList &getSaveSequenceParameters
             &getSaveSettingsParameters &getPrimer3CompleteParameters &getZeroReplacements
             &getPrimerCheckParameters &getServerParameterFiles &getServerParameterFilesList 
             &setMessage &getMessages);

$VERSION = "1.00";

# Here it stores all the messages from within the functions
my @messages;
my $doNotPick = 0;

##############################################################################
# ---------------- Installer Modifiable Variables -------------------------- #
# You may wish to modify the following variables to suit your installation.  #
##############################################################################

my %machineSettings = (
  # Who the end user will complain to:
  "MAINTAINER" =>"user&#host.com",

  # The location of the primer3_core executable.
  "PRIMER_BIN" =>  "primer3.exe",     # for Windows
# "PRIMER_BIN" =>  "./primer3_core",   # for Linux

  # Parameters which are handed in with the programm call.
  "PRIMER_RUNTIME" =>  " -strict_tags",     # for Windows
# "PRIMER_RUNTIME" =>  " -strict_tags",     # for Linux

  # The URL were to find the template HTML
  "URL_HTML_TEMPLATE" => "HtmlTemplate.html",

  # The URL for the form action (which will normally be the filename of the script)
  "URL_FORM_ACTION" => "primer3plus.cgi",

  # The URL for the result action (were primer3plus sends his found primer to) 
  "URL_PRIMER_MANAGER" => "primer3manager.cgi",

  # The URL for help regarding this screen (which will normally
  # be in the same directory as the this script)
  "URL_HELP" => "primer3plusHelp.cgi",

  # The path were primer3Manager finds its parameter files for server stored settings.
  "USER_PARAMETERS_FILES_PATH" =>  "./parameter_files/",

  # The URL to download the source code 
  "URL_SOURCE" => "http://sourceforge.net/projects/primer3/",

  # The URL for information about Primer3plus
  "URL_ABOUT" => "primer3plusAbout.cgi",

  # The URL were to find the FAQ
  "URL_FAQ" => "primer3plusHelp.cgi#FAQ",

  # The URL to use BLAST
  "URL_BLAST" => qq{<a href="http://0-www.ncbi.nlm.nih.gov.library.vu.edu.au:80/blast/Blast.cgi?ALIGNMENTS=50&amp;ALIGNMENT_VIEW=Pairwise&amp;AUTO_FORMAT=Semiauto&amp;CLIENT=web&amp;DATABASE=nr&amp;DESCRIPTIONS=100&amp;ENTREZ_QUERY=All+organisms&amp;EXPECT=1000&amp;FORMAT_BLOCK_ON_RESPAGE=None&amp;FORMAT_ENTREZ_QUERY=All+organisms&amp;FORMAT_OBJECT=Alignment&amp;FORMAT_TYPE=HTML&amp;FULL_DBNAME=nr&amp;GET_SEQUENCE=on&amp;HITLIST_SIZE=100&amp;JOB_TITLE=Nucleotide+sequence+(26+letters)&amp;LAYOUT=TwoWindows&amp;MASK_CHAR=2&amp;MASK_COLOR=1&amp;MYNCBI_USER=4308031382&amp;NEW_VIEW=on&amp;NUM_OVERVIEW=100&amp;PAGE=Nucleotides&amp;PROGRAM=blastn&amp;QUERY=&amp;QUERY_LENGTH=&amp;REPEATS=repeat_9606&amp;RID=1167160749-11191-152202340454.BLASTQ1&amp;RTOE=9&amp;SEARCH_NAME=short_bn&amp;SERVICE=plain&amp;SET_DEFAULTS.x=48&amp;SET_DEFAULTS.y=8&amp;SHOW_LINKOUT=on&amp;SHOW_OVERVIEW=on&amp;USER_TYPE=2&amp;WORD_SIZE=7&amp;dbtype=hc&amp;END_OF_HTTPGET=Yes" target="_blank">BLAST!</a>},

  # The maximal primer size primer3 can handle.
  "MAX_PRIMER_SIZE" =>  "36",

  # The maximal number of primers primer3Manager can handle.
  "MAX_NUMBER_PRIMER_MANAGER" =>  "100000",

  # The path were Primer3Manager stores its cached files.
  "USER_CACHE_FILES_PATH" =>  "./cached_data/",

  # The maximal time Primer3Manager will store the cached files / cookies.
  "MAX_STORAGE_TIME" =>  "+4d",

  # If you make any substantial modifications give this code a new
  # version designation.
  "CGI_VERSION" => "(primer3_www_results.cgi v 0.5)"
);

  # Add mispriming / mishybing libraries; 
my %misLibrary = (
  "NONE"              => "",
  "HUMAN"             => "humrep_and_simple.txt",
  "RODENT_AND_SIMPLE" => "rodrep_and_simple.txt",
  "RODENT"            => "rodent_ref.txt",
  "DROSOPHILA"        => "drosophila.w.transposons.txt"
  # Put more repeat libraries here. Add them also to the array.
);

my @misLibraryList = (
  "NONE",
  "HUMAN",
  "RODENT_AND_SIMPLE",
  "RODENT",
  "DROSOPHILA",
  # Put more repeat libraries here. Add them also to the hash.
);

  # Add sever stored Setting Files here: 
my %serverParameterFiles = (
  "Default"               => "",
  "qPCR"                  => "qPCR.txt",
  "Probe"                 => "probe.txt",
  "Use_Product_Size"      => "productSize.txt",
  # Put more Setting Files here. Add them also to the array.
);

my @serverParameterFilesList = (
  "Default",
  "qPCR",
  "Probe",  
  "Use_Product_Size",
  # Put more Setting Files here. Add them also to the hash.
);

##############################################################################
# ---------------- End Installer Modifiable Variables ---------------------- #
##############################################################################

# List of parameters to run primer3
my %defaultSettings = (
# Begin Primer3 Input Parameters
# Primer3 "Sequence" Input Tags
  "PRIMER_SEQUENCE_ID"                       => "",
  "SEQUENCE"                                 => "", 
  "INCLUDED_REGION"                          => "",
  "TARGET"                                   => "",
  "EXCLUDED_REGION"                          => "",
  "PRIMER_SEQUENCE_QUALITY"                  => "",
  "PRIMER_LEFT_INPUT"                        => "",
  "PRIMER_RIGHT_INPUT"                       => "",
  "PRIMER_START_CODON_POSITION"              => "", #"-1000000",
# Primer3 "Global" Input Tags
  "PRIMER_TASK"                              => "pick_pcr_primers",
  "PRIMER_PICK_ANYWAY"	                    => "1",
  "PRIMER_EXPLAIN_FLAG"                      => "0",
  "PRIMER_FILE_FLAG"                         => "0",
  "PRIMER_MISPRIMING_LIBRARY"                => "NONE",
  "PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS"     => "1",
  "PRIMER_MAX_MISPRIMING"                    => "12.00",
  "PRIMER_MAX_TEMPLATE_MISPRIMING"           => "-1.00",
  "PRIMER_PAIR_MAX_MISPRIMING"               => "24.00",
  "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING"      => "-1.00",
  "PRIMER_PRODUCT_MIN_TM"                    => "", #"-1000000.0",
  "PRIMER_PRODUCT_OPT_TM"                    => "",
  "PRIMER_PRODUCT_MAX_TM"                    => "", #"1000000.0",
  "PRIMER_PRODUCT_OPT_SIZE"                  => "0",
  "PRIMER_PRODUCT_SIZE_RANGE"                => "150-250 100-300 301-400 401-500 501-600 601-700 701-850 851-1000",  
  "PRIMER_GC_CLAMP"                          => "0",
  "PRIMER_OPT_SIZE"                          => "20",
  "PRIMER_MIN_SIZE"                          => "18",
  "PRIMER_MAX_SIZE"                          => "27",
  "PRIMER_OPT_TM"                            => "60.0",
  "PRIMER_MIN_TM"                            => "57.0",
  "PRIMER_MAX_TM"                            => "63.0",
  "PRIMER_MAX_DIFF_TM"                       => "100.0",
  "PRIMER_MIN_GC"                            => "20.0",
  "PRIMER_OPT_GC_PERCENT"                    => "",
  "PRIMER_MAX_GC"                            => "80.0",
  "PRIMER_SALT_CONC"                         => "50.0",
  "PRIMER_DNA_CONC"                          => "50.0",
  "PRIMER_NUM_NS_ACCEPTED"                   => "0",
  "PRIMER_SELF_ANY"                          => "8.00",
  "PRIMER_SELF_END"                          => "3.00",
  "PRIMER_MAX_POLY_X"                        => "5",
  "PRIMER_LIBERAL_BASE"                      => "1",
  "PRIMER_NUM_RETURN"                        => "5",
  "PRIMER_FIRST_BASE_INDEX"                  => "1",
  "PRIMER_MIN_QUALITY"                       => "0",
  "PRIMER_MIN_END_QUALITY"                   => "0",
  "PRIMER_QUALITY_RANGE_MIN"                 => "0",
  "PRIMER_QUALITY_RANGE_MAX"                 => "100",
  "PRIMER_INSIDE_PENALTY"                    => "",
  "PRIMER_OUTSIDE_PENALTY"                   => "0",
  "PRIMER_MAX_END_STABILITY"                 => "9.0",
  "PRIMER_WT_TM_GT"                          => "1.0",
  "PRIMER_WT_TM_LT"                          => "1.0",
  "PRIMER_WT_SIZE_LT"                        => "1.0",
  "PRIMER_WT_SIZE_GT"                        => "1.0",
  "PRIMER_WT_GC_PERCENT_LT"                  => "0.0",
  "PRIMER_WT_GC_PERCENT_GT"                  => "0.0",
  "PRIMER_WT_COMPL_ANY"                      => "0.0",
  "PRIMER_WT_COMPL_END"                      => "0.0",
  "PRIMER_WT_NUM_NS"                         => "0.0",
  "PRIMER_WT_REP_SIM"                        => "0.0",
  "PRIMER_WT_SEQ_QUAL"                       => "0.0",
  "PRIMER_WT_END_QUAL"                       => "0.0",
  "PRIMER_WT_POS_PENALTY"                    => "0.0",
  "PRIMER_WT_END_STABILITY"                  => "0.0",
  "PRIMER_WT_TEMPLATE_MISPRIMING"            => "",
  "PRIMER_PAIR_WT_PR_PENALTY"                => "1.0",
  "PRIMER_PAIR_WT_IO_PENALTY"                => "0.0",
  "PRIMER_PAIR_WT_DIFF_TM"                   => "0.0",
  "PRIMER_PAIR_WT_COMPL_ANY"                 => "0.0",
  "PRIMER_PAIR_WT_COMPL_END"                 => "0.0",
  "PRIMER_PAIR_WT_PRODUCT_TM_LT"             => "0.0",
  "PRIMER_PAIR_WT_PRODUCT_TM_GT"             => "0.0",
  "PRIMER_PAIR_WT_PRODUCT_SIZE_GT"           => "0.0",
  "PRIMER_PAIR_WT_PRODUCT_SIZE_LT"           => "0.0",
  "PRIMER_PAIR_WT_REP_SIM"                   => "0.0",
  "PRIMER_PAIR_WT_TEMPLATE_MISPRIMING"       => "",
# Primer3 Internal Oligo "Sequence" Input Tags
  "PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION"    => "",
  "PRIMER_INTERNAL_OLIGO_INPUT"              => "",
# Primer3 Internal Oligo "Global" Input Tags
  "PRIMER_INTERNAL_OLIGO_OPT_SIZE"           => "20",
  "PRIMER_INTERNAL_OLIGO_MIN_SIZE"           => "18",
  "PRIMER_INTERNAL_OLIGO_MAX_SIZE"           => "27",
  "PRIMER_INTERNAL_OLIGO_OPT_TM"             => "60.0",
  "PRIMER_INTERNAL_OLIGO_MIN_TM"             => "57.0",
  "PRIMER_INTERNAL_OLIGO_MAX_TM"             => "63.0",
  "PRIMER_INTERNAL_OLIGO_MIN_GC"             => "20.0",
  "PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT"     => "",
  "PRIMER_INTERNAL_OLIGO_MAX_GC"             => "80.0",
  "PRIMER_INTERNAL_OLIGO_SALT_CONC"          => "50.0",
  "PRIMER_INTERNAL_OLIGO_DNA_CONC"           => "50.0",
  "PRIMER_INTERNAL_OLIGO_SELF_ANY"           => "12.00",
  "PRIMER_INTERNAL_OLIGO_MAX_POLY_X"         => "5",
  "PRIMER_INTERNAL_OLIGO_SELF_END"           => "12.00",
  "PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY"     => "NONE",   
  "PRIMER_INTERNAL_OLIGO_MAX_MISHYB"         => "12.00",
  "PRIMER_INTERNAL_OLIGO_MIN_QUALITY"        => "0",
  "PRIMER_INTERNAL_OLIGO_NUM_NS"             => "0",
  "PRIMER_IO_WT_TM_GT"                       => "1.0",
  "PRIMER_IO_WT_TM_LT"                       => "1.0",
  "PRIMER_IO_WT_SIZE_LT"                     => "1.0",
  "PRIMER_IO_WT_SIZE_GT"                     => "1.0", 
  "PRIMER_IO_WT_GC_PERCENT_LT"               => "0.0",
  "PRIMER_IO_WT_GC_PERCENT_GT"               => "0.0",
  "PRIMER_IO_WT_COMPL_ANY"                   => "0.0",
  "PRIMER_IO_WT_NUM_NS"                      => "0.0",
  "PRIMER_IO_WT_REP_SIM"                     => "0.0",
  "PRIMER_IO_WT_SEQ_QUAL"                    => "0.0",
                            
# End of Primer3 Input Parameters

# Script Parameters
  "SCRIPT_TASK"                              => "Detection",
  "SCRIPT_PRINT_INPUT"                       => "0",
  "SCRIPT_FIX_PRIMER_END"                    => "5",
  
  "SCRIPT_CONTAINS_JAVA_SCRIPT"              => "1",

  "SCRIPT_SEQUENCING_LEAD"                   => "50",
  "SCRIPT_SEQUENCING_SPACING"                => "500",
  "SCRIPT_SEQUENCING_REVERSE"                => "1",
  "SCRIPT_SEQUENCING_INTERVAL"               => "250",
  "SCRIPT_SEQUENCING_ACCURACY"               => "20",

  "SCRIPT_DETECTION_PICK_LEFT"               => "1",
  "SCRIPT_DETECTION_PICK_HYB_PROBE"          => "0",
  "SCRIPT_DETECTION_PICK_RIGHT"              => "1",
  "SCRIPT_DETECTION_USE_PRODUCT_SIZE"        => "0",
  "SCRIPT_DETECTION_PRODUCT_MIN_SIZE"        => "100",
  "SCRIPT_DETECTION_PRODUCT_OPT_SIZE"        => "200",
  "SCRIPT_DETECTION_PRODUCT_MAX_SIZE"        => "1000",
  
  "SERVER_PARAMETER_FILE"                    => "DEFAULT",

  "PRIMER_NAME_ACRONYM_LEFT"                 => "F",
  "PRIMER_NAME_ACRONYM_INTERNAL_OLIGO"       => "IN",
  "PRIMER_NAME_ACRONYM_RIGHT"                => "R",
  "PRIMER_NAME_ACRONYM_SPACER"               => "_"
# if you ad parameters here also add them to the respective save array
);

# Array for the tasks Primer3plus can do (needed to build the HTML)
my @scriptTasks = (
  "Detection",
  "Cloning",
  "Sequencing",
  "Primer_List",
  "Primer_Check"
);

# All the parameters saved as a sequence file
my @sequenceParametersToSaveInFile =  qw(PRIMER_SEQUENCE_ID
SEQUENCE
INCLUDED_REGION
TARGET
EXCLUDED_REGION
PRIMER_SEQUENCE_QUALITY
PRIMER_LEFT_INPUT
PRIMER_RIGHT_INPUT
PRIMER_START_CODON_POSITION
PRIMER_INTERNAL_OLIGO_INPUT
PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION);

# All the parameters saved in a settings file
my @settingsParametersToSaveInFile = qw(PRIMER_TASK
PRIMER_PICK_ANYWAY
PRIMER_MISPRIMING_LIBRARY
PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS
PRIMER_MAX_MISPRIMING
PRIMER_MAX_TEMPLATE_MISPRIMING
PRIMER_PAIR_MAX_MISPRIMING
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING
PRIMER_PRODUCT_MIN_TM
PRIMER_PRODUCT_OPT_TM
PRIMER_PRODUCT_MAX_TM
PRIMER_PRODUCT_OPT_SIZE
PRIMER_PRODUCT_SIZE_RANGE
PRIMER_GC_CLAMP
PRIMER_OPT_SIZE
PRIMER_MIN_SIZE
PRIMER_MAX_SIZE
PRIMER_OPT_TM
PRIMER_MIN_TM
PRIMER_MAX_TM
PRIMER_MAX_DIFF_TM
PRIMER_MIN_GC
PRIMER_OPT_GC_PERCENT
PRIMER_MAX_GC
PRIMER_SALT_CONC
PRIMER_DNA_CONC
PRIMER_NUM_NS_ACCEPTED
PRIMER_SELF_ANY
PRIMER_SELF_END
PRIMER_MAX_POLY_X
PRIMER_LIBERAL_BASE
PRIMER_NUM_RETURN
PRIMER_FIRST_BASE_INDEX
PRIMER_MIN_QUALITY
PRIMER_MIN_END_QUALITY
PRIMER_QUALITY_RANGE_MIN
PRIMER_QUALITY_RANGE_MAX
PRIMER_INSIDE_PENALTY
PRIMER_OUTSIDE_PENALTY
PRIMER_MAX_END_STABILITY
PRIMER_WT_TM_GT
PRIMER_WT_TM_LT
PRIMER_WT_SIZE_LT
PRIMER_WT_SIZE_GT
PRIMER_WT_GC_PERCENT_LT
PRIMER_WT_GC_PERCENT_GT
PRIMER_WT_COMPL_ANY
PRIMER_WT_COMPL_END
PRIMER_WT_NUM_NS
PRIMER_WT_REP_SIM
PRIMER_WT_SEQ_QUAL
PRIMER_WT_END_QUAL
PRIMER_WT_POS_PENALTY
PRIMER_WT_END_STABILITY
PRIMER_WT_TEMPLATE_MISPRIMING
PRIMER_PAIR_WT_PR_PENALTY
PRIMER_PAIR_WT_IO_PENALTY
PRIMER_PAIR_WT_DIFF_TM
PRIMER_PAIR_WT_COMPL_ANY
PRIMER_PAIR_WT_COMPL_END
PRIMER_PAIR_WT_PRODUCT_TM_LT
PRIMER_PAIR_WT_PRODUCT_TM_GT
PRIMER_PAIR_WT_PRODUCT_SIZE_GT
PRIMER_PAIR_WT_PRODUCT_SIZE_LT
PRIMER_PAIR_WT_REP_SIM
PRIMER_PAIR_WT_TEMPLATE_MISPRIMING
PRIMER_INTERNAL_OLIGO_OPT_SIZE
PRIMER_INTERNAL_OLIGO_MIN_SIZE
PRIMER_INTERNAL_OLIGO_MAX_SIZE
PRIMER_INTERNAL_OLIGO_OPT_TM
PRIMER_INTERNAL_OLIGO_MIN_TM
PRIMER_INTERNAL_OLIGO_MAX_TM
PRIMER_INTERNAL_OLIGO_MIN_GC
PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT
PRIMER_INTERNAL_OLIGO_MAX_GC
PRIMER_INTERNAL_OLIGO_SALT_CONC
PRIMER_INTERNAL_OLIGO_DNA_CONC
PRIMER_INTERNAL_OLIGO_SELF_ANY
PRIMER_INTERNAL_OLIGO_MAX_POLY_X
PRIMER_INTERNAL_OLIGO_SELF_END
PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY
PRIMER_INTERNAL_OLIGO_MAX_MISHYB
PRIMER_INTERNAL_OLIGO_MIN_QUALITY
PRIMER_INTERNAL_OLIGO_NUM_NS
PRIMER_IO_WT_TM_GT
PRIMER_IO_WT_TM_LT
PRIMER_IO_WT_SIZE_LT
PRIMER_IO_WT_SIZE_GT
PRIMER_IO_WT_GC_PERCENT_LT
PRIMER_IO_WT_GC_PERCENT_GT
PRIMER_IO_WT_COMPL_ANY
PRIMER_IO_WT_NUM_NS
PRIMER_IO_WT_REP_SIM
PRIMER_IO_WT_SEQ_QUAL
SCRIPT_TASK
SCRIPT_PRINT_INPUT
SCRIPT_FIX_PRIMER_END
SCRIPT_CONTAINS_JAVA_SCRIPT
SCRIPT_SEQUENCING_LEAD
SCRIPT_SEQUENCING_SPACING
SCRIPT_SEQUENCING_REVERSE
SCRIPT_SEQUENCING_INTERVAL
SCRIPT_SEQUENCING_ACCURACY
SCRIPT_DETECTION_PICK_LEFT
SCRIPT_DETECTION_PICK_HYB_PROBE
SCRIPT_DETECTION_PICK_RIGHT
SCRIPT_DETECTION_USE_PRODUCT_SIZE
SCRIPT_DETECTION_PRODUCT_MIN_SIZE
SCRIPT_DETECTION_PRODUCT_OPT_SIZE
SCRIPT_DETECTION_PRODUCT_MAX_SIZE
SERVER_PARAMETER_FILE
PRIMER_NAME_ACRONYM_LEFT
PRIMER_NAME_ACRONYM_INTERNAL_OLIGO
PRIMER_NAME_ACRONYM_RIGHT
PRIMER_NAME_ACRONYM_SPACER);


# All the parameters primer3 can handle
my @primer3CompleteParameters =  qw(PRIMER_SEQUENCE_ID
SEQUENCE
INCLUDED_REGION
TARGET
EXCLUDED_REGION
PRIMER_SEQUENCE_QUALITY
PRIMER_LEFT_INPUT
PRIMER_RIGHT_INPUT
PRIMER_START_CODON_POSITION
PRIMER_TASK
PRIMER_PICK_ANYWAY
PRIMER_EXPLAIN_FLAG
PRIMER_FILE_FLAG
PRIMER_MISPRIMING_LIBRARY
PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS
PRIMER_MAX_MISPRIMING
PRIMER_MAX_TEMPLATE_MISPRIMING
PRIMER_PAIR_MAX_MISPRIMING
PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING
PRIMER_PRODUCT_MIN_TM
PRIMER_PRODUCT_OPT_TM
PRIMER_PRODUCT_MAX_TM
PRIMER_PRODUCT_OPT_SIZE
PRIMER_PRODUCT_SIZE_RANGE
PRIMER_GC_CLAMP
PRIMER_OPT_SIZE
PRIMER_MIN_SIZE
PRIMER_MAX_SIZE
PRIMER_OPT_TM
PRIMER_MIN_TM
PRIMER_MAX_TM
PRIMER_MAX_DIFF_TM
PRIMER_MIN_GC
PRIMER_OPT_GC_PERCENT
PRIMER_MAX_GC
PRIMER_SALT_CONC
PRIMER_DNA_CONC
PRIMER_NUM_NS_ACCEPTED
PRIMER_SELF_ANY
PRIMER_SELF_END
PRIMER_MAX_POLY_X
PRIMER_LIBERAL_BASE
PRIMER_NUM_RETURN
PRIMER_FIRST_BASE_INDEX
PRIMER_MIN_QUALITY
PRIMER_MIN_END_QUALITY
PRIMER_QUALITY_RANGE_MIN
PRIMER_QUALITY_RANGE_MAX
PRIMER_INSIDE_PENALTY
PRIMER_OUTSIDE_PENALTY
PRIMER_MAX_END_STABILITY
PRIMER_WT_TM_GT
PRIMER_WT_TM_LT
PRIMER_WT_SIZE_LT
PRIMER_WT_SIZE_GT
PRIMER_WT_GC_PERCENT_LT
PRIMER_WT_GC_PERCENT_GT
PRIMER_WT_COMPL_ANY
PRIMER_WT_COMPL_END
PRIMER_WT_NUM_NS
PRIMER_WT_REP_SIM
PRIMER_WT_SEQ_QUAL
PRIMER_WT_END_QUAL
PRIMER_WT_POS_PENALTY
PRIMER_WT_END_STABILITY
PRIMER_WT_TEMPLATE_MISPRIMING
PRIMER_PAIR_WT_PR_PENALTY
PRIMER_PAIR_WT_IO_PENALTY
PRIMER_PAIR_WT_DIFF_TM
PRIMER_PAIR_WT_COMPL_ANY
PRIMER_PAIR_WT_COMPL_END
PRIMER_PAIR_WT_PRODUCT_TM_LT
PRIMER_PAIR_WT_PRODUCT_TM_GT
PRIMER_PAIR_WT_PRODUCT_SIZE_GT
PRIMER_PAIR_WT_PRODUCT_SIZE_LT
PRIMER_PAIR_WT_REP_SIM
PRIMER_PAIR_WT_TEMPLATE_MISPRIMING
PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION
PRIMER_INTERNAL_OLIGO_INPUT
PRIMER_INTERNAL_OLIGO_OPT_SIZE
PRIMER_INTERNAL_OLIGO_MIN_SIZE
PRIMER_INTERNAL_OLIGO_MAX_SIZE
PRIMER_INTERNAL_OLIGO_OPT_TM
PRIMER_INTERNAL_OLIGO_MIN_TM
PRIMER_INTERNAL_OLIGO_MAX_TM
PRIMER_INTERNAL_OLIGO_MIN_GC
PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT
PRIMER_INTERNAL_OLIGO_MAX_GC
PRIMER_INTERNAL_OLIGO_SALT_CONC
PRIMER_INTERNAL_OLIGO_DNA_CONC
PRIMER_INTERNAL_OLIGO_SELF_ANY
PRIMER_INTERNAL_OLIGO_MAX_POLY_X
PRIMER_INTERNAL_OLIGO_SELF_END
PRIMER_INTERNAL_OLIGO_MISHYB_LIBRARY
PRIMER_INTERNAL_OLIGO_MAX_MISHYB
PRIMER_INTERNAL_OLIGO_MIN_QUALITY
PRIMER_INTERNAL_OLIGO_NUM_NS
PRIMER_IO_WT_TM_GT
PRIMER_IO_WT_TM_LT
PRIMER_IO_WT_SIZE_LT
PRIMER_IO_WT_SIZE_GT
PRIMER_IO_WT_GC_PERCENT_LT
PRIMER_IO_WT_GC_PERCENT_GT
PRIMER_IO_WT_COMPL_ANY
PRIMER_IO_WT_NUM_NS
PRIMER_IO_WT_REP_SIM
PRIMER_IO_WT_SEQ_QUAL);

# Parameters send to primer3 using primer_check
my @primer3PrimerCheckParameters =  qw(PRIMER_SEQUENCE_ID
SEQUENCE
PRIMER_SEQUENCE_QUALITY
PRIMER_LEFT_INPUT
PRIMER_TASK
PRIMER_PICK_ANYWAY
PRIMER_EXPLAIN_FLAG
PRIMER_FILE_FLAG
PRIMER_MISPRIMING_LIBRARY
PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS
PRIMER_MAX_MISPRIMING
PRIMER_MAX_TEMPLATE_MISPRIMING
PRIMER_GC_CLAMP
PRIMER_OPT_SIZE
PRIMER_MIN_SIZE
PRIMER_MAX_SIZE
PRIMER_OPT_TM
PRIMER_MIN_TM
PRIMER_MAX_TM
PRIMER_MIN_GC
PRIMER_OPT_GC_PERCENT
PRIMER_MAX_GC
PRIMER_SALT_CONC
PRIMER_DNA_CONC
PRIMER_NUM_NS_ACCEPTED
PRIMER_SELF_ANY
PRIMER_SELF_END
PRIMER_MAX_POLY_X
PRIMER_LIBERAL_BASE
PRIMER_NUM_RETURN
PRIMER_MIN_QUALITY
PRIMER_MIN_END_QUALITY
PRIMER_QUALITY_RANGE_MIN
PRIMER_QUALITY_RANGE_MAX
PRIMER_MAX_END_STABILITY
PRIMER_WT_TM_GT
PRIMER_WT_TM_LT
PRIMER_WT_SIZE_LT
PRIMER_WT_SIZE_GT
PRIMER_WT_GC_PERCENT_LT
PRIMER_WT_GC_PERCENT_GT
PRIMER_WT_COMPL_ANY
PRIMER_WT_COMPL_END
PRIMER_WT_NUM_NS
PRIMER_WT_REP_SIM
PRIMER_WT_SEQ_QUAL
PRIMER_WT_END_QUAL
PRIMER_WT_POS_PENALTY
PRIMER_WT_END_STABILITY);

# Primer3 does not add a 0 to the first parameters
# this are the parameters to replace:

my %zeroReplacements= (
  "PRIMER_LEFT"                           => "PRIMER_LEFT_0",
  "PRIMER_RIGHT"                          => "PRIMER_RIGHT_0",
  "PRIMER_INTERNAL_OLIGO"                 => "PRIMER_INTERNAL_OLIGO_0",
  "PRIMER_PRODUCT_SIZE"                   => "PRIMER_PRODUCT_SIZE_0",
  "PRIMER_PAIR_PENALTY"                   => "PRIMER_PAIR_0_PENALTY",
  "PRIMER_LEFT_PENALTY"                   => "PRIMER_LEFT_0_PENALTY",
  "PRIMER_RIGHT_PENALTY"                  => "PRIMER_RIGHT_0_PENALTY",
  "PRIMER_INTERNAL_OLIGO_PENALTY"         => "PRIMER_INTERNAL_OLIGO_0_PENALTY",
  "PRIMER_LEFT_SEQUENCE"                  => "PRIMER_LEFT_0_SEQUENCE",
  "PRIMER_RIGHT_SEQUENCE"                 => "PRIMER_RIGHT_0_SEQUENCE",
  "PRIMER_INTERNAL_OLIGO_SEQUENCE"        => "PRIMER_INTERNAL_OLIGO_0_SEQUENCE",
  "PRIMER_LEFT_TM"                        => "PRIMER_LEFT_0_TM",
  "PRIMER_RIGHT_TM"                       => "PRIMER_RIGHT_0_TM",
  "PRIMER_INTERNAL_OLIGO_TM"              => "PRIMER_INTERNAL_OLIGO_0_TM",
  "PRIMER_LEFT_GC_PERCENT"                => "PRIMER_LEFT_0_GC_PERCENT",
  "PRIMER_RIGHT_GC_PERCENT"               => "PRIMER_RIGHT_0_GC_PERCENT",
  "PRIMER_INTERNAL_OLIGO_GC_PERCENT"      => "PRIMER_INTERNAL_OLIGO_0_GC_PERCENT",
  "PRIMER_LEFT_SELF_ANY"                  => "PRIMER_LEFT_0_SELF_ANY",
  "PRIMER_RIGHT_SELF_ANY"                 => "PRIMER_RIGHT_0_SELF_ANY",
  "PRIMER_INTERNAL_OLIGO_SELF_ANY"        => "PRIMER_INTERNAL_OLIGO_0_SELF_ANY",
  "PRIMER_LEFT_SELF_END"                  => "PRIMER_LEFT_0_SELF_END",
  "PRIMER_RIGHT_SELF_END"                 => "PRIMER_RIGHT_0_SELF_END",
  "PRIMER_INTERNAL_OLIGO_SELF_END"        => "PRIMER_INTERNAL_OLIGO_0_SELF_END",
  "PRIMER_PAIR_COMPL_ANY"                 => "PRIMER_PAIR_0_COMPL_ANY",
  "PRIMER_PAIR_COMPL_END"                 => "PRIMER_PAIR_0_COMPL_END",
  "PRIMER_LEFT_MISPRIMING_SCORE"          => "PRIMER_LEFT_0_MISPRIMING_SCORE",
  "PRIMER_RIGHT_MISPRIMING_SCORE"         => "PRIMER_RIGHT_0_MISPRIMING_SCORE",
  "PRIMER_PAIR_MISPRIMING_SCORE"          => "PRIMER_PAIR_0_MISPRIMING_SCORE",
  "PRIMER_LEFT_TEMPLATE_MISPRIMING"       => "PRIMER_LEFT_0_TEMPLATE_MISPRIMING",
  "PRIMER_RIGHT_TEMPLATE_MISPRIMING"      => "PRIMER_RIGHT_0_TEMPLATE_MISPRIMING",
  "PRIMER_PAIR_TEMPLATE_MISPRIMING"       => "PRIMER_PAIR_0_TEMPLATE_MISPRIMING",
  "PRIMER_PRODUCT_TM"                     => "PRIMER_PRODUCT_0_TM",
  "PRIMER_PRODUCT_TM_OLIGO_TM_DIF"        => "PRIMER_PRODUCT_0_TM_OLIGO_TM_DIF",
  "PRIMER_PAIR_T_OPT_A"                   => "PRIMER_PAIR_0_T_OPT_A",
  "PRIMER_INTERNAL_OLIGO_MISHYB_SCORE"    => "PRIMER_INTERNAL_OLIGO_0_MISHYB_SCORE",
  "PRIMER_LEFT_MIN_SEQ_QUALITY"           => "PRIMER_LEFT_0_MIN_SEQ_QUALITY",
  "PRIMER_RIGHT_MIN_SEQ_QUALITY"          => "PRIMER_RIGHT_0_MIN_SEQ_QUALITY",
  "PRIMER_INTERNAL_OLIGO_MIN_SEQ_QUALITY" => "PRIMER_INTERNAL_OLIGO_0_MIN_SEQ_QUALITY",
  "PRIMER_LEFT_END_STABILITY"             => "PRIMER_LEFT_0_END_STABILITY",
  "PRIMER_RIGHT_END_STABILITY"            => "PRIMER_RIGHT_0_END_STABILITY" 
);

################################################
# Methods for getting the hashes or the arrays #
################################################

sub getMachineSetting {
  my ($key, $value);
  $key = shift;
  if ($key eq "") {
  	setMessage("Programming-error: Function \"getMachineSetting\" needs a Key");
  }
  else {
  	$value = $machineSettings{$key};
  }

  return $value;
}

sub getMachineSettings {
  return %machineSettings;
}

sub getMisLibrary {
  return %misLibrary;
}

sub getLibraryList {
  return @misLibraryList;
}

sub getServerParameterFiles {
  return %serverParameterFiles;
}

sub getServerParameterFilesList {
  return @serverParameterFilesList;
}

sub getDefaultSettings {
  return %defaultSettings;
}

sub getScriptTask {
  return @scriptTasks;
}

sub getSaveSequenceParameters {
  return @sequenceParametersToSaveInFile;
}

sub getSaveSettingsParameters {
  return @settingsParametersToSaveInFile;
}

sub getPrimer3CompleteParameters {
  return @primer3CompleteParameters;
}

sub getPrimerCheckParameters {
  return @primer3PrimerCheckParameters;
}

sub getZeroReplacements {
  return %zeroReplacements;
}

############################################
# Handling from Messages and Errormessages #
############################################

sub setMessage {
  my $note = shift;
  push @messages, $note;

  return;
}      

sub getMessages {
  return @messages;
}

#################################
# To Stop Primer3plus at errors #
#################################
sub setDoNotPick {
  my $value;
  $value = shift;
  if ($value eq 0) {
  	$doNotPick = 0;
  }
  else {
  	$doNotPick = 1;
  }

  return;
}

sub getDoNotPick {
  return $doNotPick;
}


1;
