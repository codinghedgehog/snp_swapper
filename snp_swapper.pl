#!/usr/bin/perl
#
# SNP Swapper
#
# This script reconstructs the base sequence using a base reference sequence in FASTA format, and a SNP loci file and indel file
# produced by prephix.  It then does substituting, inserting, or deleting of the base sequence bases at the given locations in 
# the provided reference base sequence. It writes out a separate regenerated base sequence for each strain in the SNP input file.
#
# The output is suitable for use by mlstar (i.e. is a FASTA formatted reference base sequence).
#
# The reference base sequence input file should be in FASTA format.  It is assumed to to start at loci 1 with the first base.
#
# The SNP file input should contain three TAB-delimited columns and no headers:
# SNP_ID (i.e. strain id) [TAB] Base position [TAB] Base
#
# So something like:
# A12	1045	G
# A12	4056	A
# A12 13004	T
# A35	4	A
# A35	401	C
#
# This is the same format of the snp files produced by the prephix program.
#
# The indel input file is in the format  produced by the prephix program and contain modified VAAL4 K28 and NUCMER lines.  
# The modification include strain ID and either k28 or nuc is in the first and second columns.
# The indel input files have no header files, e.g.:
#
# STRAIN_ID k28 0 316 left=CAGGTATTTGACATATAGAG sample=A ref=G right=ACTGAAAAAGTATAATTGTG
# STRAIN_ID k28 0 419 left=CTGTGCATAACTAATAAGCA sample= ref=ACG right=GATAAAGTTATCCACCGATT
# STRAIN_ID k28 0 929 left=GACACTTTTGTAATCGGACC sample= ref=C right=GGTAACCGCTTTCCACATGC
# STRAIN_ID k28 0 953 left=AACCGCTTTCCACATGCAGC sample=A ref= right=AGTTTAGCTGTGGCCGAAGC
# STRAIN_ID k28 0 965 left=CATGCAGCGAGTTTAGCTGT sample=AAT ref= right=GCCGAAGCACCAGCCAAAGC
# STRAIN_ID k28 0 1013 left=CCATTATTTATCTATGGAGG sample=G ref= right=GTTGGTTTAGGAAAAACCCA
# STRAIN_ID nuc 759437  A G 732302  1 732302  1 1 NC007793  JKD6159
# STRAIN_ID nuc 759441  T A 732306  3 732306  1 1 NC007793  JKD6159
# STRAIN_ID nuc 759444  T C 732309  3 732309  1 1 NC007793  JKD6159
# STRAIN_ID nuc 759456  G A 732321  6 732321  1 1 NC007793  JKD6159
# STRAIN_ID nuc 759462  A T 732327  6 732327  1 1 NC007793  JKD6159
# STRAIN_ID nuc 759504  A G 732369  36  732369  1 1 NC007793  JKD6159
# STRAIN_ID nuc 759540  T C 732405  15  732405  1 1 NC007793  JKD6159
#
# Usage: $0 <reference base file> <prephix SNP loci input file> [prephix indel file]
#
# 5/18/2012 AJP Updated to process multiple strains in snp file, and also updated to read new indel input file format
# from prephix (which was also modified to write this new format).  New format is the strain_id and indel line type column
# followed by either the k28 or nucmer indel file (discarded by prephix).

use strict;

my $VERSION="2.1 Rev 1";

print "\nSNP Swapper v$VERSION\n\n";

if ($#ARGV < 1){
  print " Usage: $0 <reference base file> <prephix SNP loci input file> [prephix indel file] [-debug] -quiet]\n";
  exit 1;
}

my $snpfile;
my $reffile;
my $outfile;
my $logfile;
my %baseRefTable;
my %outputTable;
my %snpTable;
my $i=0;
my $j=0;
my $currentStrain="";
my $loci;
my $strainCount=0;
my $charCount=0;
my $warnings=0;
my $debug="N";
my $quiet="N";
my $arg_num=0;
my $indelFilename;
my $indelfile;
my $realLoci;
my $in_strain;
my $in_loci;
my $in_base;

open($reffile,"<",$ARGV[0]) or die "Unable to open reference file $ARGV[0] for reading!  $!\n";
open($snpfile,"<",$ARGV[1]) or die "Unable to open snp input file $ARGV[1] for reading!  $!\n";
open($outfile,">","$ARGV[0].swapped") or die "unable to open output file $ARGV[0].swapped for writing! $!\n";
open($logfile,">","$ARGV[0].swapped.log") or die "unable to open log file $ARGV[0].swapped.log for writing! $!\n";

# Process optional parameters
if ($#ARGV > 1){
  $arg_num=2;
  while ($arg_num <= $#ARGV){
    # Using if-else logic because given/when isn't compatible with< 5.10 and Switch module was removed from >5.14.
    if ($ARGV[$arg_num] eq "-debug"){
				$debug="Y";
        print "Producing debug output.\n";
    }
    elsif ($ARGV[$arg_num] eq "-quiet") { 
        print "Producing quiet (no stdout) output.  Log file is still generated.\n";
				$quiet="Y";
    }
		elsif (-r $ARGV[$arg_num]){
		    $indelFilename = $ARGV[$arg_num];
		}
    else{
        print "*** ERROR: Unknown command line argument or unreadable file $ARGV[$arg_num].\n";
        exit 1;
    }
    $arg_num++;
  }
}  

print_all("Loading loci entries from reference FASTA file...\n");
# Read in and store reference data in a hash table keyed off of loci number.
# Note, data is actually an array consisting of base sequence and a bit to indicate if it has been modified by
# snp/indel processing.  This is to allow warnings to be generated if a snp or indel file references the same loci
# in the reference base sequence.

my $refID;
my $baseLoci=0;
my $line=0;
while (<$reffile>){
  $line++;
  chomp;

	# Save off the ID from the first line, which should be prefixed with >.
	if (/^>(.+)/){
		$refID=$1;
		print_debug("Found reference genome id $refID\n");
		next;
	}

	# The rest of the file should be lines of base sequences
	if (/^[AGCT]+$/){
    # Split into individual bases and load into reference has table, using auto-incremented counter as loci.
    my $refBase;
    foreach $refBase (split(//)){
      $baseLoci++;
			$baseRefTable{$baseLoci}=[$refBase,"not_modified"];
    }
  }
	else{
		print_all("*** ERROR: Invalid format in reference file on line $line (invalid base pair character?): $_\n");
		exit 1;
	}
}
print "\n";
close($reffile);
print_all("Read $baseLoci bases from reference file.\n");

print "Processing SNP file $ARGV[1]\n";
# For each strain in the input snp file, alter the base reference at the given loci to match the snp, then process its
# corresponding indel file (if any) and write out the swapped base sequence.
$i = 0;
%outputTable = %baseRefTable;
while (<$snpfile>){
  $i++;

  chomp;

  if (!(/^(\S+)\t(\d+)\t([A-Z]+)$/)){
    print_all("*** ERROR: Bad data format on line $i of SNP loci file.  Expect three columns [Strain ID (Alphanumeric)] [Loci (Numeric)] [Base (alphabetic)]\nGot \"$_\" instead.\n");
		print "Failed.\n";
    exit 1;
  }

  ($in_strain,$in_loci,$in_base) = split /\t/;

  print_debug("Diag: Read in $in_strain, $in_loci, $in_base\n");

	if ($currentStrain eq ""){
	  $currentStrain=$in_strain;
	}

	# If a new strain was detected, process for any indels for this strain then write it out.
	if ($currentStrain ne $in_strain){

    if ($indelFilename ne "" ){
			# For each matching strain entry in indel file, alter the base reference at the given loci to match the indel.
			$j = 0;
			my $loci_offset=0; # If reading indel file generated by prephix, the loci offset is 0. (k28 out file loci offset is 1).
			my $lineType="";
			my $indelStrain="";
			open ($indelfile,"<",$indelFilename) or die "Unable to open $indelFilename for processing! $!\n";
			while (<$indelfile>){
				$j++;

				chomp;
				print_debug("At line $_\n");

				# Below is recycled code from prephix for reading k28 files, which is basically the same
				# format the indel files that prephix outputs. 
				# Prephix adds two columns -- strain id and nuc|k28 to indicate the associated strain's indel and
				# the format of the line.  The rest of the line is just dumped out as-is from the input file passed
				# to prephix.

				# Get strain ID and line type from current line.  See if it matches the current strain ID.
				# If not, keep looking.
				if (/^(.+?)\t(nuc|k28)/){
					$indelStrain=$1;
					$lineType=$2;
					
					print_debug("Found strain id of $indelStrain in file $indelFilename\n");
					if ($indelStrain ne $currentStrain){
						print_debug("Skipping since it doesn't match $currentStrain\n");
						next;
					}

					# If it matches, then process it depending on whether it is k28.out or NUCMER format.
					if ($lineType eq "k28"){
						print_debug("Line $j appears to be k28.out/VAAL format.\n");
						$loci_offset=1;
						# Assuming indel input body data to be in the format:
						# strain_id k28 0 <snp_locus> <left flank seq> <sample> <ref> <right flank seq>
						# Only care about the locus and sample columns.
						#
						# Regex note:
						#
						# With indels, either sample= or ref= may have no value, so match for [ATCG] and check for length 1.
						if (/^\S+\s+k28\s+[0-9]+\s+([0-9]+)\s+left=[ATCG]*\s+sample=([ATCG]*)\s+ref=([ATCG]*)\s+right=[ATCG]*$/){
							$realLoci=$1;
							$realLoci += $loci_offset; # VAAL k28.out file loci is offset by +1.

							# If sample=something and ref=nothing, then this is an insert, so insert the ref base sequence at the given loci.
							if (( length($2) >= 1 ) && ( length($3) == 0) ){
								print_debug("Inserting $2 at $realLoci\n");
								if (exists($outputTable{$realLoci})){
									$outputTable{$realLoci} = ["$outputTable{$realLoci}[0]$2","modified"];
								}
								else{
									$outputTable{$realLoci} = [$2,"modified"];
								}
							}
							# If sample=nothing and ref=something, then this is a deletion, so delete the base at the loci in the ref base
							# sequence, and the following bases (as many as there are in the deletion sequence in the sample).
							# I.e. if ref=GGG then need to delete base at the given loci as well as the next two loci.
							elsif (( length($2) == 0 ) && ( length($3) >= 1) ){
								my $del_count=0;
								while ($del_count < length($3)){
									# Sanity check that the deleted sequence does indeed start at this loci (and subsequent loci matches the deleted bases).
									if ($outputTable{$realLoci + $del_count}[0] ne substr($3,$del_count,1)){
										print_all("*** ERROR: Encountered indel loci " . ($realLoci + $del_count) . ", for deletion which does not match expected base in ref base! Expected base " . substr($3,$del_count,1) . ", found $outputTable{$realLoci + $del_count}[0] instead. (line $j of indel input file $indelFilename). Aborting...\n");
										exit 1;
									}
									elsif ($outputTable{$realLoci + $del_count}[1] ne "not_modified"){
										print_all("*** WARNING: Encountered indel loci " . ($realLoci + $del_count) . ", which is already modified in base ref sequence! (line $j of indel input file $indelFilename) Skipping...\n");
										$warnings++;
									}
									else{
										print_debug("Deleting " . $outputTable{$realLoci + $del_count}[0]  . " at " . ($realLoci + $del_count) . "\n");
										$outputTable{$realLoci + $del_count} = ["","modified"];
									}
									$del_count++;
								}
							}
							else{
								print_all("*** ERROR: Indel file line format not recognized. Got \"$_\" at line $j. Can't figure if it is deletion or insertion (substitutions should be defined in the snp file).\n");
								print "Failed.\n";
								exit 1;
							}
						}
						else{
							print_all("*** ERROR: Bad k28 indel line format! $_\n");
							exit 1;
						}
					}
					elsif ($lineType eq "nuc"){
						print_debug("Line $j appears to be NUCMER format.\n");
						$loci_offset=0;
						# Assuming indel input body data to be in the format:
						# strain_id nuc [P1]  [SUB] [SUB] [P2]  [BUFF]  [DIST]  [LEN R] [LEN Q] [FRM] [TAGS]
						# Only care about the P1 (ref base loci), the first [SUB] which is the ref base at that loci,
						# and the second [SUB] which is the SNP base at that loci.
						#
						if (/^\S+\s+nuc\s+([0-9]+)\s+([ATCG]*)\s+([ATCG]*)\s+.+$/){
							$realLoci=$1;

							# If [SUB2] is something and [SUB1] is nothing, then this is an insert, so insert the ref base sequence at the given loci.
							if (( length($3) >= 1 ) && ( length($2) == 0) ){
								print_debug("Inserting $3 at $realLoci\n");
								if (exists($outputTable{$realLoci})){
									$outputTable{$realLoci} = ["$outputTable{$realLoci}[0]$3","modified"];
								}
								else{
									$outputTable{$realLoci} = [$3,"modified"];
								}
							}
							# If [SUB2] is nothing and [SUB1] is something, then this is a deletion, so delete the base at the loci in the ref base
							# sequence, and the following bases (as many as there are in the deletion sequence in the sample).
							# I.e. if ref=GGG then need to delete base at the given loci as well as the next two loci.
							elsif (( length($3) == 0 ) && ( length($2) >= 1) ){
								my $del_count=0;
								while ($del_count < length($2)){
									# Sanity check that the deleted sequence does indeed start at this loci (and subsequent loci matches the deleted bases).
									if ($outputTable{$realLoci + $del_count}[0] ne substr($2,$del_count,1)){
										print_all("*** ERROR: Encountered indel loci " . ($realLoci + $del_count) . ", for deletion which does not match expected base in ref base! Expected base " . substr($2,$del_count,1) . ", found $outputTable{$realLoci + $del_count}[0] instead. (line $j of indel input file $indelFilename). Aborting...\n");
										exit 1;
									}
									elsif ($outputTable{$realLoci + $del_count}[1] ne "not_modified"){
										print_all("*** WARNING: Encountered indel loci " . ($realLoci + $del_count) . ", which is already modified in base ref sequence! (line $j of indel input file $indelFilename) Skipping...\n");
										$warnings++;
									}
									else{
										print_debug("Deleting " . $outputTable{$realLoci + $del_count}[0]  . " at " . ($realLoci + $del_count) . "\n");
										$outputTable{$realLoci + $del_count} = ["","modified"];
									}
									$del_count++;
								}
							}
							else{
								print_all("*** ERROR: Indel file line format not recognized. Got \"$_\" at line $j. Can't figure if it is deletion or insertion (substitutions should be defined in the snp file).\n");
								print "Failed.\n";
								exit 1;
							}
						}
						else{
							print_all("*** ERROR: Bad NUMCER indel line format! $_\n");
							exit 1;
						}

					}
					else{
						print_all("Cannot determine type of indel line: $_ (in $indelFilename, line $j)\n");
						exit 1;
					}
				}
				else{
					print_all("Cannot parse indel line: $_ (in $indelFilename, line $j)\n");
					exit 1;
				}
			}
			close($indelfile);
    }
		# Dump out the modified base ref sequence in FASTA format.
		$charCount = 0;
		print $outfile ">$currentStrain\n";
		foreach $loci (sort {$a <=> $b} keys %outputTable){

			$charCount += length($outputTable{$loci}[0]);

			# Limit to 70 chars per line, for FASTA compliance.
			if ($charCount > 70){
				print $outfile "\n";
				$charCount = 0;
			}
			print $outfile "$outputTable{$loci}[0]";
		}
		print $outfile "\n";

		# Reset for new strain.
    %outputTable = %baseRefTable;
    $i = 1;
		$currentStrain=$in_strain;
	}
	
  if ($outputTable{$in_loci}[1] ne "not_modified"){
    print_all("*** WARNING: Encountered SNP loci $in_loci, which is already modified in base ref sequence! (line $i of SNP input file) Skipping...\n");
    $warnings++;
  }
  else{
		$outputTable{$in_loci} = [$in_base,"modified"];
  }
}

# Process any indels for the final strain.
# For each matching strain entry in indel file, alter the base reference at the given loci to match the indel.
if ($indelFilename ne ""){
	$j = 0;
	my $loci_offset=0; # If reading indel file generated by prephix, the loci offset is 0. (k28 out file loci offset is 1).
	my $lineType="";
	my $indelStrain="";
	open ($indelfile,"<",$indelFilename) or die "Unable to open $indelFilename for processing! $!\n";
	while (<$indelfile>){
		$j++;

		chomp;
		print_debug("At line $_\n");

		# Below is recycled code from prephix for reading k28 files, which is basically the same
		# format the indel files that prephix outputs. 
		# Prephix adds two columns -- strain id and nuc|k28 to indicate the associated strain's indel and
		# the format of the line.  The rest of the line is just dumped out as-is from the input file passed
		# to prephix.

		# Get strain ID and line type from current line.  See if it matches the current strain ID.
		# If not, keep looking.
		if (/^(.+?)\t(nuc|k28)/){
			$indelStrain=$1;
			$lineType=$2;
			
			print_debug("Found strain id of $indelStrain in file $indelFilename\n");
			if ($indelStrain ne $currentStrain){
				print_debug("Skipping since it doesn't match $currentStrain\n");
				next;
			}

			# If it matches, then process it depending on whether it is k28.out or NUCMER format.
			if ($lineType eq "k28"){
				print_debug("Line $j appears to be k28.out/VAAL format.\n");
				# Assuming indel input body data to be in the format:
				# strain_id k28 0 <snp_locus> <left flank seq> <sample> <ref> <right flank seq>
				# Only care about the locus and sample columns.
				#
				# Regex note:
				#
				# With indels, either sample= or ref= may have no value, so match for [ATCG] and check for length 1.
				if (/^\S+\s+k28\s+[0-9]+\s+([0-9]+)\s+left=[ATCG]*\s+sample=([ATCG]*)\s+ref=([ATCG]*)\s+right=[ATCG]*$/){
					$realLoci=$1;
					$realLoci += $loci_offset; # VAAL k28.out file loci is offset by +1, prephix indel file loci offset is 0.

					# If sample=something and ref=nothing, then this is an insert, so insert the ref base sequence at the given loci.
					if (( length($2) >= 1 ) && ( length($3) == 0) ){
						print_debug("Inserting $2 at $realLoci\n");
						if (exists($outputTable{$realLoci})){
							$outputTable{$realLoci} = ["$outputTable{$realLoci}[0]$2","modified"];
						}
						else{
							$outputTable{$realLoci} = [$2,"modified"];
						}
					}
					# If sample=nothing and ref=something, then this is a deletion, so delete the base at the loci in the ref base
					# sequence, and the following bases (as many as there are in the deletion sequence in the sample).
					# I.e. if ref=GGG then need to delete base at the given loci as well as the next two loci.
					elsif (( length($2) == 0 ) && ( length($3) >= 1) ){
						my $del_count=0;
						while ($del_count < length($3)){
							# Sanity check that the deleted sequence does indeed start at this loci (and subsequent loci matches the deleted bases).
							if ($outputTable{$realLoci + $del_count}[0] ne substr($3,$del_count,1)){
								print_all("*** ERROR: Encountered indel loci " . ($realLoci + $del_count) . ", for deletion which does not match expected base in ref base! Expected base " . substr($3,$del_count,1) . ", found $outputTable{$realLoci + $del_count}[0] instead. (line $j of indel input file $indelFilename). Aborting...\n");
								exit 1;
							}
							elsif ($outputTable{$realLoci + $del_count}[1] ne "not_modified"){
								print_all("*** WARNING: Encountered indel loci " . ($realLoci + $del_count) . ", which is already modified in base ref sequence! (line $j of indel input file $indelFilename) Skipping...\n");
								$warnings++;
							}
							else{
								print_debug("Deleting " . $outputTable{$realLoci + $del_count}[0]  . " at " . ($realLoci + $del_count) . "\n");
								$outputTable{$realLoci + $del_count} = ["","modified"];
							}
							$del_count++;
						}
					}
					else{
						print_all("*** ERROR: Indel file line format not recognized. Got \"$_\" at line $j. Can't figure if it is deletion or insertion (substitutions should be defined in the snp file).\n");
						print "Failed.\n";
						exit 1;
					}
				}
				else{
					print_all("*** ERROR: Bad k28 indel line format! $_\n");
					exit 1;
				}
			}
			elsif ($lineType eq "nuc"){
				print_debug("Line $j appears to be NUCMER format.\n");

				# Assuming indel input body data to be in the format:
				# strain_id nuc [P1]  [SUB] [SUB] [P2]  [BUFF]  [DIST]  [LEN R] [LEN Q] [FRM] [TAGS]
				# Only care about the P1 (ref base loci), the first [SUB] which is the ref base at that loci,
				# and the second [SUB] which is the SNP base at that loci.
				#
				if (/^\S+\s+nuc\s+([0-9]+)\s+([ATCG]*)\s+([ATCG]*)\s+.+$/){
					$realLoci=$1;

					# If [SUB2] is something and [SUB1] is nothing, then this is an insert, so insert the ref base sequence at the given loci.
					if (( length($3) >= 1 ) && ( length($2) == 0) ){
						print_debug("Inserting $3 at $realLoci\n");
						if (exists($outputTable{$realLoci})){
							$outputTable{$realLoci} = ["$outputTable{$realLoci}[0]$3","modified"];
						}
						else{
							$outputTable{$realLoci} = [$3,"modified"];
						}
					}
					# If [SUB2] is nothing and [SUB1] is something, then this is a deletion, so delete the base at the loci in the ref base
					# sequence, and the following bases (as many as there are in the deletion sequence in the sample).
					# I.e. if ref=GGG then need to delete base at the given loci as well as the next two loci.
					elsif (( length($3) == 0 ) && ( length($2) >= 1) ){
						my $del_count=0;
						while ($del_count < length($2)){
							# Sanity check that the deleted sequence does indeed start at this loci (and subsequent loci matches the deleted bases).
							if ($outputTable{$realLoci + $del_count}[0] ne substr($2,$del_count,1)){
								print_all("*** ERROR: Encountered indel loci " . ($realLoci + $del_count) . ", for deletion which does not match expected base in ref base! Expected base " . substr($2,$del_count,1) . ", found $outputTable{$realLoci + $del_count}[0] instead. (line $j of indel input file $indelFilename). Aborting...\n");
								exit 1;
							}
							elsif ($outputTable{$realLoci + $del_count}[1] ne "not_modified"){
								print_all("*** WARNING: Encountered indel loci " . ($realLoci + $del_count) . ", which is already modified in base ref sequence! (line $j of indel input file $indelFilename) Skipping...\n");
								$warnings++;
							}
							else{
								print_debug("Deleting " . $outputTable{$realLoci + $del_count}[0]  . " at " . ($realLoci + $del_count) . "\n");
								$outputTable{$realLoci + $del_count} = ["","modified"];
							}
							$del_count++;
						}
					}
					else{
						print_all("*** ERROR: Indel file line format not recognized. Got \"$_\" at line $j. Can't figure if it is deletion or insertion (substitutions should be defined in the snp file).\n");
						print "Failed.\n";
						exit 1;
					}
				}
				else{
					print_all("*** ERROR: Bad NUMCER indel line format! $_\n");
					exit 1;
				}

			}
			else{
				print_all("Cannot determine type of indel line: $_ (in $indelFilename, line $j)\n");
				exit 1;
			}
		}
		else{
			print_all("Cannot parse indel line: $_ (in $indelFilename, line $j)\n");
			exit 1;
		}
	}
	close($indelfile);
}

# Dump out the modified base ref sequence in FASTA format for the final strain.
$charCount = 0;
print $outfile ">$currentStrain\n";
foreach $loci (sort {$a <=> $b} keys %outputTable){

  $charCount += length($outputTable{$loci}[0]);

  # Limit to 70 chars per line, for FASTA compliance.
  if ($charCount > 70){
    print $outfile "\n";
    $charCount = 0;
  }
  print $outfile "$outputTable{$loci}[0]";
}
print $outfile "\n";

print_all("\n$warnings warnings were detected.\n");
print "Log file for this run can be found in $ARGV[1].swapped.log \n";

close($snpfile);
close($outfile);

#############
# FUNCTIONS
#############
sub print_debug(){
  # Prints output to STDOUT and log file if debug flag is set.  Otherwise nothing.
  # If the quiet function is also set, then only log to file.
  #
  # Parameters:
  # $1 = Text to print.
  #
  if ($debug eq "Y"){
    if ($quiet eq "N"){
      print "$_[0]";
    }
    print $logfile "$_[0]";
  }
}

sub print_all(){
  # Silly MUX'ed function to write output to both standard out and logfile in one call.
  # If the quiet function is set, then only log to file.
  #
  # Parameters:
  # $1 = Text to print.
  #
  if ($quiet eq "N"){
    print "$_[0]";
  }
  print $logfile "$_[0]";
}

