#!/usr/bin/perl
use warnings;
use strict;

#This tool pulls the index value out of the first line of each read of a FASTQ file.
#The values are stored and then reported.
#This tool is especially useful for evalulating "Undetermined Indices" FASTQ files 
#after the completion of the demultiplexing process.

#This program takes one argument, which is mandatory. The argument specifies a
#"pre_align" path wherein Undetermined FASTQ files might reside.

my $prealignArgument = $ARGV[0];
#Add a forward slash to the end of the basecall path for subsequent operation
unless ($prealignArgument =~ /\/$/){
	$prealignArgument .= "/";
	print "Added a slash to the end of the filepath.\n";
}

#Store the Illumina indexes
my %illuminaIndexes;
foreach ("ATTCCT", "ACTGAT"," GAGTGG", "GTGGCC", "GTGAAA", "GTCCGC", "CCGTCC", "ATGTCA", "AGTTCC", "AGTCAA", "CTTGTA", "GGCTAC" , "TAGCTT", "GATCAG", "ACTTGA", "CAGATC", "GCCAAT", "ACAGTG", "TGACCA", "TTAGGC", "CGATGT", "ATCACG"){
	$illuminaIndexes{$_} = 1;
}

chdir $prealignArgument or die "Can't chdir to $prealignArgument. Shucks...\n";
my @pathContents = <*>; #Store the directory contents in an array for subsequent analysis
my $validPath = 0;
my $undetPath; #holder for the undetermined path

foreach my $testItem (@pathContents) {
	print "Examining $testItem...\n";
	if ($testItem =~ /^Undetermined_indices/){
		$validPath = 1;
		$undetPath = "${prealignArgument}${testItem}/"; #Construct a path to the undetermined-index FASTQ files.
		print "Found the undetermined indices path: $undetPath\n";
		last;
	}
}

unless($validPath){
	die "I can't find a path containing FASTQ files with undetermined indices!\n";
} else {
	print "Valid directory found, let's roll.\n";
}

chdir $undetPath or die "Can't reach the path for the Undetermined Indices FASTQ Files :(\n";
my @samplePathsTest = <*>;
my %sampleDirs;

foreach my $sampleTest (@samplePathsTest){
	if($sampleTest =~ /^Sample/i){
		$sampleDirs{$sampleTest} = "${undetPath}${sampleTest}/"; #Add a sample path to the sampleDirs hash
		print "A new sample path is: $sampleDirs{$sampleTest}\n";
	}
}

foreach my $sampleDir (keys %sampleDirs){
	chdir $sampleDirs{$sampleDir} or die "Can't reach the following path:\n$sampleDir\n";
	my @currentPathFiles = <*>;
	foreach my $testFile (@currentPathFiles) {
		my $unzipFileName;
		if ($testFile =~ /^(lane.*)\.gz$/){	#file is gzipped
			$unzipFileName = $1;
			print "The file $unzipFileName is zipped. Let's unzip it first...\n";
			print "/usr/bin/pigz -d -p 8 $testFile\n";
			print `/usr/bin/pigz -d -p 8 $testFile`;
			print "File is unzipped. Proceeding with the analysis...\n";
			&computeStats($unzipFileName);
			print "Rezipping to leave the file as we found it...\n";
			print "/usr/bin/pigz -p 8 $unzipFileName\n";
			print `/usr/bin/pigz -p 8 $unzipFileName`;
		} elsif ($testFile =~ /^(lane.*)/){ #FASTQ file is not compressed
			$unzipFileName = $1;
			print "It is safe to proceed with the analysis of $unzipFileName\n";
			&computeStats($unzipFileName);
		} else { #we're not looking at a FASTQ and so we'll skip it
			print "Skipping $testFile\n";
		}
	}
}

sub computeStats(){
	my $computeFile = shift(@_); #Only one argument should be passed to this function, and it should be a filename
	my $statsFileName; #placeholder variable for the stats file name
	my %indexes; #Reset the indexes hash so the totals are cumulative
	if($computeFile =~ /^(lane[\d]{1}).*R([12]{1})/){ #capture the lane and read numbers
		$statsFileName = "stats_$1_read$2.txt";
		print "Stats file name will be: $statsFileName\n";
	}
	else{
		die("Regex check of $computeFile failed! \nIt looks like the compute function was passed a bogus file!\n");
	}

	open FASTQFILE, "<$computeFile" or die "Can't open the specified fastq file!\n";
	print "File opened. Writing to RAM...\n";
	my @fileLines = <FASTQFILE>;
	my $numLines = @fileLines;
	my $lineOfFile = 0;
	close FASTQFILE;

	#Validate the FASTQ file
	unless ($numLines % 4 == 0){
		die "Something's up with the FASTQ file. It doesn't have a total line count divsible by four. Aborting...\n";
	} else {
		print "File passed validation\n";
	}
	
	#Iterate through the file and pull the index from each read
	while ($lineOfFile < $numLines){
		if($lineOfFile % 4 == 0){
			if ($fileLines[$lineOfFile] =~ /:([ACTGN]+)$/){
				my $index = $1;
				unless ($indexes{$index}){
					$indexes{$index} = 1;
					print "A new index is: $index\n";
				}
				my $tempCount = $indexes{$index};
				$tempCount++;
				$indexes{$index} = $tempCount;
			}
			$lineOfFile++;
		} else {
			$lineOfFile++;
		}
	}
	
	#Now to create the stats file
	open STATSFILE, ">$statsFileName" or die "Can't create file handle!\n$!\n";
	my $numFailedIndexReads = 0;
	my $numSuccessfulIndexReads = 0;
	my $numFailedIndexes = 0;
	my $numSuccessfulIndexes = 0;
	foreach my $anIndex (sort {$indexes{$b} <=> $indexes{$a}} keys %indexes){
		print STATSFILE "The count for $anIndex is: $indexes{$anIndex}";
		if (exists $illuminaIndexes{$anIndex}){
			print STATSFILE "   <--This is also the sequence of an Illumina Index\n";
		} else {
			print STATSFILE "\n";
		}
		if (length($anIndex) != 6){ #An index must be exactly six characters or it's a "failed read."
			$numFailedIndexReads += $indexes{$anIndex};
			$numFailedIndexes++;
		} elsif ($anIndex =~ /N{2,}/) { #An index must not contain more than one "N" read.
			$numFailedIndexReads += $indexes{$anIndex};
			$numFailedIndexes++;
		} else {
			$numSuccessfulIndexReads += $indexes{$anIndex}; #The index was a "successful read."
			$numSuccessfulIndexes++;
		}
	}
	
	my $numTotalReads = $numSuccessfulIndexReads + $numFailedIndexReads;
	my $percentFailed = sprintf("%.2f%%", ($numFailedIndexReads/$numTotalReads) * 100);
	my $percentSuccessful = sprintf("%.2f%%", ($numSuccessfulIndexReads/$numTotalReads) * 100);
	
	print STATSFILE "Number of total unique successful indexes: $numSuccessfulIndexes\n";
	print STATSFILE "Number of total unique failed indexes: $numFailedIndexes\n";
	print STATSFILE "Number of total successful index reads: $numSuccessfulIndexReads\n";
	print STATSFILE "Percentage of total reads: $percentSuccessful\n";
	print STATSFILE "Number of total failed index reads: $numFailedIndexReads\n";
	print STATSFILE "Percentage of total reads: $percentFailed\n";
	close STATSFILE;
}

sub usageExit(){
	die ("fastqIndexCounter.pl\nUsage: fastqIndexCounter.pl <pathToCASAVAPrealign> A tool for counting which indexes appeared in the \"Undetermined_indices\" FASTQ files for a flowcell output by CASAVA.\n Each lane in a flowcell can have its own Undetermined_indices FASTQ file. fastqIndexCounter.pl will unzip each of these files, count the indexes, and create a stats textfile.\nIf you have questions, contact Dan Dorset at Vanderbilt VANTAGE: daniel.dorset@vanderbilt.edu")
}