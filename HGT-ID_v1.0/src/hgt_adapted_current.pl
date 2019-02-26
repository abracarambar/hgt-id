#!/usr/bin/perl

#### Required bedtools & samtools to be in path
use Getopt::Long;
use strict;
use File::Basename;
use FindBin qw($Bin);
use Cwd;
use Cwd qw(chdir);
use File::Copy;
use File::Path qw(make_path remove_tree rmtree);
my ($BAMFILE,$SAMPLE,$OUTPUT,$CONFIG,$verbose,$debug);

#Declare variables
GetOptions(
	'b=s' => \$BAMFILE,
	's=s' => \$SAMPLE,
	'o=s' => \$OUTPUT,
	'c=s' => \$CONFIG,
	'v' => \$verbose,
	'd' => \$debug,
	"help|h|?"	=> \&usage);

print "Usage = hgt.pl -b $BAMFILE -c $CONFIG \n\n";
sub usage {
	print "\nusage: hgt.pl [-sov] -b <BAMFILE> -c <CONFIG> \n";
	print "\t-s\t\tsample Name [Extracted from BAM header]\n";
	print "\t-o\t\tOutput Directory [cwd]\n";
	print "\t-v\t\tverbose flag [no]\n";
	print "\t-d\t\tdebug flag [no]\n";
	exit 1;
	}
if(defined($BAMFILE)){$BAMFILE=$BAMFILE} else {print usage();die "Where is the BAM file?\n\n"}
if(defined($CONFIG)){$CONFIG=$CONFIG} else {print usage();die "Where is the configuration file?\n\n"}

##remove it
#rmtree("$OUTPUT");
#if(defined($OUTPUT)){$OUTPUT=$OUTPUT;if (-d "$OUTPUT" || -e "$OUTPUT"){die "Output folder already exist\n";} elsif ( ! -d "$OUTPUT") {mkdir $OUTPUT}} else {$OUTPUT=getcwd;}
if(defined($OUTPUT)){$OUTPUT=$OUTPUT; if ( ! -d "$OUTPUT") {mkdir $OUTPUT}} else {$OUTPUT=getcwd;}
chdir ($OUTPUT) or die "cannot change: $!\n";

### creating the directory structure
my $logs=$OUTPUT ."/.logs";
my $autocode=$OUTPUT ."/.scripts";
mkdir $logs;
mkdir $autocode;
my $human_mapping=$OUTPUT ."/.human";
my $human_mapping_again=$OUTPUT ."/.human_again";
my $viral_mapping=$OUTPUT ."/.virus";
my $scoring=$OUTPUT ."/.scoring";
mkdir $human_mapping;
mkdir $human_mapping_again;
mkdir $viral_mapping;
mkdir $scoring;

### reading the configuration file
copy($CONFIG,"$OUTPUT/config.txt");
$CONFIG="$OUTPUT/config.txt";
my $config_vars=read_files_var($CONFIG);
my $SAMTOOLS=$config_vars->{'SAMTOOLS'};
my $BEDTOOLS=$config_vars->{'BEDTOOLS'};
my $CYTOBAND=$config_vars->{'HUMAN_CYTOBAND'};
my $BWA=$config_vars->{'BWA'};
my $PICARD=$config_vars->{'PICARD'};
my $FANCYBOX=$config_vars->{'FANCYBOX'};
my $RLIB=$config_vars->{'RLIB'};
### add tools to the path
$ENV{'PATH'} = $SAMTOOLS . ':' . $BEDTOOLS . ':' . $BWA . ':' . $PICARD . ':' . $ENV{'PATH'};

### get the path to all the scripts
my $SCRIPT_DIR=$Bin;

### references
my $USER_HUMAN_database=$config_vars->{'USER_HUMAN_database'};
my $HUMAN_database=$config_vars->{'HUMAN_database'};
my $HUMAN_database_Index=$config_vars->{'HUMAN_database_Index'};
my $VIRUS_database=$config_vars->{'VIRUS_database'};
my $VIRUS_database_Index=$config_vars->{'VIRUS_database_Index'};
my $VIRUS_HUMAN_database=$config_vars->{'VIRUS_HUMAN_database'};
my $VIRUS_HUMAN_database_Index=$config_vars->{'VIRUS_HUMAN_database_Index'};
my $REF_FLAT=$config_vars->{'REF_FLAT'};
my $MINRP=$config_vars->{'MINRP'};
my $MINSOFT=$config_vars->{'MINSOFT'};
### parameters
my $THREADS=$config_vars->{'THREADS'};
my $SEQUENCE_COMPLEXITY=$config_vars->{'SEQUENCE_COMPLEXITY'};
my $DEPTH=$config_vars->{'DEPTH'};

### make sure the index file is available for the BAM file
my ($fn,$pathname) = fileparse($BAMFILE,".bam");
my $index=`ls $pathname/$fn*bai|head -1`;
if(!$index){die "\n\nERROR: you need to index your BAM file using samtools\n\n"}

### get current time
print "Start Time : " . &spGetCurDateTime() . "\n";
my $now = time;
my $command = "";
my $OUTNAME ="";

#skip for now
#extract sample name
#select row name with RG info, loop through columns and find SM tag, extract the sample name, only retain the first line

#if ( ! $SAMPLE)	{
#	$SAMPLE=`samtools view -H $BAMFILE |awk '{if(\$1~/^\@RG/){for(i;i<=NF;i++){if (\$i ~ /^SM:/){col=i}};sub("SM:","",\$col);print \$col}}' | head -1`;
#	$SAMPLE=~s/\n//g;
#}

if($SAMPLE ne ""){$OUTNAME=$OUTPUT . "/" . $SAMPLE.".txt";}
else {$OUTNAME=$OUTNAME=$OUTPUT . "/output.txt";$SAMPLE="dummy";}
print "Writing results to $OUTNAME\n";

#this holds the job 
#extract Inferred insert SIZE
#retain f2 flag
#retain q20 quality
#filter cigar column for Skipped region from the reference
#retain only pairs of sequences that are on the same reference: Mate Reference sequence NaMe (= if same as RNAME)

#my $libsize=`samtools view -q 20 -f2 $BAMFILE | awk '\$6 !~ /N/' | awk '\$7 ~ /=/' | cut -f9| awk '\$1<1000' | head -10000 | awk '{if (\$1<0){\$1=-\$1}else{\$1=\$1} sum+=\$1;} END {printf(\"\%3.0f\\n",sum/NR)}'`;
#chomp $libsize;
#if ( $libsize =~ m/^\D/)        {
#	$libsize=450;
#}
my $libsize="309";
print "library size is inferred as : $libsize\n";

#derive read length
#my $readlength=`samtools view $BAMFILE | head -1 | awk '{print length(\$10)}'`;
#chomp $readlength;

my $readlength="101";
print "$readlength\n";


## to deal with TCGA sample/not mayo samples
#my $chrflag=`samtools view -H $BAMFILE | grep '^\@SQ' | head -n1 | awk '{if(\$2 ~ /^SN:chr/) {print \"yes\"} else {print \"no\"}}'`;
#chomp $chrflag;
#print "$chrflag\n";

print "preprocessing done and calculations completed!\n";

print "extracting unmapped, split and soft-clipped reads from human aligned bam...\n";

&extractcandidatesplitreadsfromhumanbam($BAMFILE, $human_mapping, $scoring, $THREADS, $PICARD, $logs);

&extractcandidateunmappedreadsfromhumanbam($BAMFILE, $human_mapping, $scoring, $THREADS, $PICARD, $logs);

&extractsoftclippedreadsfromhumanbam($BAMFILE, $human_mapping, $human_mapping_again, $THREADS, $PICARD, $SCRIPT_DIR, $HUMAN_database_Index, $logs);

&mapagain($BAMFILE, $human_mapping, $human_mapping_again, $THREADS, $PICARD, $logs);

#&mergeallbamsandextractreads($BAMFILE, $human_mapping, $human_mapping_again, $THREADS, $PICARD, $logs);

&mappingtoviralgenomes($THREADS, $SAMPLE, $VIRUS_database_Index, $human_mapping_again, $viral_mapping, $scoring, $SCRIPT_DIR);

&findHGTcandidates($SCRIPT_DIR, $human_mapping_again, $viral_mapping, $OUTPUT, $logs, $VIRUS_HUMAN_database, $THREADS, $SEQUENCE_COMPLEXITY, $PICARD, $DEPTH);

&extract_regions($readlength, $libsize, $OUTPUT, $VIRUS_HUMAN_database,$HUMAN_database, $BAMFILE, $human_mapping_again, $viral_mapping, $THREADS, $CYTOBAND);

print  "find the integration point...\n";
$command=join ("","perl ",$SCRIPT_DIR,"/integration.pl -b ",$OUTPUT,"/",$SAMPLE,".forcalling.cyto.bam -f ",$VIRUS_HUMAN_database, " -o ",$OUTNAME, " -d ",$libsize, " -m ",$MINRP," -l ",$MINSOFT," -r ",$HUMAN_database," -t -v ");
submit($command);

#### filter the file to remove the reads very close to 5KB centromere using the cytoband file
$command=join("","cat ",$OUTNAME, " | awk '{if(NR==1){print \"#chr\tstart\tend\t\"\$0} else {print \$1\"\\t\"\$2\"\\t\"\$2\"\\t\"\$0}}' > ",$OUTNAME,".bed");
submit($command);
$command=join ("","zcat ",$CYTOBAND," | grep acen | awk '{if(\$4 ~ /^p/) {print \$1\"\\t\"\$2-50000\"\\t\"\$3} else {print \$1\"\\t\"\$2\"\\t\"\$3+50000}}' | intersectBed -a ",$OUTNAME,".bed -b stdin -v -header | cut -f4- > ",$OUTNAME);
`intersectBed -a $OUTNAME.bed -b $CYTOBAND -v -header | cut -f4- > $OUTNAME`;
submit($command);

#### map it to gene
print "create a gene bed file to annotate HGT...\n";
$command=join ("","zcat ",$REF_FLAT," | grep -v random | grep -v chrUn | grep -v hap | awk '{print \$3\"\\t\"\$5\"\\t\"\$6\"\\t\"\$1}' | sortBed -i stdin | perl ",$SCRIPT_DIR,"/uniqgene.pl > ",$OUTPUT,"/gene.bed");
submit($command);
$command=join ("","cat ",$OUTNAME," | awk 'NR>1{print \$1\"\\t\"\$2\"\\t\"\$2}' | closestBed -t first -a stdin -d -b ",$OUTPUT,"/gene.bed  | awk '{print \$1\"\\t\"\$2\"\\t\"\$(NF-1)\"\\t\"\$NF}' > ",$OUTNAME,".gene.txt");
submit($command);
$command=join ("","perl ",$SCRIPT_DIR,"/map.hgt.annot.pl ",$OUTNAME,".gene.txt ",$OUTNAME," | grep -v NA > ",$OUTNAME,".annot.txt");
submit($command);
$command=join ("","mv ",$OUTNAME,".annot.txt ",$OUTNAME);
submit($command);

#### filter the results if it is both the human pairs
open VIRUSDB, "$VIRUS_database.fai" or die "can't open the VIRUS index file\n";
my %virus;
while(my $l=<VIRUSDB>)	{
	chomp $l;
	my @virus_contig=split(/\t/,$l);
	$virus{$virus_contig[0]}=1;
}
close VIRUSDB;
open OUT, "$OUTNAME" or die "can't open the output file\n";
open FLT, ">$OUTNAME.tmp.txt" or die "can't open the output folder\n";
while(my $l = <OUT>)	{
	print FLT $l if($.==1);
	chomp $l;
	my @call=split(/\t/,$l);
	if(defined $virus{$call[2]} || defined $virus{$call[7]})	{
		print FLT "$l\n";
	}
}
close OUT;
close FLT;
$command=join ("","mv ",$OUTNAME,".tmp.txt ",$OUTNAME);
submit($command);


&extract_regions($readlength, $libsize, $OUTPUT, $VIRUS_HUMAN_database,$HUMAN_database, $BAMFILE, $human_mapping_again, $viral_mapping, $THREADS, $CYTOBAND);




#fix !!!!
#&getviralcoverage($coverage, $OUTNAME, $OUTPUT, $autocode, $VIRUS_database.fai, $SCRIPT_DIR, $viral_mapping, $scoring);

### primer design
#$command=join("","sh ",$SCRIPT_DIR,"/PrimerDesign_HGT_ID.sh ",$OUTNAME," ",$OUTPUT," ",$CONFIG);
#submit($command);

#### create CIRCOS plots
#my $plot=$OUTPUT ."/circosplot";
#mkdir $plot;
#$command=join("","Rscript ",$SCRIPT_DIR,"/circos.Rscript ",$OUTNAME," ",$OUTPUT," ",$plot," ",$RLIB);
#submit($command);

#### get the scoring numbers
#get the Human Softclipped Reads

open FH, ">$autocode/scores.sh" or die "can't open the script to write\n";
print FH "echo \"HumanSoftClipping\" > $scoring/HSoft.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f5-7 | awk '{print \$1\":\"\$2-$readlength\"-\"\$3+$readlength}'`;\ndo\nsamtools view -f2 $OUTPUT/$SAMPLE.forcalling.cyto.bam \$i | awk '\$6 ~/S/' | awk '\$6 !~ /D/'  | awk '\$6 !~ /I/' | perl -ane 'my \@CIGAR = split(/([0-9]+[SMIDNHXP])/, \$F[5]); my \$hash; map { push(\@{\$hash->{\$2}}, \$1) if (/(\\d+)([SMIDNHXP])/) } \@CIGAR; \$count=0;foreach my \$softclip (\@{\$hash->{S}}) {if(\$softclip >= $MINSOFT){\$count++}}; print \"\$count\\n\"' | awk '\$1>0' | wc -l >> $scoring/HSoft.txt;\ndone\n";
print FH "echo \"ViralProperMappedReads\" > $scoring/VProper.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f8-10 | awk '{print \$1\":\"\$2\"-\"\$3}'`;\ndo\nsamtools view -f2 $scoring/virus.fromUnmapped.sort.bam \$i | cut -f1 | sort | uniq -c | wc -l >> $scoring/VProper.txt;\ndone\n";
print FH "echo \"ViralSoftClipping\" > $scoring/VSoft.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f8-10 | awk '{print \$1\":\"\$2\"-\"\$3}'`;\ndo\nsamtools view -f2 $scoring/virus.fromUnmapped.sort.bam \$i | awk '\$6 ~/S/' | awk '\$6 !~ /D/'  | awk '\$6 !~ /I/' | perl -ane 'my \@CIGAR = split(/([0-9]+[SMIDNHXP])/, \$F[5]); my \$hash; map { push(\@{\$hash->{\$2}}, \$1) if (/(\\d+)([SMIDNHXP])/) } \@CIGAR; \$count=0;foreach my \$softclip (\@{\$hash->{S}}) {if(\$softclip >= $MINSOFT){\$count++}}; print \"\$count\\n\"' | awk '\$1>0' | wc -l >> $scoring/VSoft.txt;\ndone\n";
print FH "echo \"HumanProperMappedReads\" > $scoring/HProper.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f5-7 | awk '{print \$1\":\"\$2-$readlength\"-\"\$3+$readlength}'`; do samtools view -f2 $OUTPUT/$SAMPLE.forcalling.cyto.bam \$i |cut -f1 | sort | uniq -c | wc -l  >> $scoring/HProper.txt;\ndone\n";
#print FH "if [ $chrflag == \"yes\" ]\n";
print FH "then\n";
print FH "echo \"TotalCoverage\" > $scoring/TotalCoverageH.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f5-7 | awk '{print \$1\":\"\$2-$libsize\"-\"\$3+$libsize}'`; do samtools view $BAMFILE \$i | wc -l >> $scoring/TotalCoverageH.txt;\ndone\n";
print FH "else\n";
print FH "echo \"TotalCoverage\" > $scoring/TotalCoverageH.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f5-7 | awk '{print \$1\":\"\$2-$libsize\"-\"\$3+$libsize}' | sed -e 's/chr//g'`; do samtools view $BAMFILE \$i | wc -l >> $scoring/TotalCoverageH.txt;\ndone\n";
print FH "fi\n";
print FH "echo \"TotalCoverageViral\" > $scoring/TotalCoverageV.txt;\nfor i in `cat $OUTNAME | awk 'NR>1' | cut -f8-10 | awk '{print \$1\":\"\$2\"-\"\$3}'`; do samtools view $scoring/virus.fromUnmapped.sort.bam \$i | wc -l >> $scoring/TotalCoverageV.txt;\ndone\n";
print FH "paste $scoring/TotalCoverageV.txt $scoring/TotalCoverageH.txt | awk 'NR>1' | awk 'BEGIN{print \"Coverage\"} {print \$1+\$2}' > $scoring/TotalCoverage.txt\n";

print FH "cat $OUTNAME  | cut -f1-11 | paste - $scoring/HSoft.txt $scoring/HProper.txt $scoring/VProper.txt $scoring/VSoft.txt $scoring/TotalCoverage.txt  | awk '{if(NR==1){print \$0\"\\tScore\"} else {print \$0\"\\t\"(\$11+\$12+\$15-(((\$13+\$12)\*(\$14+\$15))/\$16))}}' > $OUTNAME.tmp.txt\n";
print FH "mv $OUTNAME.tmp.txt $OUTNAME\n";
$command=join("","chmod 777 ",$autocode,"/scores.sh");
submit($command);
$command=join ("","sh ",$autocode,"/scores.sh");
submit($command);

#### create an HTML page for the report and circos plot
#my $fancybox=$OUTPUT ."/fancybox";
#mkdir $fancybox;
#`cp -Rf $FANCYBOX/source $fancybox`;
#$command=join("","perl ",$SCRIPT_DIR,"/create_html.pl ",$OUTNAME," ",$OUTPUT,"/output.primers ",$OUTPUT," > ",$OUTPUT,"/results.html");
#submit($command);

#print "Finish Time : " . &spGetCurDateTime() . "\n";
#$now = time - $now;
#printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60),
#int($now % 60));

exit;
#########################
#### completed all the steps
############################################################

#######################
### subroutines
#######################


### parse the configuration file
sub read_files_var{
	my ($filename) = @_;
		open INFILE, "<$filename" or die "$filename $!\n";
	my %variables;
	LINE:while (<INFILE>) {
	next LINE if /^#|^\s/;
	chomp;

	my ($var, $value) = split (/=/,$_);
	$value =~ tr/\"//d;
	$value =~ s/\s+$//;
	$variables{$var}=$value;
	}
	return (\%variables);
}
### get time
sub spGetCurDateTime {
	my ($sec, $min, $hour, $mday, $mon, $year) = localtime();
	my $curDateTime = sprintf "%4d-%02d-%02d %02d:%02d:%02d",
	$year+1900, $mon+1, $mday, $hour, $min, $sec;
	return ($curDateTime);
}

sub submit{
	#if you want to sort the output
	$command=shift;
	if($verbose){	print "$command\n";}
	system("$command");
}

sub extract_regions {
	my ($readlength, $libsize, $OUTPUT, $VIRUS_HUMAN_database,$HUMAN_database, $BAMFILE, $human_mapping_again, $viral_mapping, $THREADS, $CYTOBAND) = @_;
	print "#### creating regions to go back to original BAM file to extract reads to find soft clipping...\n";
	my $region2Capture=$readlength+$libsize;
	$command=join ("","bamToBed -i ",$OUTPUT,"/VIRUS.HUMAN.flt.sort.bam | sortBed -i stdin | mergeBed -i stdin | slopBed -i stdin -g ",$VIRUS_HUMAN_database,".fai -l ",$region2Capture," -r ",$region2Capture," | mergeBed -i stdin > ",$OUTPUT,"/regions.bed");
	submit($command);
	my $chrs2keep=`cat $OUTPUT/chromosomes2keep.txt | tr \"\\n\" \"|\" | sed -e 's/\\(\.\*\\)./\\1/'`;
	$command=join("","cat ",$OUTPUT,"/regions.bed | grep -w -E '",$chrs2keep,"' > ",$OUTPUT,"/regions.2keep.bed");
	submit($command);

	if (-z "$OUTPUT/regions.2keep.bed"){
		print "***********************************\n";
		print "Bad Luck !! NO HGT candidates found\n";
		print "***********************************\n";
		open OUT, ">$OUTNAME" or die "failed to open the file\n";
		print OUT "HumanChr\tHIntegrationPos\tHChr\tHStart\tHPosition\tVirusID\tVStart\tVirusEnd\tDiscordantReads\tSoftclippingReads\n";
		close OUT;
		exit 0;
	}
	my $whichchr="";
	my $chr2look=`cat $HUMAN_database.fai | cut -f1 |  tr \"\\n\" \"|\" | sed -e 's/\\(\.\*\\)./\\1/'`;
	$whichchr=`cat $OUTPUT/regions.2keep.bed  | cut -f1 | sort | uniq | grep -w -E '$chr2look' | tr "\n"  " "| sed \'s/\\s*\$\/\/g\'`;
	chomp $whichchr;
	print "go back to original BAM file to extract nearby reads...\n";

	#org for original
	`samtools view -b -F 256 -F 1024 -L $OUTPUT/regions.2keep.bed $BAMFILE $whichchr > "$OUTPUT/human.org.bam`;
	$command=join ("","samtools view -h ",$OUTPUT,"/human.org.bam  | sed -e 's/SN:chr\\([0-9XYG]\\)/SN:\\1/' -e 's/SN:chrMT/SN:MT/' -e 's/SN:chrM_rCRS/SN:MT/' | samtools reheader - ",$OUTPUT,"/human.org.bam > ",$OUTPUT,"/human.org.tmp.bam");
	submit($command);
	$command=join ("","mv ",$OUTPUT,"/human.org.tmp.bam ",$OUTPUT,"/human.org.bam");
	submit($command);
	#}
	#2 = proper pair
	#filter 256 not primary alignment
	#filter 1024 read is PCR or optical duplicate
	`samtools view -b -f2 -F 256 -F 1024 -L $OUTPUT/regions.2keep.bed $human_mapping_again/2_human.bam > $OUTPUT/human.again.bam`;
	`samtools view -b -f2 -F 256 -F 1024 -L $OUTPUT/regions.2keep.bed $viral_mapping/virus.bam > $OUTPUT/virus.org.bam`;

	print "#### merge the BAM file to find the integration point...\n";
	`samtools merge -f $OUTPUT/merged.bam $OUTPUT/VIRUS.HUMAN.flt.sort.bam $OUTPUT/human.org.bam $OUTPUT/human.again.bam $OUTPUT/virus.org.bam`;
	`samtools sort -@ $THREADS $OUTPUT/merged.bam $OUTPUT/$SAMPLE.forcalling`;
	### filter the BAM file to remove the reads very close to 5KB centromere using the cytoband file
	###small centromere start + 5000 = end <=> long centromere start  = end + 5000
	#$command=join ("","zcat ",$CYTOBAND," | grep acen | awk '{if(\$4 ~ /^p/) {print \$1\"\\t\"\$2-50000\"\\t\"\$3} else {print \$1\"\\t\"\$2\"\\t\"\$3+50000}}' | intersectBed -abam ",$OUTPUT,"/",$SAMPLE,".forcalling.bam -b stdin -v > ",$OUTPUT,"/",$SAMPLE,".forcalling.cyto.bam" );
	`intersectBed -abam $OUTPUT/$SAMPLE.forcalling.bam -b $CYTOBAND -v > $OUTPUT/$SAMPLE.forcalling.cyto.bam`;
	#$command=join ("","mv ",$OUTPUT,"/",$SAMPLE,".forcalling.cyto.bam ", $OUTPUT,"/",$SAMPLE,".forcalling.bam");
	`samtools index $OUTPUT/$SAMPLE.forcalling.cyto.bam`;
	print "#### completed extracting regions!\n";
}

 sub mappingtoviralgenomes {
	my ($THREADS, $SAMPLE, $VIRUS_database_Index, $human_mapping_again, $viral_mapping, $scoring, $SCRIPT_DIR) = @_;
	print "#### maping the partially mapping reads to viral genomes...\n";
 	`bwa mem -t $THREADS -M -R \'\@RG\\tID:$SAMPLE\\tSM:$SAMPLE\' $VIRUS_database_Index $human_mapping_again/read1.fq $human_mapping_again/read2.fq | samtools view -u -bS - >  $viral_mapping/virus.bam`;
	print "#### completed mapping to viral genome!\n";
	
	print "#### filter split reads into oneEndMapped.bam and oneEndUnMapped.bam...\n"; 
    `samtools view -u -b -f 8 -F 260 $viral_mapping/virus.bam > $viral_mapping/oneEndMapped.bam`;
    `samtools view -u -b -f 4 -F 264 $viral_mapping/virus.bam > $viral_mapping/oneEndUnMapped.bam`;
    `samtools merge -f -u -n $viral_mapping/2_virus.bam $viral_mapping/oneEndMapped.bam $viral_mapping/oneEndUnMapped.bam`;
    `samtools sort -@ $THREADS -n $viral_mapping/2_virus.bam $viral_mapping/virus.sort`;
    #this removes any white space characters \s
    `perl $SCRIPT_DIR/keepreads.pl $viral_mapping/virus.sort.bam | sed '/^\\s*\$/d' | samtools view -u -bS - > $viral_mapping/virus.fix.sort.bam`;
    
 	print "#### mapping the unmapped reads to viral genomes...\n";
   	`bwa mem -t $THREADS -M -R \'\@RG\\tID:$SAMPLE\\tSM:$SAMPLE\' $VIRUS_database_Index $scoring/read1.fq $scoring/read2.fq | samtools view -u -bS - > $scoring/virus.fromUnmapped.bam`;    
    `samtools sort -@ $THREADS $scoring/virus.fromUnmapped.bam $scoring/virus.fromUnmapped.sort`;
    `samtools index $scoring/virus.fromUnmapped.sort.bam`;  
}

sub findHGTcandidates {
	my ($SCRIPT_DIR, $human_mapping_again, $viral_mapping, $OUTPUT, $logs, $VIRUS_HUMAN_database, $THREADS, $SEQUENCE_COMPLEXITY, $PICARD, $DEPTH) = @_;
	print "#### finding the initial candidate list of HGT...\n";
	#### create report for virus and human mapping for each HGT candidate
	`perl $SCRIPT_DIR/map.virus.human.pl $human_mapping_again/human.fix.sort.bam $viral_mapping/virus.fix.sort.bam $OUTPUT/HGT.candidates.txt $OUTPUT/VIRUS.HUMAN.sam > $logs/5.map.virus.human.log`;
	#delete all blank lines from a file
	`cat $OUTPUT/VIRUS.HUMAN.sam | sed '/^\$/d' | samtools view -bt $VIRUS_HUMAN_database.fai - > $OUTPUT/VIRUS.HUMAN.bam`;
	#`cat $OUTPUT/VIRUS.HUMAN.sam | sed '/^$/d' | samtools view -bt $VIRUS_HUMAN_database_Index - > $OUTPUT/VIRUS.HUMAN.bam`;
	`samtools sort -@ $THREADS $OUTPUT/VIRUS.HUMAN.bam $OUTPUT/VIRUS.HUMAN.sort`;
	`samtools index $OUTPUT/VIRUS.HUMAN.sort.bam`;
	print "filtering candidate list using sequence complexity fiter ...\n";
	`cat $OUTPUT/HGT.candidates.txt | perl $SCRIPT_DIR/filter.pl  | awk 'NR>1' | awk -v complex=$SEQUENCE_COMPLEXITY '\$NF>=complex' > $OUTPUT/HGT.candidates.flt.txt`;
	`java -Xmx6g -Xms3g -jar $PICARD/FilterSamReads.jar I=$OUTPUT/VIRUS.HUMAN.sort.bam RLF=$OUTPUT/HGT.candidates.flt.txt FILTER=includeReadList O=$OUTPUT/VIRUS.HUMAN.flt.sort.bam WRITE_READS_FILES=false VALIDATION_STRINGENCY=SILENT TMP_DIR=$OUTPUT/tmp_dir CREATE_INDEX=TRUE > $logs/6.FilterSamReads.log`;
	`samtools idxstats $OUTPUT/VIRUS.HUMAN.flt.sort.bam | awk -v depth=$DEPTH '\$NF+\$(NF-1)>depth' | cut -f1 > $OUTPUT/chromosomes2keep.txt`;
	print "#### completed finding the initial candidate list of HGT!\n";
}

sub extractcandidatesplitreadsfromhumanbam {
	my ($BAMFILE, $human_mapping, $scoring, $THREADS, $PICARD, $logs) = @_;
	print "#### filter split reads into oneEndMapped.bam and oneEndUnMapped.bam...\n";
	print "#### extracting reads from human aligned bam for first pass mapping back to human...\n";
	print $BAMFILE;
	# R1 mapped, R2 unmapped
	`samtools view -u -b -f 8 -F 260 $BAMFILE > $human_mapping/oneEndMapped.bam`;
	# R1 unmapped, R2 mapped
	`samtools view -u -b -f 4 -F 264 $BAMFILE > $human_mapping/oneEndUnMapped.bam`;
}

sub extractcandidateunmappedreadsfromhumanbam {
	my ($BAMFILE, $human_mapping, $scoring, $THREADS, $PICARD, $logs) = @_;
	print "#### extracting unmapped reads...\n";
	# R1 & R2 unmapped
	### both reads are unmapped
	`samtools view -u -b -f 12 $BAMFILE > $scoring/UnMapped.bam`;
	#### get the FASTQ from BAM file for unmapped reads
	$command=join ("","java -XX:ParallelGCThreads=",$THREADS," -Xmx6g -Xms3g -jar ",$PICARD,"/SamToFastq.jar INPUT=",$scoring,"/UnMapped.bam FASTQ=",$scoring,"/read1.fq SECOND_END_FASTQ=",$scoring,"/read2.fq VALIDATION_STRINGENCY=SILENT > ",$logs,"/1a.bam2fastq.log 2>&1");
	submit($command);
}

sub extractsoftclippedreadsfromhumanbam {
	my ($BAMFILE, $human_mapping, $human_mapping_again, $THREADS, $PICARD, $SCRIPT_DIR, $HUMAN_database_Index, $logs) = @_;
	#soft-clipped reads where softclip segment is more or equal to half the length of sequencing read
	#perl -e allows to run perl command on the command line
	#perl -n adds loops around the -e code, here processes the filtered sam file that is piped, a line at a time
	#perl -a option turns on autosplit mode. In this mode, each input record is split and the resulting list of elements is stored in an array, here called @F
	`samtools view -F 8 -f 2 -F 4 -f 64 $BAMFILE \\
	 | awk '\$6 ~ /S/' | awk '\$6 !~ /I/ ' | awk '\$6 !~ /D/' | awk '\$6 ~ /S/' \\
	 | perl -ane '\$len=length(\$F[9]);\$len=\$len/2;\@CIGAR = split(/([0-9]+[SMIDNHXP])/, \$F[5]);my \$hash; \\
	map { push(\@{\$hash->{\$2}}, \$1) if (/(\\d+)([SMIDNHXP])/) } \@CIGAR; foreach my \$softclip (\@{\$hash->{S}}) {if(\$softclip >= \$len){ print }}' > $human_mapping/read1`;

	`samtools view -F 8 -f 2 -F 4 -f 128 $BAMFILE \\
	 | awk '\$6 ~ /S/' | awk '\$6 !~ /I/ ' | awk '\$6 !~ /D/' | awk '\$6 ~ /S/' \\
	 | perl -ane '\$len=length(\$F[9]);\$len=\$len/2; \@CIGAR = split(/([0-9]+[SMIDNHXP])/, \$F[5]); \\
	 my \$hash;map { push(\@{\$hash->{\$2}}, \$1) if (/(\\d+)([SMIDNHXP])/) } \@CIGAR; foreach my \$softclip (\@{\$hash->{S}}) {if(\$softclip >= \$len){ print }}'> $human_mapping/read2`;

	#get the read ids for the extracted softclipped IDs
	print "#### get all the soft clipped read IDs\n";
	`cat $human_mapping/read1 $human_mapping/read2  | cut -f1 | sort -T \$PWD | uniq > $human_mapping/IDS`;
	`if [ -s $human_mapping/IDS ];\n then\n java -XX:ParallelGCThreads=$THREADS -Xmx6g -Xms3g -jar $PICARD/FilterSamReads.jar I=$BAMFILE WRITE_READS_FILES=false RLF=$human_mapping/IDS FILTER=includeReadList O=$human_mapping/softclip.bam VALIDATION_STRINGENCY=SILENT SORT_ORDER=queryname TMP_DIR=$human_mapping/tmp_dir  > $logs/0.FilterSamReads.log 2>&1;\n fi`;

	print "#### split reads with soft clips\n";
	#returns error for chromosomes that are not part of their reference
	`if [ -s $human_mapping/softclip.bam ];\n then\n samtools view $human_mapping/softclip.bam | perl $SCRIPT_DIR/splitReads.pl | samtools view -bt /g/data1a/jp48/scripts/hgtid/HGT-ID_v1.0/resources/human.fa.fai - >  $human_mapping/soft.bam;\n fi`;
	
	#added this step, get back to it
	#`samtools sort -\@ $THREADS -n $human_mapping/soft.bam $human_mapping/soft.sorted`;
	
	#merge all bams obatined in a preliminary human bam
	`samtools merge -u -n $human_mapping/human.bam $human_mapping/oneEndMapped.bam $human_mapping/oneEndUnMapped.bam $human_mapping/soft.bam`;
	`samtools sort -\@ $THREADS -n $human_mapping/human.bam $human_mapping/human.sort`;
	`samtools view -h $human_mapping/human.sort.bam | grep -v XF | samtools view -bS - > $human_mapping/tmp.bam`;
	`cp $human_mapping/tmp.bam $human_mapping/human.sort.bam`;
	`perl $SCRIPT_DIR/keepreads.pl $human_mapping/human.sort.bam | sed '/^\\s*\$/d' | samtools view -u -bS - > $human_mapping/human.fix.sort.bam`;
	#sed '/^\\s*\$/d' 
	print "#### converting the BAM file to fastq...\n";
	#$command=join ("","java -XX:ParallelGCThreads=",$THREADS," -Xmx6g -Xms3g -jar ",$PICARD,"/SamToFastq.jar INPUT=",$human_mapping,"/soft.sorted.bam FASTQ=",$human_mapping,"/read1.fq SECOND_END_FASTQ=",$human_mapping,"/read2.fq VALIDATION_STRINGENCY=SILENT > ",$logs,"/1.bam2fastq.log 2>&1");
	#submit($command);
	$command=join ("","java -XX:ParallelGCThreads=",$THREADS," -Xmx6g -Xms3g -jar ",$PICARD,"/SamToFastq.jar INPUT=",$human_mapping,"/soft.fix.sort.bam FASTQ=",$human_mapping,"/read1.fq SECOND_END_FASTQ=",$human_mapping,"/read2.fq VALIDATION_STRINGENCY=SILENT > ",$logs,"/1.bam2fastq.log 2>&1");
	submit($command);
	
	#`bwa mem -t $THREADS -M -R \'\@RG\\tID:$SAMPLE\\tSM:$SAMPLE\' $HUMAN_database_Index $human_mapping/read1.fq $human_mapping/read2.fq | samtools view -u -bS - > \\
 	# $human_mapping_again/human_soft.bam`;

	#print "#### extracting reads from human aligned bam for viral mapping...\n";
	#`samtools view -u -b -f 8 -F 260 $human_mapping_again/human_soft.bam > $human_mapping_again/oneEndMapped_soft.bam`;
	#`samtools view -u -b -f 4 -F 264 $human_mapping_again/human_soft.bam > $human_mapping_again/oneEndUnMapped_soft.bam`;
	#### get the FASTQ from BAM file for unmapped reads
	$command=join ("","java -XX:ParallelGCThreads=",$THREADS," -Xmx6g -Xms3g -jar ",$PICARD,"/SamToFastq.jar INPUT=",$scoring,"/UnMapped.bam FASTQ=",$scoring,"/read1.fq SECOND_END_FASTQ=",$scoring,"/read2.fq VALIDATION_STRINGENCY=SILENT > ",$logs,"/1a.bam2fastq.log 2>&1");
	submit($command);
}

sub mapagain{
	my ($BAMFILE, $human_mapping, $human_mapping_again, $THREADS, $PICARD, $logs)= @_;
	print "mapping the reads back to human again...\n";
	`bwa mem -t $THREADS -M -R \'\@RG\\tID:$SAMPLE\\tSM:$SAMPLE\' $HUMAN_database_Index $human_mapping/read1.fq $human_mapping/read2.fq | samtools view -u -bS - >  $human_mapping_again/human.bam`;
	`samtools view -u -b -f 8 -F 260 $human_mapping_again/human.bam > $human_mapping_again/oneEndMapped.bam`;
	`samtools view -u -b -f 4 -F 264 $human_mapping_again/human.bam > $human_mapping_again/oneEndUnMapped.bam`;
	print "#### converting the BAM file to fastq...\n";
	`perl $SCRIPT_DIR/keepreads.pl $human_mapping_again/human.sort.bam | sed '/^\s*$/d' | samtools view -u -bS - > $human_mapping_again/human.fix.sort.bam`;
	$command=join ("","java -Xmx6g -Xms3g -jar ",$PICARD,"/SamToFastq.jar INPUT=",$human_mapping_again,"/human.fix.sort.bam FASTQ=",$human_mapping_again,"/read1.fq SECOND_END_FASTQ=",$human_mapping_again,"/read2.fq VALIDATION_STRINGENCY=SILENT > ",$logs,"/3.bam2fastq.log 2>&1");
	submit($command);
}

sub mergeallbamsandextractreads{
	my ($BAMFILE, $human_mapping, $human_mapping_again, $THREADS, $PICARD, $logs) = @_;
	print "### merging and sorting the Human BAM file...\n";
	`samtools merge -u -n $human_mapping_again/2_human.bam $human_mapping_again/oneEndMapped_soft.bam $human_mapping_again/oneEndUnMapped_soft.bam $human_mapping/oneEndMapped.bam $human_mapping/oneEndUnMapped.bam`;
	`samtools sort -@ $THREADS -n $human_mapping_again/2_human.bam $human_mapping_again/human.sort`;
	#this removes any white space characters \s
	`perl $SCRIPT_DIR/keepreads.pl $human_mapping_again/human.sort.bam | sed '/^\s*$/d' | samtools view -u -bS - > $human_mapping_again/human.fix.sort.bam`;
	
	print "#### converting the BAM file to fastq...\n";
	$command=join ("","java -Xmx6g -Xms3g -jar ",$PICARD,"/SamToFastq.jar INPUT=",$human_mapping_again,"/human.fix.sort.bam FASTQ=",$human_mapping_again,"/read1.fq SECOND_END_FASTQ=",$human_mapping_again,"/read2.fq VALIDATION_STRINGENCY=SILENT > ",$logs,"/3.bam2fastq.log 2>&1");
	submit($command);
}

sub getviralcoverage{
	my ($coverage, $OUTNAME, $OUTPUT, $autocode, $VIRUS_database, $SCRIPT_DIR, $viral_mapping, $scoring) = @_;
    print "#### get viral coverage...";
	my $coverage=$OUTPUT ."/.coverage";
	mkdir $coverage;
	my $count = `wc -l < $OUTNAME`;
	sleep 10;
	
	if ( $count > 1 )	{
		open FH, ">$autocode/coverage.sh" or die "can't open the script to write\n";
		print FH "for i in `cat $OUTNAME | awk 'NR>1' | cut -f8 | sort | uniq`\n";
		print FH "do\n";
		print FH "cat $VIRUS_database.fai | grep -w \$i | cut -f1,2 | perl $SCRIPT_DIR/coverage_split.pl | coverageBed -abam $viral_mapping/virus.fix.sort.bam -b stdin > $coverage/\$i.coverage.out \&\n";
		print FH "pid=\$!\n";
		print FH "cat $VIRUS_database.fai | grep -w \$i | cut -f1,2 | perl $SCRIPT_DIR/coverage_split.pl | coverageBed -abam $scoring/virus.fromUnmapped.sort.bam -b stdin > $coverage/\$i.proper.coverage.out \&\n";
		print FH "pid1=\$!\n";
		print FH "wait \$pid \$pid1\n";
		print FH "done\n";
		close FH;
		$command=join("","chmod 777 ",$autocode,"/coverage.sh");
		submit($command);
		$command=join ("","sh ",$autocode,"/coverage.sh");
		submit($command);
		sleep 10;
		submit($command);
	}
	sleep 10;
}

 


