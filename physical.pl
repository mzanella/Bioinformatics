#!/usr/bin/perl
die"Sintax: physical.pl samfile genomelen" unless $#ARGV>0; # if no arguments print instructions
@genome_change = (0) x $ARGV[1];                            # initialize genome_change, all to 0
 
open(INFILE, $ARGV[0]);                                     # open sam file: INFILE as file handler
while($line = <INFILE>)                                     # put in $line the next line of the file
{   if($line =~ /^@/) {next;}                               # if $line starts with @ then skip it
    @item = split("\t", $line);                             # split $line at tab, put results in @item
    if((($item[1] & 3) == 3) && ($item[8] > 0) )            # both reads align correctly and length > 0
    {  $genome_change[$item[3]]++;                          # increment start position
       $genome_change[ $item[7] + length($item[9]) ]--;     # decrement end position
    } 
}
close INFILE;                                               # the sam file has been fully read

print("fixedStep chrom=genome start=1 step=1 span=1\n");    # print the heading of the wiggle file
$current_coverage=0;                                        
for($i=0; $i<$ARGV[1]; $i++)                                # output the results
{   $current_coverage += $genome_change[$i]; 
    print "$current_coverage\n"; 
}

