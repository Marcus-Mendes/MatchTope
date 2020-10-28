#!/usr/bin/perl -w

# The program converts the distance matrix created by pipsa into a matrix that can be used in phylip program.
# The result is output in STDOUT file. Program should be run from source directory. 

open(NAMES, "<./names") or die "Can't open file ./names.";

while(<NAMES>)
  {
    chop;
#    s/(.*)\.pdb.*/$1/;
    s/(.{1,40}).*/$1/g;
    $l = length;
    if ($l<40) { $_.=" "x(40-$l);}
    push(@names, $_);
  }
close(NAMES);

@matrix=<STDIN>;
($l = @matrix+0) or die "No entries present in file sims.mat";
print "\t$l\n";
print STDERR "The following entries are saved: \n";
foreach (@matrix)
  {
    $name = shift(@names);
    s/([ ]*)(.*)/$2/;
    print "$name$_";
    print STDERR $name;
    
  }
print STDERR "\n";
# close(PHYL);
#close(MAT);
