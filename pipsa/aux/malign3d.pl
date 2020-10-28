#!/usr/bin/perl -w

# The program runs modeller 6 to create 3D multiple alignment of all the proteins in pdbs directory 
# The new, aligned proteins will be present in ../malign3d/ directory.
# In principle, it does the same thing, as malign3d.sh, but using modeller, not whatif.
# Important: The .pdb files should have "\n" at the end of the line or modeller will complain with strange messages.
# Important: The modeller seems not to work on SGI O2s. But I think it works without problems on Indigo2. Not sure if it's a problem with OS or the processor.



chdir "../pdbs/" or die("Can't change directory to ../pdbs/");
# Remove the old .fit files... 

`rm -f *.fit`;


@files = glob("*.pdb");
die("There are no .pdb files in ../pdbs/ directory") if (@files+0==0);


# Firstly, let's create the script file for modeller...
open(OUT, ">/tmp/modeller.top") or die("Can't open the file modeller.in for writing. Script for the modeller can't be created");


print "Constructing modeller input file...\n";
$_=shift(@files);
s/(.*)\.pdb/$1/g;
print OUT "# This is a file created by malign3d.pl that will be used 
# as input for the modeller to align sequences.
# The next line is for debugging purposes
# SET OUTPUT_CONTROL = 1 1 1 1 0
SET MAXRES = 1200
READ_MODEL FILE = '$_'
SEQUENCE_TO_ALI ATOM_FILES = '$_', ALIGN_CODES = '$_'
";


while (@files)
  {
    $_=shift(@files);
    s/(.*)\.pdb/$1/g;
    print OUT "READ_MODEL FILE = '$_'
SEQUENCE_TO_ALI ADD_SEQUENCE = on, ATOM_FILES = ATOM_FILES '$_', ALIGN_CODES = ALIGN_CODES '$_'
";
  }
print OUT "MALIGN 
MALIGN3D FIT = on, WRITE_FIT = on,  WRITE_WHOLE_PDB = on
";
close(OUT);

# And now let's run modeller...

print "Running modeller...\n";
`mod6a /tmp/modeller.top`;
die("Modeller exited with an error. Check log file.") unless ($?==0);


print "Moving the aligned files (.fit files) to ../malign3d/ directory...\n";
mkdir("../malign3d/", 0755);
`rm -f ../malign3d/*.pdb`; # one of them should do the stuff.

@fitfiles=glob("*fit.pdb");

while(@fitfiles)
  {
    $_ = ($fit = shift(@fitfiles));
    s/(.*)_fit.pdb/$1.pdb/g;
    `mv $fit ../malign3d/$_`;
  }
`cp /tmp/modeller.log ../malign3d/`; # copy the log file.

print "Aligned files can be found in the ../malign3d/ directory\n";
print "Done.\n";
