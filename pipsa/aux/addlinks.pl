#!/usr/bin/perl -w

# This program will add a link to every line that has '(.*) show' in a .ps file. 
# The first passed argument should be the .ps file name
# THe second arguemnt should be the link that would be substituted. 
# So if you write addlinks.pl infile.ps 'http://www.ble.com/, it will add links in the form http://www.ble.com/$1, where $1 is the text in the brackets.
# This program was used to obtain a tree output from Phylip that would have links

die("Usage: addlinks.pl infile.ps link") if  (@ARGV+0!=2);


open(INFILE, "<$ARGV[0]") or die("Can't open the file for reading: $ARGV[0]");
$linkbase=$ARGV[1];


# This piece of code is used, so that the PS device will not implement the pdfmark operator. (Consult pdfmark Reference Manual) (http://partners.adobe.com/asn/developer/acrosdk/DOCS/pdfmark.pdf)

$pdfmark="/pdfmark where
{pop} {userdict /pdfmark /cleartomark load put} ifelse
";

$perc=0;
# This is the link code (Consult pdfmark Reference Manual)




while(<INFILE>)
  {
    if (!/^%%/ && $perc==0) {print $pdfmark; $perc=1;}
    if (/\((.+)\) show/)
      {
	$link=$linkbase.$1;
print "[/Rect [ 0 0 10 5 ] 
/Action << /Subtype /URI /URI (".$link.") >>
/Border [ 0 0 1 [ 1 2 ] ]
/Color [ 1 0 0 ]
/Subtype /Link
/ANN pdfmark
";

      }
    print;
  }
