#!/usr/bin/perl

#marked some comments with my nick name 'bug', check those!

use Getopt::Std;

#@dist_methods=('root','euclidean', 'manhattan', 'minkowski_p2', 'maximum');
#dist_methods=('minkowski_p2')best results so far;
@dist_methods=('maximum');
#@clust_methods=('ward', 'average', 'complete', 'single' ,'median', 'mcquitty', 'centroid');
@clust_methods=('complete');
#clust_methods=('median', 'single')best results so far;
@method_dist=('correlation');
#clust_methods=('complete', 'abscor')best results so far;

getopts("s:t:o:dgebka:m:n:h:f");
if (length($opt_d)>0) {
	$maxCarbo = 0;
	$maxMax = 2;
	$minCarbo = 2;
	$minMin = 0;
} else {
	$maxCarbo = -1.0;
	$minCarbo = 1.0;
	$maxMax=1.0;
	$minMin=-1.0;
}
if (length($opt_n)>0) {
    if ($opt_n!~/^[\d]+$/) {
        $opt_n=2;
    }
} else {
    $opt_n=2;
}
unless (-s $opt_s >0) {
	help();
	die "File $opt_s not found";
}
unless ($opt_t eq 'h' or $opt_t eq 'c' or $opt_t eq 'cs' or $opt_t eq 'hs' or $opt_t eq 'cse' or $opt_t eq 'd' or $opt_t eq 'e') {
	help();
	die "Valid options for the type of the analysis are 'c' for carbo and 'h' for hodgkin";
}
if (defined($opt_g)) {
	print STDERR "Generating heatmap requested\n";
	unless ($opt_m) {
		print STDERR "Heatmap can only be specfied if an output matrix file name is given (option -m)\n";
		exit();
	}
	if (defined($opt_b)) {
		print STDERR "Heatmap will be in black and white\n";
	}
}
if ($opt_t eq 'h') {
	print STDERR "Hodgkin type analyisis requested\n";
}
if ($opt_t eq 'c') {
	print STDERR "Carbo type analyisis requested\n";
}
if ($opt_t eq 'd') {
	print STDERR "Potential differences analysis requested\n";
	@dist_methods=('asis');
}
if ($opt_t eq 'e') {
	print STDERR "Exp. potential differences analysis requested\n";
	@dist_methods=('asis');
}
if (defined($opt_a)) {
	if ($opt_a ne 'diana' and $opt_a ne 'agnes') {
		help();
		exit();
	}
} else {
	$opt_a="hclust";
}
my $data = parseLog($opt_s);
my %data = %$data;
my @proteins = ();
if (length($opt_o)>0) {
	print STDERR "Order according to '$opt_o'\n";
	@proteins = parseTreeFile();
} else {
	@proteins = keys %data;
}
print STDERR "Writing matrix to ${opt_m}\n";
if ($opt_m) {
	open FILE,">$opt_m";
	select(FILE);
}
print "\t";
foreach my $prot1 (@proteins) {
	print $prot1,"\t";
}
print "\n";
foreach my $prot1 (@proteins) {
	print $prot1,"\t";
	foreach my $prot (@proteins) {
		print $data{$prot}{$prot1},"\t";
	}
	print "\n";
}
if ($opt_m) {
	select(STDOUT);
	close(FILE);
}
if ($opt_d) {
	print STDERR "Distances range from $minCarbo to $maxCarbo (maximum range is from 0 to 2)\n";
} else {
	print STDERR "Similarities range from $minCarbo to $maxCarbo (maximum range is from -1 to 1)\n";
}
if (defined($opt_g)) {
	generateHeatMap();
}

sub help()
{
	print "Usage: pipsa2R.pl -s <sims.log> -t <h|c> -o <treefile>\n";
	print "\t-s\tsims.log: Pipsa output file\n";
	print "\t-t\th: Hodgkin index, c: Carbo index, cse: geometric mean carbo shape and carbo electrostatics\n\t\t\tcs: carbo index of shape only, hs: hodgkin index of shape only\n";
	print "\t-o\ttreefile: To specify the order of the entries in the matrix.\n";
	print "\t-d\tCalculate distance from the similarities (sqrt(1-sim)).\n";
	print "\t-m\tmatrixfile: Specify the output file for the matrix written, defaults to STDOUT\n";
	print "\t-g\tProduce heatmap (-m must be given)\n";
	print "\t-b\tProduce black and white heatmap (-g and -m must be given)\n";
	print "\t-n\tn\tNumber of clusters to create in the heatmap\n";
	print "\t-h\ttitle\tInserts the string as heading in the heatmap\n";
	print "\t-a\tagnes or diana\tuse agnes or diana clustering method, without -a hclust is used\n";
	print "\t-k\t\tkeep original names, do not remove _templatename for swissprot entries\n";
        print "\t-f\t\tfixed range (0-2) for the color code\n";
}

# sample format of sims.log file, shortened to the first rows, starting in column 6 after "|" (columns 13 - 15 only with spheres?)
#   |#	12
#   |TPIS_TRYCR TPIS_TRYCR
#   |  1.000  1.000  0.2880725E+06  0.2880725E+06  0.2880725E+06  0.2880725E+06  0.2880725E+06  1.000  1.000  34287  34287  34287  0.0000000E+00			NAN  0.0000000E+00
#   |TPI1_GIALA TPIS_TRYCR
#   |  0.215  0.217  0.4649634E+05  0.3569861E+05  0.8849652E+04  0.5656893E+05  0.2880725E+06  0.722  0.722  34373  34287  24796 -0.2338139E+00			INF  0.1875099E-01
#   |TPI1_GIALA TPI1_GIALA
#   |  1.000  1.000  0.5656893E+05  0.5656893E+05  0.5656893E+05  0.5656893E+05  0.5656893E+05  1.000  1.000  34373  34373  34373  0.0000000E+00			NAN  0.0000000E+00
#   |TPIS_CHICK TPIS_TRYCR
#   |  0.107  0.138  0.1697193E+05  0.7448458E+05  0.4889585E+04  0.2368744E+05  0.2880725E+06  0.750  0.750  33020  34287  25242 -0.4189886E+00		   -INF -0.3251152E+02
#   |TPIS_CHICK TPI1_GIALA
#   |  0.446  0.456  0.1850112E+05  0.2812830E+05  0.1039900E+05  0.2368744E+05  0.5656893E+05  0.712  0.712  33020  34373  23999 -0.1375701E+00 -0.5252399E+02  0.4709018E+00
# regexp matches:
#	  -$1--  -$2--  -----$3------  -----$4------  -----$5------  -----$6------  -----$7------  -$8--  -$9--  -$10-  -$11-  -$12- -----$13?----- -----$14?----- -----$15?-----

sub parseLog()
{
	my $file = shift;												# get passed parameter
	open(FILE, $file) or die "Can't open '$file': $!";				# open sims.log file
	my $n=0;
	my $carbo;
	my %prots;
	my $prot1 = undef;
	my $prot2 = undef;
	my %data;														# output
	while (<FILE>)													# go through file
	{
		if (/^#\s+(\d+)/) {											# matches first line of log-file, e.g. "#	12"
			  $n = $1;												# then n is number of compared structures, NEVER READ (bug)
		}
		if (/^\s+([-.\d]+)\s+([-.\d]+)\s+([-+.\deE]+)\s+([-+.\deE]+)\s+([-+.\deE]+)\s+([-+.\deE]+)\s+([-+.\deE]+)\s+([-.\d]+)\s+([-.\d]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([-+.\deE]+)?\s*([-+.\deE]+)?\s*([-+.\deE]+)?/)	# list of  numbers
		{
			if ($opt_t eq 'h') {									# if hodgkin
				$carbo = $1;										# carbo is first entry
			}
			if ($opt_t eq 'c') {									# if carbo
				$carbo = $2;										# carbo is second entry
			}
			if ($opt_t eq 'cse') {
				$carbo = sprintf("%.5f",$2+$9-$2*$9);
#				if ($2 != 0.0) {
#					$sign = $2/abs($2);
#					$carbo = sprintf("%.5f", $sign*sqrt(abs($2)*$9));  #geometric mean
#				} else {
#					$carbo = 0;
#				}
			}
			if ($opt_t eq 'cs') {
				$carbo = $9;
			}
			if ($opt_t eq 'hs')	{
				$carbo = $8;
			}
			if ($opt_t eq 'e') {
				$carbo = exp(abs($13));
			}
			if (length($opt_d)>0) {									# mit max output?
				if ($opt_t eq 'd') {
					$carbo = abs($carbo);
				} else {
					if ($opt_t eq 'e') {
						$carbo = exp(abs($carbo));					#bug: see above, now $carbo = exp(abs(exp(abs($13)))); on purpose?
						                                            #yes, this was an option later removed in the interface, only available if type set to 'e' 
					} else {
						$carbo = sprintf("%.5f", sqrt(2.0-2.0*$carbo));
					}
				}
			}
			if ($carbo < $minCarbo and  $carbo>$minMin) {			# refresh min
				$minCarbo = $carbo;
			}
			if ($carbo > $maxCarbo  and  $carbo<$maxMax) {			# refresh max
				$maxCarbo = $carbo;
			}
			if ($prot1 ne undef and $prot2 ne undef)				# if matrix entry != undef
			{
				$data{$prot1}{$prot2} = $carbo;						# set ab = ba
				if ($opt_t eq 'd') {								# bug: same statement in if and else
					$data{$prot2}{$prot1} = $carbo;                 # yes here the problem was that Razif had introduced this 
				} else {                                            #measure of the average difference of electrostatic potential which also could be negative
					$data{$prot2}{$prot1} = $carbo;                 #depending on the direction of comparion. But this turned out to ruin the 
				}                                                   #distance matrix. Therefore here the absolute distance (see above) is used.
		 	}
		 	$prot1=undef;
		 	$prot2=undef;
		}

		#if (/\s*([a-zA-Z]\w+)\s+([a-zA-Z]\w+)/)					# two times one letter, then alpha num
																	# JUST SWISSPROT SYNTAX!
		#if (/\s*(\w+)\s+(\w+)/) {

		#if (/^([\-a-zA-Z_\.\d]+)\s+([\-a-zA-Z_\.\d]+)$/) {
		if (/^([-.\w]+)\s+([-.\w]+)$/) {							#same as above
			$prot1 = $1;											# prot 1 as first entry
			$prot2 = $2;											# prot 2 as second entry
			$prot1 =~ s/addH//g;									# remove addH
			$prot2 =~ s/addH//g;
			if (!$opt_k) {
				$prot1 =~ s/([\w\d]{6,6})_([\w\d]{4,5})/$1/g;
				$prot2 =~ s/([\w\d]{6,6})_([\w\d]{4,5})/$1/g;
			}
#			 if(length $prot1 > 6)
#			 {
#				 $prot1 = substr($prot1,0,6);
#			 }
#			 if(length $prot2 > 6)
#			 {
#				 $prot2 = substr($prot2,0,6);
#			 }
			 $prots{$prot1}=1;
			 $prots{$prot2}=1;
		}
	}
	return \%data;
}

sub parseTreeFile()
{
	open (TREE, $opt_o) or die "Cannot open $opt_o";
	my $tree = "";
	my $dum = "";
	my @prots=();
	while (<TREE>) {
		chomp($_);
		$tree .= $_;
	}
	#while($tree =~ /[^:\.\)\d]([A-Za-z0-9]{4,6})/g)
	#while($tree =~ /[^:\.\)\d](\w+)/g)
	while ($tree =~ /([^\(\),]+):/g) {
		push @prots,($1);
	}
	return @prots;
}

sub generateHeatMap() {
	my ($dist_method,$clust_method);
	if (length(@clust_methods)<1) {
		@clust_methods=('ward');
	}
	if (length(@dist_methods<1)) {
		@dist_methods=("root");
	}
	my $a = "";
#	if (-s "${opt_m}_${opt_n}_${opt_t}.pdf" and -s "${opt_m}_${opt_n}_${opt_t}.png" ) {
#                print "PDF or PNG file already exists, NO NEW ONE IS GENERATED\n";
#	 	return;					#bug: why commented? calculates mat-files again, but does not have to draw the whole thing again
#		                        #because at the last moment i had to introduce the possiblitiy to draw no tree. ...
#	}
	if ($opt_t eq 'h') {
		$type = "Hodgkin";
	} else {
		$type = "Carbo";
	}
	open (RSCRIPT,">heatmap.R");
	print RSCRIPT "library(gplots)\n";
	print RSCRIPT "library(cluster)\n";
	print RSCRIPT "library(pvclust)\n";
#	print RSCRIPT "pdf(file=\"${opt_m}_${opt_t}_${opt_a}${opt_n}.pdf\",width=8,height=8,paper=\"a4\")\n";
	print RSCRIPT "pdf(file=\"${opt_m}_${opt_n}_${opt_t}.pdf\",width=9,height=9,paper=\"a4\")\n";
	print RSCRIPT "n=$opt_n\n";
	foreach $dist_method (@dist_methods) {			# bug: array only contains one entry, either "root" or "asis"
                print STDERR "Using distance: '$dist_method'\n";
		foreach $clust_method (@clust_methods) {	# bug: array only contains one entry, "ward"
                       print STDERR "Using clustering method '$clust_method'\n";
			print RSCRIPT "x=as.matrix(read.table(\"$opt_m\"))\n";
			print RSCRIPT "write.csv(x, file='table.csv')\n";
#			my $title='expression(paste("Electrostatic Distance (",D[list(a,b)]==sqrt(2-2*SI[list(a,b)]),")" sep=""))' ; #"\"$opt_h clust. method: $clust_method,$dist_method sim. index: $type, no. cluster: ${opt_n}\"";
			my $title = '';
			if ($opt_d) {
				print RSCRIPT "x.dist=as.dist(x)\n";
				if ($dist_method eq 'root') {
					$title = 'expression(paste("Electrostatic Distance ",D[list(a,b)]==sqrt(2-2*SI[list(a,b)])))';
				} elsif ($dist_method eq 'asis') {
					if ($opt_t eq 'd') {
			 			$title = "\"Absolute average electrostatic potential difference\"";			 			
					} else {
						$title = "\"Exp average electrostatic potential difference\"";
		  			}
				}
			print RSCRIPT "results <- pvclust(order, method.hclust='average', method.dist='cor', nboot=1000, parallel=FALSE)\n";
			} elsif ($dist_method ne 'root') {
				$title = "\"Electrostatic Distance $dist_method\"";
				if ($dist_method eq 'minkowski_p2') {
					print RSCRIPT "x.dist=dist(x,method=\"minkowski\",p=2)\n";
				} elsif ($dist_method eq 'asis') {
					$title = "\"Average electrostatic potential difference\"";
					print RSCRIPT "x.dist=as.dist(sqrt(abs(x)))\n";
				} else {
					print RSCRIPT "x.dist=dist(x,method=\"$dist_method\")\n";
				}
			} else {
				print RSCRIPT "x.dist=as.dist(sqrt(2-2*x))\n";
				$title = 'expression(paste("Electrostatic Distance ",D[list(a,b)]==sqrt(2-2*SI[list(a,b)])))';
			}
			print RSCRIPT "title=$title\n";
			if ($opt_a eq 'agnes') {
				print RSCRIPT "x.agnes=agnes(x.dist,diss=TRUE,method=\"$clust_method\")\n";
				print RSCRIPT "x.clust=as.hclust(x.agnes)\n";
			} elsif ($opt_a eq 'diana') {
				print RSCRIPT "x.diana=diana(x.dist,diss=TRUE)\n";
				print RSCRIPT "x.clust=as.hclust(x.diana)\n";
	 		} else {
				print RSCRIPT "x.clust=hclust(x.dist,method=\"$clust_method\")\n";
			}
			print RSCRIPT "x.dend=as.dendrogram(x.clust)\n";
			print RSCRIPT "x.clust.members=cutree(x.clust,k=n)\n";
			if ($opt_b) {
				print RSCRIPT "x.clust.coltypes=gray(1:n/(n+1))\n";
			} else {
				print RSCRIPT "x.clust.coltypes=rainbow(n,s=0.4,start=0/6,end=5/6)\n";
			}
			print RSCRIPT "x.clust.colors=1:ncol(x)\n";
			print RSCRIPT "order=order.dendrogram(x.dend)\n";
			print RSCRIPT "for(count in 1:ncol(x)) {\n";
			print RSCRIPT "   x.clust.colors[count]=x.clust.coltypes[x.clust.members[count]]\n";
			print RSCRIPT "}\n";
			print RSCRIPT "##### Snippet taken from SLmisc package to calculate constant colors for data\n";
			if ($opt_b) {
				print RSCRIPT "col=rev(gray(1:512 / 550))\n";
			} else {
#		  		print RSCRIPT "col=rev(heat.colors(256))\n";
		  		print RSCRIPT "col=rainbow(512,s=0.8,v=1,start=0/6,end=5/6)\n";
			}
                        if ($opt_f) {
                            print STDERR "Using abolute scaling\n";
			    print RSCRIPT "nrcol = length(col)\n";
			    print RSCRIPT "x.dist.range = range(as.matrix(x.dist))\n";
			    print RSCRIPT "ulim = 2\n";
			    print RSCRIPT "llim = 0\n";
			    print RSCRIPT "reps1 = ceiling(nrcol*(llim-x.dist.range[1])/(2*ulim))\n";
			    print RSCRIPT "reps2 = ceiling(nrcol*(x.dist.range[2]-llim)/(ulim))\n";
			    print RSCRIPT "coln = col[1:reps2]# , col, rep(col[nrcol],reps2))\n";
                        } else {
                            print STDERR "Using relative scaling\n";
  			    print RSCRIPT "coln = col\n";
                        }
			if ($opt_k) {
				print RSCRIPT "margin=c(12,12)\n";
			} else {
				print RSCRIPT "margin=c(5,5)\n";
			}
			print RSCRIPT "#####\n";
			if ($opt_e) {
				print RSCRIPT "heatmap.2(as.matrix(x.dist),dendrogram=\"none\",margins=c(14,14),col=coln,scale=\"none\",distfun=function(n) n,hclustfun=function(x) hclust(x.dist,method=\"$clust_method\"),trace=\"none\",symkey=FALSE,denscol=\"white\",main=$title)\n";
			} else {
				print RSCRIPT "heatmap.2(as.matrix(x.dist),margins=c(14,14),col=coln,scale=\"none\",distfun=function(n) n,Colv=x.dend,Rowv=x.dend,hclustfun=function(x) hclust(x.dist,method=\"$clust_method\"),trace=\"none\",revC=TRUE,symkey=FALSE,denscol=\"white\",density.info='density',densadj=1,ColSideColors=x.clust.colors,main=$title)\n";
			}
			print RSCRIPT "rm(x)\n";
		}
	}
	print RSCRIPT "results <- pvclust(read.table(\"$opt_m\"), method.hclust='complete', method.dist='cor', nboot=10000, parallel=FALSE)\n";
	print RSCRIPT "plot(results)\n";
	print RSCRIPT "output <- pvrect(results, alpha=0.97)\n";
	print RSCRIPT "dev.off()\n";
	close (RSCRIPT);
#	print STDERR `R --vanilla < heatmap.R`;
	`R --vanilla < heatmap.R`;
	`pstoimg -crop all -out ${opt_m}_${opt_n}_${opt_t}.png ${opt_m}_${opt_n}_${opt_t}.pdf`;
#	`ps2epsi ${opt_m}_${opt_n}.ps`;
#	`convert ${opt_m}_${opt_n}.epsi ${opt_m}_${opt_n}.png`;
}
