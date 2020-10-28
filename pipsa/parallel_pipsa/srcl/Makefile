all:	
	g77  -O -static 2potsim_noskin.f -o ../binl/2potsim_noskin
	g77  -O -static 2potsim_skin.f -o ../binl/2potsim_skin
	g77  -O -static 2potsim_skin_parts.f -o ../binl/2potsim_skin_parts
	g77  -O -static ccenter.f -o ../binl/ccenter 
	g77  -O -static grid_asc2bin.f -o ../binl/grid_asc2bin
	g77  -O -static mkdismx.f -o ../binl/mkdismx
	g77  -O -static mkkin.f -o ../binl/mkkin
	g77  -O -static mkgridin.f -o ../binl/mkgridin
	g77  -O -static mkapbsin.f -o ../binl/mkapbsin 
	g77  -O -static mkuhbdin.f -o ../binl/mkuhbdin 
	g77  -O -static modeller2grin.f -o ../binl/modeller2grin
	g77  -O -static n1potsim.f -o ../binl/n1potsim
	g77  -O -static npotsim.f -o ../binl/npotsim
	g77  -O -static qdipsim.f -o ../binl/qdipsim
	g77  -O -static uhbd_asc2bin.f -o ../binl/uhbd_asc2bin
	g77  -O -static whatif2uhbd.f -o ../binl/whatif2uhbd
	cp mat2phyl.pl ../binl/
	cp mat2phyl40.pl ../binl/
allPGF:	
	pgf77  -fast -Bstatic 2potsim_noskin.f -o ../binl/2potsim_noskinPGF
	pgf77  -fast -Bstatic 2potsim_skin.f -o ../binl/2potsim_skinPGF
	pgf77  -fast -Bstatic 2potsim_skin_parts.f -o ../binl/2potsim_skin_partsPGF
	pgf77  -fast -Bstatic ccenter.f -o ../binl/ccenterPGF
	pgf77  -fast -Bstatic grid_asc2bin.f -o ../binl/grid_asc2binPGF
	pgf77  -fast -Bstatic mkdismx.f -o ../binl/mkdismxPGF
	pgf77  -fast -Bstatic mkkin.f -o ../binl/mkkinPGF
	pgf77  -fast -Bstatic mkgridin.f -o ../binl/mkgridinPGF
	pgf77  -fast -Bstatic mkapbsin.f -o ../binl/mkapbsinPGF
	pgf77  -fast -Bstatic mkuhbdin.f -o ../binl/mkuhbdinPGF
	pgf77  -fast -Bstatic modeller2grin.f -o ../binl/modeller2grinPGF
	pgf77  -fast -Bstatic n1potsim.f -o ../binl/n1potsimPGF
	pgf77  -fast -Bstatic npotsim.f -o ../binl/npotsimPGF
	pgf77  -fast -Bstatic qdipsim.f -o ../binl/qdipsimPGF
	pgf77  -fast -Bstatic uhbd_asc2bin.f -o ../binl/uhbd_asc2binPGF
	pgf77  -fast -Bstatic whatif2uhbd.f -o ../binl/whatif2uhbdPGF
	cp mat2phyl.pl ../binl/
	cp mat2phyl40.pl ../binl/
