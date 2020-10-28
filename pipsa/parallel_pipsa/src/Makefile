F77=g77
all:	
	$(F77)  -O -static 2potsim_noskin.f -o ../bin/2potsim_noskin
	$(F77)  -O -static 2potsim_skin.f -o ../bin/2potsim_skin
	$(F77)  -O -static 2potsim_skin_parts.f -o ../bin/2potsim_skin_parts
	$(F77)  -O -static ccenter.f -o ../bin/ccenter 
	$(F77)  -O -static grid_asc2bin.f -o ../bin/grid_asc2bin
	$(F77)  -O -static mkdismx.f -o ../bin/mkdismx
	$(F77)  -O -static mkkin.f -o ../bin/mkkin
	$(F77)  -O -static mkgridin.f -o ../bin/mkgridin
	$(F77)  -O -static mkapbsin.f -o ../bin/mkapbsin 
	$(F77)  -O -static mkuhbdin.f -o ../bin/mkuhbdin 
	$(F77)  -O -static modeller2grin.f -o ../bin/modeller2grin
	$(F77)  -O -static n1potsim.f -o ../bin/n1potsim
	$(F77)  -O -static npotsim.f -o ../bin/npotsim
	$(F77)  -O -static qdipsim.f -o ../bin/qdipsim
	$(F77)  -O -static uhbd_asc2bin.f -o ../bin/uhbd_asc2bin
	$(F77)  -O -static whatif2uhbd.f -o ../bin/whatif2uhbd
	cp mat2phyl.pl ../bin/
	cp mat2phyl40.pl ../bin/
