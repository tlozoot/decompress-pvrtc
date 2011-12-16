all: clean decompress-2

decompress-2:
	gcc -ggdb -Wall -std=c99 -o decompress decompress.c fortseige_pvr.c
	./decompress 2
	open img/firefox-2bpp.bmp

decompress-4:
	gcc -ggdb -Wall -std=c99 -o decompress decompress.c fortseige_pvr.c
	./decompress 4
	open img/firefox-4bpp.bmp

# Not yet working...
mypvr:
	gcc -ggdb -Wall -std=c99 -o decompress decompress.c my_pvr.c
	./decompress 4
	open img/firefox-4bpp.bmp

clean:
	rm -rf decompress *.o *.dSYM

