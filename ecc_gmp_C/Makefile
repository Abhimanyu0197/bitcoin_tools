default:
	gcc -O3 -c sha256/sha256.c -o sha256.o
	gcc -O3 -c base58/base58.c -o base58.o
	gcc -O3 -c rmd160/rmd160.c -o rmd160.o
	gcc -O3 -c bech32m/segwit_addr.c -o segwit_addr.o
	gcc -O3 -c util.c -o util.o
	gcc -O3 -o test_gmp test_gmp.c util.o base58.o sha256.o rmd160.o segwit_addr.o -lgmp
	#rm *.o
clean:
	#rm -r *.o
