default:
	gcc -O3 -c sha256/sha256.c -o sha256.o
	gcc -O3 -c base58/base58.c -o base58.o
	gcc -O3 -c rmd160/rmd160.c -o rmd160.o
	gcc -O3 -c bech32/segwit_addr.c -o segwit_addr.o
	gcc -O3 -c secp256k1/secp256k1.c -o secp256k1.o
	gcc -O3 -c util.c -o util.o
	gcc -O3 -o gmp_test main.c  secp256k1.o util.o base58.o sha256.o rmd160.o segwit_addr.o -lgmp
	#rm *.o
clean:
	#rm -r *.o
