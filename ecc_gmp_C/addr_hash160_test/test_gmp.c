#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
//#include <time.h>
#include <sys/time.h>
#include <gmp.h>

#include "base58/libbase58.h"
#include "rmd160/rmd160.h"
#include "sha256/sha256.h"
//#include "bech32/bech32.h"
#include "bech32m/segwit_addr.h"
#include "util.h"

struct Elliptic_Curve {
    mpz_t a;
    mpz_t b;
    mpz_t p;
    mpz_t n;
};

struct Point {
    mpz_t x;
    mpz_t y;
};

struct Elliptic_Curve EC; 
struct Point Curve_G;

void Point_Doubling(struct Point P, struct Point *R) {
    mpz_t slope, temp;
    mpz_init(temp); mpz_init(slope);
    
    mpz_mul_ui(temp, P.y, 2);
    mpz_invert(temp, temp, EC.p);
    mpz_mul(slope, P.x, P.x);
    mpz_mul_ui(slope, slope, 3);
    mpz_add(slope, slope, EC.a);
    mpz_mul(slope, slope, temp);
    mpz_mod(slope, slope, EC.p);
    
    mpz_mul(R->x, slope, slope);
    mpz_sub(R->x, R->x, P.x);
    mpz_sub(R->x, R->x, P.x);
    mpz_mod(R->x, R->x, EC.p);
    
    mpz_sub(temp, P.x, R->x);
    mpz_mul(R->y, slope, temp);
    mpz_sub(R->y, R->y, P.y);
    mpz_mod(R->y, R->y, EC.p);
    
    mpz_clear(temp); mpz_clear(slope);
}

void Point_Addition(struct Point P, struct Point Q, struct Point *R) {
    struct Point S;
    mpz_init(S.x); mpz_init(S.y);
    mpz_set(S.x, P.x); mpz_set(S.y, P.y);
    
    if(mpz_cmp_ui(S.x, 0) == 0 && mpz_cmp_ui(S.y, 0) == 0) {
        mpz_set(R->x, Q.x);
        mpz_set(R->y, Q.y);
        return;
    }
    if(mpz_cmp_ui(Q.x, 0) == 0 && mpz_cmp_ui(Q.y, 0) == 0) {
        mpz_set(R->x, S.x);
        mpz_set(R->y, S.y);
        return;
    }
    
    mpz_t temp; mpz_init(temp);
    
    if(mpz_cmp_ui(Q.y, 0) != 0) {
        mpz_sub(temp, EC.p, Q.y);
        mpz_mod(temp, temp, EC.p);
    } else {
        mpz_set_ui(temp, 0);
    }
    
    if(mpz_cmp(S.y, temp) == 0 && mpz_cmp(S.x, Q.x) == 0) {
        mpz_set_ui(R->x, 0);
        mpz_set_ui(R->y, 0);
        mpz_clear(temp);
        return;
    }
    if(mpz_cmp(S.x, Q.x) == 0 && mpz_cmp(S.y, Q.y) == 0) {
        Point_Doubling(S, R);
        mpz_clear(temp);
        return;
    } else {
        mpz_t slope;
        mpz_init_set_ui(slope, 0);
        
        mpz_sub(temp, S.x, Q.x);
        mpz_mod(temp, temp, EC.p);
        
        mpz_invert(temp, temp, EC.p);
        mpz_sub(slope, S.y, Q.y);
        mpz_mul(slope, slope, temp);
        mpz_mod(slope, slope, EC.p);
        
        mpz_mul(R->x, slope, slope);
        mpz_sub(R->x, R->x, S.x);
        mpz_sub(R->x, R->x, Q.x);
        mpz_mod(R->x, R->x, EC.p);
        
        mpz_sub(temp, S.x, R->x);
        mpz_mul(R->y, slope, temp);
        mpz_sub(R->y, R->y, S.y);
        mpz_mod(R->y, R->y, EC.p);
        
        mpz_clear(temp);
        mpz_clear(slope);
        mpz_clear(S.x);
        mpz_clear(S.y);
        return;
    }
}

void Scalar_Multiplication(struct Point *R, mpz_t m) {
    struct Point P;
    mpz_init(P.x); mpz_init(P.y);
    mpz_set(P.x, Curve_G.x); mpz_set(P.y, Curve_G.y);
    struct Point Q, T;
    mpz_init(Q.x); mpz_init(Q.y); mpz_init(T.x); mpz_init(T.y);
    long no_of_bits, loop;
    no_of_bits = mpz_sizeinbase(m, 2);
    mpz_set_ui(R->x, 0); mpz_set_ui(R->y, 0);
    if(mpz_cmp_ui(m, 0) == 0) return;
    
    mpz_set(Q.x, P.x); mpz_set(Q.y, P.y);
    
    if(mpz_tstbit(m, 0) == 1) {
        mpz_set(R->x, P.x);
        mpz_set(R->y, P.y);
    }
    
    for(loop = 1; loop < no_of_bits; loop++) {
        mpz_set_ui(T.x, 0); mpz_set_ui(T.y, 0);
        Point_Doubling(Q, &T);
        mpz_set(Q.x, T.x); mpz_set(Q.y, T.y);
        mpz_set(T.x, R->x); mpz_set(T.y, R->y);
        if(mpz_tstbit(m, loop)) Point_Addition(T, Q, R);
	}
    
    mpz_clear(P.x); mpz_clear(P.y);
    mpz_clear(Q.x); mpz_clear(Q.y);
    mpz_clear(T.x); mpz_clear(T.y);    
}

void Point_Multiplication(struct Point P, struct Point *R, mpz_t m) {
    struct Point Q, T;
    mpz_init(Q.x); mpz_init(Q.y); mpz_init(T.x); mpz_init(T.y);
    long no_of_bits, loop;
    no_of_bits = mpz_sizeinbase(m, 2);
    mpz_set_ui(R->x, 0); mpz_set_ui(R->y, 0);
    if(mpz_cmp_ui(m, 0) == 0) return;
    
    mpz_set(Q.x, P.x); mpz_set(Q.y, P.y);
    
    if(mpz_tstbit(m, 0) == 1) {
        mpz_set(R->x, P.x);
        mpz_set(R->y, P.y);
    }
    
    for(loop = 1; loop < no_of_bits; loop++) {
        mpz_set_ui(T.x, 0); mpz_set_ui(T.y, 0);
        Point_Doubling(Q, &T);
        mpz_set(Q.x, T.x); mpz_set(Q.y, T.y);
        mpz_set(T.x, R->x); mpz_set(T.y, R->y);
		if(mpz_tstbit(m, loop)) Point_Addition(T, Q, R);
    }
    
    mpz_clear(Q.x); mpz_clear(Q.y);
    mpz_clear(T.x); mpz_clear(T.y);
}

void Point_Negation(struct Point *A) {
    mpz_sub(A->y, EC.p, A->y);
}

void Point_Subtraction(struct Point A, struct Point *B, struct Point *R) {
    if (mpz_cmp(A.x, B->x) == 0 && mpz_cmp(A.y, B->y) == 0) {
        mpz_set_ui(R->x, 0);
        mpz_set_ui(R->y, 0);
        return;
    } else {
        struct Point Q;
        Point_Negation(B);
        mpz_init(Q.x); mpz_init(Q.y);
        mpz_set(Q.x, B->x); mpz_set(Q.y, B->y);
        Point_Addition(A, Q, R);
        mpz_clear(Q.x); mpz_clear(Q.y);
        return;      
   }
}

void Point_Division(struct Point D, struct Point *S, mpz_t k) {
    mpz_t m; mpz_init(m);
    mpz_invert(m, k, EC.n);
    Point_Multiplication(D, S, m);
    mpz_clear(m);
}

bool Point_On_Curve(struct Point A) {
    mpz_t X, Y;
    mpz_init(X); mpz_init(Y);
    mpz_pow_ui(X, A.x, 3);
    mpz_add(X, X, EC.b);
    mpz_mod(X, X, EC.p);
    mpz_pow_ui(Y, A.y, 2);
    mpz_mod(Y, Y, EC.p);
    if (mpz_cmp(X, Y) == 0) {
        mpz_clear(X); mpz_clear(Y);
        return true;
    } else {
        mpz_clear(X); mpz_clear(Y);
        return false;
    }
}

static char upub[132];
static char cpub[68];
static char address[50];
static char hash160[42];
static char bech32_output[128];
static char wif[54];
static char bin_publickey[65];
static char bin_sha256[32];
static char bin_digest[60];

const char * Point_To_Upub(struct Point A) {
    gmp_snprintf(upub, 132, "04%0.64Zx%0.64Zx", A.x, A.y);
    return upub;
}

const char * Point_To_Cpub(struct Point A) {
    if(mpz_tstbit(A.y, 0) == 0) { 
        gmp_snprintf(cpub, 67, "02%0.64Zx", A.x); 
    } else { 
        gmp_snprintf(cpub, 67, "03%0.64Zx", A.x);
    }
    return cpub;
}

const char * Point_To_Hash160(struct Point pubkey, bool compressed) {
    char bin_rmd160[20];
    if(compressed) {
        if(mpz_tstbit(pubkey.y, 0) == 0) {
            gmp_snprintf (cpub, 67, "02%0.64Zx", pubkey.x);
        } else {
            gmp_snprintf(cpub, 67, "03%0.64Zx", pubkey.x);
		}
        hexs2bin(cpub, bin_publickey);
        sha256(bin_publickey, 33, bin_sha256);
	}
    else {
        gmp_snprintf(upub, 132, "04%0.64Zx%0.64Zx", pubkey.x, pubkey.y);
        hexs2bin(upub, bin_publickey);
        sha256(bin_publickey, 65, bin_sha256);
	}
    RMD160Data((const unsigned char*)bin_sha256, 32, bin_rmd160);
    tohex_dst(bin_rmd160, 20, hash160);
    return hash160;    
}

const char * Point_To_Legacy_Address(struct Point pubkey, bool compressed) {
    size_t pubaddress_size = 50;
    if(compressed) {
        if(mpz_tstbit(pubkey.y, 0) == 0) {
            gmp_snprintf(cpub, 68, "02%0.64Zx", pubkey.x);
		} else {
            gmp_snprintf(cpub, 68, "03%0.64Zx", pubkey.x);
        }
        hexs2bin(cpub, bin_publickey);
        sha256(bin_publickey, 33, bin_sha256);
    } else {
        gmp_snprintf(upub, 132, "04%0.64Zx%0.64Zx", pubkey.x, pubkey.y);
        hexs2bin(upub, bin_publickey);
        sha256(bin_publickey, 65, bin_sha256);
    }
    RMD160Data(bin_sha256, 32, bin_digest + 1);
    bin_digest[0] = 0;
    sha256(bin_digest, 21, bin_digest + 21);
    sha256(bin_digest + 21, 32, bin_digest + 21);
    b58enc(address, &pubaddress_size, bin_digest, 25);
    return address;
}

const char * Point_To_P2SH_Address(struct Point pubkey) {
    size_t pubaddress_size = 50;
    if(mpz_tstbit(pubkey.y, 0) == 0) {
        gmp_snprintf(cpub, 68, "02%0.64Zx", pubkey.x);
    }
    else {
        gmp_snprintf(cpub, 68, "03%0.64Zx", pubkey.x);
    }
    hexs2bin(cpub, bin_publickey);
    sha256(bin_publickey, 33, bin_sha256);
    RMD160Data(bin_sha256, 32, bin_digest + 2);
    bin_digest[0] = 0x00; bin_digest[1] = 0x14;
    sha256(bin_digest, 22, bin_sha256);
    RMD160Data(bin_sha256, 32, bin_digest + 1);
    bin_digest[0] = 0x05;
    sha256(bin_digest, 21, bin_sha256);
    sha256(bin_sha256, 32, bin_sha256);
    bin_digest[21] = bin_sha256[0];
    bin_digest[22] = bin_sha256[1];
    bin_digest[23] = bin_sha256[2];
    bin_digest[24] = bin_sha256[3];
    b58enc(address, &pubaddress_size, bin_digest, 25);
    return address;
}

const char * Point_To_Bech32_P2WPKH_Address(struct Point pubkey) {
    if(mpz_tstbit(pubkey.y, 0) == 0) {
        gmp_snprintf(cpub, 68, "02%0.64Zx", pubkey.x);
    }
    else {
        gmp_snprintf(cpub, 68, "03%0.64Zx", pubkey.x);
    }
    hexs2bin(cpub, bin_publickey);    
    sha256(bin_publickey, 33, bin_sha256);
    RMD160Data(bin_sha256, 32, bin_digest);
    segwit_addr_encode(bech32_output, "bc", 0, bin_digest, 20);
    return bech32_output;
}

const char * Point_To_Bech32_P2WSH_Address(struct Point pubkey) {
    if(mpz_tstbit(pubkey.y, 0) == 0) {
        gmp_snprintf(cpub, 68, "02%0.64Zx", pubkey.x);
    }
    else {
        gmp_snprintf(cpub, 68, "03%0.64Zx", pubkey.x);
    }
    hexs2bin(cpub, bin_publickey);
    sha256(bin_publickey, 33, bin_sha256);
    segwit_addr_encode(bech32_output, "bc", 0, bin_sha256, 32);
    return bech32_output;
}

char tap_tweak[] = "TapTweak";

const char * Tagged_Hash(mpz_t x) {
    char hex_mpz[66];
    char pub_bin[33];
    char taghash[33];
    char taghash_final[99];
    gmp_snprintf(hex_mpz, 65, "%0.64Zx", x);
    hexs2bin(hex_mpz, pub_bin);
    sha256(tap_tweak, strlen(tap_tweak), taghash);
    memmove(taghash_final, taghash, sizeof(taghash));
    memmove(taghash_final + 32, taghash, sizeof(taghash));
    memmove(taghash_final + 64, pub_bin, sizeof(pub_bin));
    sha256(taghash_final, 96, bin_sha256);
    return tohex(bin_sha256, 32);
}

const char * Taproot_Tweak_PrivKey(mpz_t k) {
    mpz_t seckey, t, s;
    mpz_init(seckey); mpz_init(t); mpz_init(s);
    struct Point Q;
    mpz_init(Q.x);
    mpz_init(Q.y);
    Scalar_Multiplication(&Q, k);
    if(mpz_tstbit(Q.y, 0) == 0) {
        mpz_set(seckey, k);
    }
    else {
        mpz_sub(s, EC.n, k);
        mpz_set(seckey, s);
    }
    mpz_set_str(t, Tagged_Hash(Q.x), 16);
    mpz_add(seckey, seckey, t);
    mpz_mod(seckey, seckey, EC.n);
    char * number_str = mpz_get_str(NULL, 16, seckey);
    mpz_clear(t); mpz_clear(s); mpz_clear(seckey);
    mpz_clear(Q.x); mpz_clear(Q.y);
    return number_str;
}

void Public_Key_X_Coordinate_To_Taproot_Tweaked_Pubkey() {// to be continued not complete
}

const char * Point_To_Bech32m_P2TR_Address() { // to be continued not complete
    char pub[] = "da4710964f7852695de2da025290e24af6d8c281de5a0b902b7135fd9fd74d21";
    gmp_snprintf(cpub, 68, "%s", pub);
    hexs2bin(cpub, bin_publickey);
    segwit_addr_encode(bech32_output, "bc", 1, bin_publickey, 32);
    return bech32_output;
}

const char * Private_Key_To_WIF(mpz_t pk, bool compressed) {
    size_t wif_size = 54;
    char privatekey[66];
    char bin_privatekey[39];
    bin_privatekey[0] = 0x80;
    gmp_snprintf(privatekey, 65, "%0.64Zx", pk);
    hexs2bin(privatekey, bin_privatekey+1);
    if (compressed) {
        bin_privatekey[33] = 0x01;
        sha256(bin_privatekey, 34, bin_sha256);
        sha256(bin_sha256, 32, bin_sha256);    
        bin_privatekey[34] = bin_sha256[0];
        bin_privatekey[35] = bin_sha256[1];
        bin_privatekey[36] = bin_sha256[2];
        bin_privatekey[37] = bin_sha256[3];
        b58enc(wif, &wif_size, bin_privatekey, 38);
        return wif;
    }
    else {
        sha256(bin_privatekey, 33, bin_sha256);
        sha256(bin_sha256, 32, bin_sha256);    
        bin_privatekey[33] = bin_sha256[0];
        bin_privatekey[34] = bin_sha256[1];
        bin_privatekey[35] = bin_sha256[2];
        bin_privatekey[36] = bin_sha256[3];
        b58enc(wif, &wif_size, bin_privatekey, 37);
        return wif;
    }
}

int main(int argc, char *argv[])
{
    mpz_init(EC.a); mpz_init(EC.b); mpz_init(EC.p); mpz_init(EC.n);
    mpz_set_str(EC.a, "0x0", 0); mpz_set_str(EC.b, "0x7", 0);
    mpz_set_str(EC.p, "0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 0);
    mpz_set_str(EC.n, "0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 0);
    
    mpz_init(Curve_G.x); mpz_init(Curve_G.y);
    mpz_set_str(Curve_G.x, "0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", 0);
    mpz_set_str(Curve_G.y, "0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 0);
    
    struct Point R;
    mpz_init(R.x); mpz_init(R.y);
    mpz_set(R.x, Curve_G.x); mpz_set(R.y, Curve_G.y);
    
    mpz_t m;
    mpz_init(m);
    
    //char pub[] = "0397c4e775d49f77c67f0ca9486d0694c8df1ab67e7d7fdb64d8413b79d8409f8c";
    puts("");
    struct timeval  tv1, tv2;
    gettimeofday(&tv1, NULL);
    
    mpz_set_str(m, "0x8b8f7d3fdd5f346b728775a60212427a0fd474607ce951c7774412fe9ca27685", 0);
    gmp_printf("Taproot Tweaked privkey: %064s\n", Taproot_Tweak_PrivKey(m));
    gmp_printf("WIF_U: %s\n", Private_Key_To_WIF(m, false)); // WIF Uncompressed
    gmp_printf("WIF_C: %s\n", Private_Key_To_WIF(m, true)); // WIF Compressed
    gmp_printf("Address_Bech32m_P2TR: %s\n", Point_To_Bech32m_P2TR_Address()); // Bech32 P2TR Address
    for (int i = 0; i < 1; i++) {
        gmp_printf("X:%0.64Zx Y:%0.64Zx\n", R.x, R.y);
        gmp_printf("Address_U: %s\n", Point_To_Legacy_Address(R, false)); // P2PKH Uncompressed Address
        gmp_printf("Address_C: %s\n", Point_To_Legacy_Address(R, true)); // P2PKH Compressed Address
        gmp_printf("Address_P2SH: %s\n", Point_To_P2SH_Address(R)); // P2SH Address
        gmp_printf("Address_Bech32_P2WPKH: %s\n", Point_To_Bech32_P2WPKH_Address(R)); // Bech32 P2WPKH Address
        gmp_printf("Address_Bech32_P2WSH: %s\n", Point_To_Bech32_P2WSH_Address(R)); // Bech32 P2WSH Address
        gmp_printf("Hash160_U: %s\n", Point_To_Hash160(R, false)); // Hash160 Uncompressed
        gmp_printf("Hash160_C: %s\n", Point_To_Hash160(R, true)); // Hash160 Compressed
        Point_Addition(R, Curve_G, &R);
        //if (strcmp(Point_To_Cpub(Q), pub) == 0) { printf("Point %s found\n", Point_To_Cpub(R)); }
    }

    puts("");
    gettimeofday(&tv2, NULL);
    printf ("Total time = %f seconds\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec));
    
    mpz_clear(EC.a); mpz_clear(EC.b); mpz_clear(EC.p); mpz_clear(EC.n); // free memory for mpz variables
    mpz_clear(Curve_G.x); mpz_clear(Curve_G.y);
    mpz_clear(R.x); mpz_clear(R.y);
    mpz_clear(m);
    return 0;
}
