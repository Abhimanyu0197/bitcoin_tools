#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
//#include <time.h>
#include <sys/time.h>
#include <gmp.h>
#include "secp256k1/secp256k1.h"


int main(int argc, char *argv[])
{
    secp256k1_context_init();

    Affine_Point R, stride_Point;
    init_AffinePoint(&R); init_AffinePoint(&stride_Point);
    mpz_t m, stride;
    mpz_inits(m, stride, NULL);
    mpz_set_str(m, "0x1", 0);
    mpz_set_str(stride, "0x1", 0);
    Scalar_Multiplication(&R, m);
    Scalar_Multiplication(&stride_Point, stride);

    Jacobian_Point JP, JQ;
    init_set_JacobianPoint(&JP, "0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", "0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", "0x1");
    init_set_JacobianPoint(&JQ, "0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", "0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", "0x1");

    puts("");
    struct timeval  tv1, tv2;
    gettimeofday(&tv1, NULL);

    for (int i = 0; i < 10; i++) {
        gmp_printf("Private key: %0.64Zx\n", m);
        gmp_printf("WIF_U: %s\n", Private_Key_To_WIF(m, false)); // WIF Uncompressed
        gmp_printf("WIF_C: %s\n", Private_Key_To_WIF(m, true)); // WIF Compressed
        gmp_printf("X:%0.64Zx Y:%0.64Zx\n", R.x, R.y);
        gmp_printf("Hash160_U: %s\n", Point_To_Hash160(&R, false)); // Hash160 Uncompressed
        gmp_printf("Hash160_C: %s\n", Point_To_Hash160(&R, true)); // Hash160 Compressed
        gmp_printf("Address_U: %s\n", Point_To_Legacy_Address(&R, false)); // P2PKH Uncompressed Address
        gmp_printf("Address_C: %s\n", Point_To_Legacy_Address(&R, true)); // P2PKH Compressed Address
        gmp_printf("Address_P2SH: %s\n", Point_To_P2SH_Address(&R)); // P2SH Address
        gmp_printf("Address_Bech32_P2WPKH: %s\n", Point_To_Bech32_P2WPKH_Address(&R)); // Bech32 P2WPKH Address
        gmp_printf("Address_Bech32_P2WSH: %s\n", Point_To_Bech32_P2WSH_Address(&R)); // Bech32 P2WSH Address
        gmp_printf("Address_Bech32m_P2TR: %s\n", Point_To_Bech32m_P2TR_Address(&R)); // P2TR Taproot Address
        gmp_printf("Taproot Tweaked privkey: %064s\n", Taproot_Tweaked_PrivKey(m)); // Taproot tweaked privkey
        gmp_printf("Taproot Tweaked pubkey: %064s\n", Point_To_Taproot_Tweaked_Pubkey(&R)); // Taproot tweaked pubkey
        mpz_add(m, m, stride);
        //Point_Addition(&R, &R, &stride_Point);
        Point_Addition2(&R, &R, &stride_Point);
        //JacobianPoint_Addition(&JP, &JP, &JQ);
        //JacobianPoint_To_AffinePoint(&JP, &R);
        //if (strcmp(Point_To_Cpub(Q), pub) == 0) { printf("Point %s found\n", Point_To_Cpub(R)); }
        puts("");
    }

    puts("");
    gettimeofday(&tv2, NULL);
    printf ("Total time = %f seconds\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec));

    clear_AffinePoint(&R);
    clear_AffinePoint(&stride_Point);
    clear_JacobianPoint(&JP);
    clear_JacobianPoint(&JQ);
    mpz_clears(m, stride, NULL);

    secp256k1_context_clear();
    return 0;
}
