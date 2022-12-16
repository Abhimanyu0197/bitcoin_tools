#include <gmp.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include "secp256k1.h"
#include "../util.h"
#include "../base58/libbase58.h"
#include "../rmd160/rmd160.h"
#include "../sha256/sha256.h"
#include "../bech32/segwit_addr.h"

Elliptic_Curve EC;
Affine_Point Curve_Generator, A1Temp, A2Temp, AffinePoint, TapPoint;
Jacobian_Point JacobianPoint, J1Temp, J2Temp;
mpz_t u1, v1, u, v, us2, vs2, vs3, usw2, vs2v2, _2vs2v2, a, vs3u2;
mpz_t z2, x2, _3x2, w, s, s2, b, _8b, _8y2s2, y2, h;
mpz_t slope, temp, tap_t, seckey, tap_s;
mpz_t X, Y;
char tap_tweak[] = "TapTweak";
static char upub[132];
static char cpub[68];
static char address[50];
static char hash160[42];
static char bech32_output[128];
static char wif[54];
static char bin_publickey[65];
static char bin_sha256[32];
static char bin_digest[60];
static char bin_rmd160[20];
static char privatekey[66];
static char bin_privatekey[39];
static char hex_mpz[66];
static char pub_bin[33];
static char taghash[33];
static char taghash_final[99];

void init_JacobianPoint(Jacobian_Point *J) { mpz_inits(J->x, J->y, J->z, NULL); return; }
void clear_JacobianPoint(Jacobian_Point *J) { mpz_clears(J->x, J->y, J->z, NULL); return; }
void init_AffinePoint(Affine_Point *A) { mpz_inits(A->x, A->y, NULL); return; }
void clear_AffinePoint(Affine_Point *A) { mpz_clears(A->x, A->y, NULL); return; }

void init_set_JacobianPoint(Jacobian_Point *J, char *x, char *y, char *z) {
    init_JacobianPoint(J);
    mpz_set_str(J->x, x, 0); mpz_set_str(J->y, y, 0); mpz_set_str(J->z, z, 0);
    return;
}

void init_set_AffinePoint(Affine_Point *A, char *x, char *y) {
    init_AffinePoint(A);
    mpz_set_str(A->x, x, 0); mpz_set_str(A->y, y, 0);
    return;
}

void print_JacobianPoint(const Jacobian_Point *J) {
    gmp_printf("X:%0.64Zx\nY:%0.64Zx\nZ:%0.64Zx\n", J->x, J->y, J->z);
    return;
}

void print_AffinePoint(const Affine_Point *A) {
    gmp_printf("X:%0.64Zx Y:%0.64Zx\n", A->x, A->y);
    return;
}

void reduce_JacobianPoint(Jacobian_Point *J) {
    mpz_invert(temp, J->z , EC.p);
    mpz_mul(J->x, J->x, temp);
    mpz_mod(J->x, J->x, EC.p);
    mpz_mul(J->y, J->y, temp);
    mpz_mod(J->y, J->y, EC.p);
    mpz_set_ui(J->z, 1);
    return;
}

void JacobianPoint_To_AffinePoint(Jacobian_Point *J, Affine_Point *A) {
    reduce_JacobianPoint(J);
    mpz_set(A->x, J->x); mpz_set(A->y, J->y);
    return;
}

void AffinePoint_To_JacobianPoint(const Affine_Point *A, Jacobian_Point *J) {
    mpz_set(J->x, A->x); mpz_set(J->y, A->y); mpz_set_ui(J->z, 1);
    return;
}

void JacobianPoint_Doubling(Jacobian_Point *J, Jacobian_Point *J1) {
    mpz_pow_ui(z2, J1->z, 2);
    mpz_mod(z2, z2, EC.p);
    mpz_set_ui(z2, 0);
    mpz_pow_ui(x2, J1->x, 2);
    mpz_mod(x2, x2, EC.p);
    mpz_add(_3x2, x2, x2);
    mpz_mod(_3x2, _3x2, EC.p);
    mpz_add(_3x2, _3x2, x2);
    mpz_mod(_3x2, _3x2, EC.p);
    mpz_add(w, z2, _3x2);
    mpz_mod(w, w, EC.p);
    mpz_mul(s, J1->y, J1->z);
    mpz_mod(s, s, EC.p);
    mpz_mul(b, J1->y, s);
    mpz_mod(b, b, EC.p);
    mpz_mul(b, b, J1->x);
    mpz_mod(b, b, EC.p);
    mpz_pow_ui(h, w, 2);
    mpz_mod(h, h, EC.p);
    mpz_add(_8b, b, b);
    mpz_mod(_8b, _8b, EC.p);
    mpz_add(_8b, _8b, _8b);
    mpz_mod(_8b, _8b, EC.p);
    mpz_add(_8b, _8b, _8b);
    mpz_mod(_8b, _8b, EC.p);
    mpz_sub(h, h, _8b);
    mpz_mod(h, h, EC.p);
    mpz_mul(J->x, h, s);
    mpz_mod(J->x, J->x, EC.p);
    mpz_add(J->x, J->x, J->x);
    mpz_mod(J->x, J->x, EC.p);
    mpz_pow_ui(s2, s, 2);
    mpz_mod(s2, s2, EC.p);
    mpz_pow_ui(y2, J1->y, 2);
    mpz_mod(y2, y2, EC.p);
    mpz_mul(_8y2s2, y2, s2);
    mpz_mod(_8y2s2, _8y2s2, EC.p);
    mpz_add(_8y2s2, _8y2s2, _8y2s2);
    mpz_mod(_8y2s2, _8y2s2, EC.p);
    mpz_add(_8y2s2, _8y2s2, _8y2s2);
    mpz_mod(_8y2s2, _8y2s2, EC.p);
    mpz_add(_8y2s2, _8y2s2, _8y2s2);
    mpz_mod(_8y2s2, _8y2s2, EC.p);
    mpz_add(J->y, b, b);
    mpz_mod(J->y, J->y, EC.p);
    mpz_add(J->y, J->y, J->y);
    mpz_mod(J->y, J->y, EC.p);
    mpz_sub(J->y, J->y, h);
    mpz_mod(J->y, J->y, EC.p);
    mpz_mul(J->y, J->y, w);
    mpz_mod(J->y, J->y, EC.p);
    mpz_sub(J->y, J->y, _8y2s2);
    mpz_mod(J->y, J->y, EC.p);
    mpz_mul(J->z, s2, s);
    mpz_mod(J->z, J->z, EC.p);
    mpz_add(J->z, J->z, J->z);
    mpz_mod(J->z, J->z, EC.p);
    mpz_add(J->z, J->z, J->z);
    mpz_mod(J->z, J->z, EC.p);
    mpz_add(J->z, J->z, J->z);
    mpz_mod(J->z, J->z, EC.p);
    return;
}

void JacobianPoint_Addition(Jacobian_Point *J, Jacobian_Point *J1, Jacobian_Point *J2) {
    if(mpz_cmp(J1->x, J2->x) == 0 && mpz_cmp(J1->y, J2->y) == 0 && mpz_cmp(J1->y, J2->y) == 0) {
        JacobianPoint_Doubling(J, J1);
        return;
    } else {
        mpz_set(J1Temp.x, J1->x); mpz_set(J1Temp.y, J1->y); mpz_set(J1Temp.z, J1->z);
        mpz_set(J2Temp.x, J2->x); mpz_set(J2Temp.y, J2->y); mpz_set(J2Temp.z, J2->z);
        mpz_mul(u1, J2->y, J1->z);
        mpz_mod(u1, u1, EC.p);
        mpz_mul(v1, J2->x, J1->z);
        mpz_mod(v1, v1, EC.p);
        mpz_sub(u, u1, J1->y);
        mpz_mod(u, u, EC.p);
        mpz_sub(v, v1, J1->x);
        mpz_mod(v, v, EC.p);
        mpz_pow_ui(us2, u, 2);
        mpz_mod(us2, us2, EC.p);
        mpz_pow_ui(vs2, v, 2);
        mpz_mod(vs2, vs2, EC.p);
        mpz_mul(vs3, vs2, v);
        mpz_mod(vs3, vs3, EC.p);
        mpz_mul(usw2, us2, J1->z);
        mpz_mod(usw2, usw2, EC.p);
        mpz_mul(vs2v2, vs2, J1->x);
        mpz_mod(vs2v2, vs2v2, EC.p);
        mpz_add(_2vs2v2, vs2v2, vs2v2);
        mpz_mod(_2vs2v2, _2vs2v2, EC.p);
        mpz_sub(a, usw2, vs3);
        mpz_mod(a, a, EC.p);
        mpz_sub(a, a, _2vs2v2);
        mpz_mod(a, a, EC.p);
        mpz_mul(J->x, v, a);
        mpz_mod(J->x, J->x, EC.p);
        mpz_mul(vs3u2, vs3, J1Temp.y);
        mpz_mod(vs3u2, vs3u2, EC.p);
        mpz_sub(J->y, vs2v2, a);
        mpz_mod(J->y, J->y, EC.p);
        mpz_mul(J->y, J->y, u);
        mpz_mod(J->y, J->y, EC.p);
        mpz_sub(J->y, J->y, vs3u2);
        mpz_mod(J->y, J->y, EC.p);
        mpz_mul(J->z, vs3, J1Temp.z);
        mpz_mod(J->z, J->z, EC.p);
        return;
    }
}

void Point_Doubling(Affine_Point *A, Affine_Point *A1) {
    mpz_set(A1Temp.x, A1->x); mpz_set(A1Temp.y, A1->y);
    mpz_mul_ui(temp, A1->y, 2);
    mpz_invert(temp, temp, EC.p);
    mpz_mul(slope, A1->x, A1->x);
    mpz_mul_ui(slope, slope, 3);
    mpz_add(slope, slope, EC.a);
    mpz_mul(slope, slope, temp);
    mpz_mod(slope, slope, EC.p);
    mpz_mul(A->x, slope, slope);
    mpz_sub(A->x, A->x, A1Temp.x);
    mpz_sub(A->x, A->x, A1Temp.x);
    mpz_mod(A->x, A->x, EC.p);
    mpz_sub(temp, A1Temp.x, A->x);
    mpz_mul(A->y, slope, temp);
    mpz_sub(A->y, A->y, A1Temp.y);
    mpz_mod(A->y, A->y, EC.p);
    return;
}

void Point_Addition(Affine_Point *A, Affine_Point *A1, Affine_Point *A2) {
    mpz_set(A1Temp.x, A1->x); mpz_set(A1Temp.y, A1->y);
    mpz_set(A2Temp.x, A2->x); mpz_set(A2Temp.y, A2->y);
    if(mpz_cmp_ui(A1->x, 0) == 0 && mpz_cmp_ui(A1->y, 0) == 0) {
        mpz_set(A->x, A2->x);
        mpz_set(A->y, A2->y);
        return;
    }
    if(mpz_cmp_ui(A2->x, 0) == 0 && mpz_cmp_ui(A2->y, 0) == 0) {
        mpz_set(A->x, A1->x);
        mpz_set(A->y, A1->y);
        return;
    }
    if(mpz_cmp(A1->x, A2->x) == 0 && mpz_cmp(A1->y, A2->y) != 0) {
        mpz_set_ui(A->x, 0);
        mpz_set_ui(A->y, 0);
        return;
    }
    if(mpz_cmp(A1->x, A2->x) == 0 && mpz_cmp(A1->y, A2->y) == 0) {
        Point_Doubling(A, A1);
        return;
    } else {
        mpz_sub(temp, A2->x, A1->x);
        mpz_mod(temp, temp, EC.p);
        mpz_invert(temp, temp, EC.p);
        mpz_sub(slope, A2->y, A1->y);
        mpz_mul(slope, slope, temp);
        mpz_mod(slope, slope, EC.p);
        mpz_mul(A->x, slope, slope);
        mpz_sub(A->x, A->x, A1Temp.x);
        mpz_sub(A->x, A->x, A2Temp.x);
        mpz_mod(A->x, A->x, EC.p);
        mpz_sub(temp, A1Temp.x, A->x);
        mpz_mul(A->y, slope, temp);
        mpz_sub(A->y, A->y, A1Temp.y);
        mpz_mod(A->y, A->y, EC.p);
        return;
    }
}

void Point_Addition2(Affine_Point *A, Affine_Point *A1, Affine_Point *A2) {
    mpz_set(A1Temp.x, A1->x); mpz_set(A1Temp.y, A1->y);
    mpz_set(A2Temp.x, A2->x); mpz_set(A2Temp.y, A2->y);
    if(mpz_cmp_ui(A1->x, 0) == 0 && mpz_cmp_ui(A1->y, 0) == 0) {
        mpz_set(A->x, A2->x);
        mpz_set(A->y, A2->y);
        return;
    }
    if(mpz_cmp_ui(A2->x, 0) == 0 && mpz_cmp_ui(A2->y, 0) == 0) {
        mpz_set(A->x, A1->x);
        mpz_set(A->y, A1->y);
        return;
    }
    if(mpz_cmp(A1->x, A2->x) == 0 && mpz_cmp(A1->y, A2->y) != 0) {
        mpz_set_ui(A->x, 0);
        mpz_set_ui(A->y, 0);
        return;
    }
    if(mpz_cmp(A1->x, A2->x) == 0 && mpz_cmp(A1->y, A2->y) == 0) {
        Point_Doubling(A, A1);
        return;
    } else {

        mpz_sub(u, A2->y, A1->y);
        mpz_sub(v, A2->x, A1->x);
        mpz_invert(v, v, EC.p);
        mpz_mul(slope, u, v);
        mpz_mod(slope, slope, EC.p);

        mpz_mul(A->x, slope, slope);
        mpz_sub(A->x, A->x, A1Temp.x);
        mpz_sub(A->x, A->x, A2Temp.x);
        mpz_mod(A->x, A->x, EC.p);

        mpz_sub(temp, A1Temp.x, A->x);
        mpz_mul(A->y, slope, temp);
        mpz_sub(A->y, A->y, A1Temp.y);
        mpz_mod(A->y, A->y, EC.p);
        return;
    }
}

void Scalar_Multiplication(Affine_Point *A, const mpz_t m) {
    int no_of_bits, loop;
    no_of_bits = mpz_sizeinbase(m, 2);
    mpz_set(A->x, Curve_Generator.x);
    mpz_set(A->y, Curve_Generator.y);
    for(loop = no_of_bits - 2; loop >= 0 ; loop--) {
        Point_Doubling(A, A);
        if(mpz_tstbit(m, loop)) {
            Point_Addition(A, A, &Curve_Generator);
        }
    }
    return;
}

void Point_Multiplication(Affine_Point *A, Affine_Point *A1, const mpz_t m) {
    int no_of_bits, loop;
    no_of_bits = mpz_sizeinbase(m, 2);
    mpz_set(A->x, A1->x);
    mpz_set(A->y, A1->y);
    for(loop = no_of_bits - 2; loop >= 0 ; loop--) {
        Point_Doubling(A, A);
        if(mpz_tstbit(m, loop)) {
            Point_Addition(A, A, A1);
        }
    }
    return;
}

void Point_Negation(Affine_Point *A) {
    mpz_sub(A->y, EC.p, A->y);
    return;
}

void Point_Subtraction(Affine_Point *A, Affine_Point *A1, Affine_Point *A2) {
    if (mpz_cmp(A1->x, A2->x) == 0 && mpz_cmp(A1->y, A2->y) == 0) {
        mpz_set_ui(A->x, 0);
        mpz_set_ui(A->y, 0);
        return;
    } else {
        mpz_set(AffinePoint.x, A2->x);
        mpz_set(AffinePoint.y, A2->y);
        Point_Negation(&AffinePoint);
        Point_Addition(A, A1, &AffinePoint);
        return;
   }
}

void Point_Division(Affine_Point *A, Affine_Point *A1, const mpz_t k) {
    mpz_invert(temp, k, EC.n);
    Point_Multiplication(A, A1, temp);
    return;
}

bool Point_On_Curve(const Affine_Point *A) {
    mpz_pow_ui(X, A->x, 3);
    mpz_add(X, X, EC.b);
    mpz_mod(X, X, EC.p);
    mpz_pow_ui(Y, A->y, 2);
    mpz_mod(Y, Y, EC.p);
    if (mpz_cmp(X, Y) == 0) {
        return true;
    } else {
        return false;
    }
}

const char * Point_To_Upub(const Affine_Point *A) {
    gmp_snprintf(upub, 132, "04%0.64Zx%0.64Zx", A->x, A->y);
    return upub;
}

const char * Point_To_Cpub(const Affine_Point *A) {
    if(mpz_tstbit(A->y, 0) == 0) {
        gmp_snprintf(cpub, 67, "02%0.64Zx", A->x);
    } else {
        gmp_snprintf(cpub, 67, "03%0.64Zx", A->x);
    }
    return cpub;
}

const char * Point_To_Hash160(const Affine_Point *pubkey, bool compressed) {
    if(compressed) {
        if(mpz_tstbit(pubkey->y, 0) == 0) {
            gmp_snprintf (cpub, 67, "02%0.64Zx", pubkey->x);
        } else {
            gmp_snprintf(cpub, 67, "03%0.64Zx", pubkey->x);
        }
        hexs2bin(cpub, bin_publickey);
        sha256(bin_publickey, 33, bin_sha256);
    } else {
        gmp_snprintf(upub, 132, "04%0.64Zx%0.64Zx", pubkey->x, pubkey->y);
        hexs2bin(upub, bin_publickey);
        sha256(bin_publickey, 65, bin_sha256);
    }
    RMD160Data((const unsigned char*)bin_sha256, 32, bin_rmd160);
    tohex_dst(bin_rmd160, 20, hash160);
    return hash160;
}

const char * Point_To_Legacy_Address(const Affine_Point *pubkey, bool compressed) {
    size_t pubaddress_size = 50;
    if(compressed) {
        if(mpz_tstbit(pubkey->y, 0) == 0) {
            gmp_snprintf(cpub, 68, "02%0.64Zx", pubkey->x);
        } else {
            gmp_snprintf(cpub, 68, "03%0.64Zx", pubkey->x);
        }
        hexs2bin(cpub, bin_publickey);
        sha256(bin_publickey, 33, bin_sha256);
    } else {
        gmp_snprintf(upub, 132, "04%0.64Zx%0.64Zx", pubkey->x, pubkey->y);
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

const char * Point_To_P2SH_Address(const Affine_Point *pubkey) {
    size_t pubaddress_size = 50;
    if(mpz_tstbit(pubkey->y, 0) == 0) {
        gmp_snprintf(cpub, 68, "02%0.64Zx", pubkey->x);
    }
    else {
        gmp_snprintf(cpub, 68, "03%0.64Zx", pubkey->x);
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

const char * Point_To_Bech32_P2WPKH_Address(const Affine_Point *pubkey) {
    if(mpz_tstbit(pubkey->y, 0) == 0) {
        gmp_snprintf(cpub, 68, "02%0.64Zx", pubkey->x);
    }
    else {
        gmp_snprintf(cpub, 68, "03%0.64Zx", pubkey->x);
    }
    hexs2bin(cpub, bin_publickey);
    sha256(bin_publickey, 33, bin_sha256);
    RMD160Data(bin_sha256, 32, bin_digest);
    segwit_addr_encode(bech32_output, "bc", 0, bin_digest, 20);
    return bech32_output;
}

const char * Point_To_Bech32_P2WSH_Address(const Affine_Point *pubkey) {
    if(mpz_tstbit(pubkey->y, 0) == 0) {
        gmp_snprintf(cpub, 68, "02%0.64Zx", pubkey->x);
    }
    else {
        gmp_snprintf(cpub, 68, "03%0.64Zx", pubkey->x);
    }
    hexs2bin(cpub, bin_publickey);
    sha256(bin_publickey, 33, bin_sha256);
    segwit_addr_encode(bech32_output, "bc", 0, bin_sha256, 32);
    return bech32_output;
}

const char * Tagged_Hash(const mpz_t x) {
    gmp_snprintf(hex_mpz, 65, "%0.64Zx", x);
    hexs2bin(hex_mpz, pub_bin);
    sha256(tap_tweak, strlen(tap_tweak), taghash);
    memmove(taghash_final, taghash, sizeof(taghash));
    memmove(taghash_final + 32, taghash, sizeof(taghash));
    memmove(taghash_final + 64, pub_bin, sizeof(pub_bin));
    sha256(taghash_final, 96, bin_sha256);
    return tohex(bin_sha256, 32);
}

const char * Taproot_Tweaked_PrivKey(const mpz_t k) {
    Scalar_Multiplication(&AffinePoint, k);
    if(mpz_tstbit(AffinePoint.y, 0) == 0) {
        mpz_set(seckey, k);
    }
    else {
        mpz_sub(tap_s, EC.n, k);
        mpz_set(seckey, tap_s);
    }
    mpz_set_str(tap_t, Tagged_Hash(AffinePoint.x), 16);
    mpz_add(seckey, seckey, tap_t);
    mpz_mod(seckey, seckey, EC.n);
    return mpz_get_str(NULL, 16, seckey);
}

const char * Point_To_Taproot_Tweaked_Pubkey(const Affine_Point *pubkey) {
    mpz_set(TapPoint.x, pubkey->x);
    mpz_set(TapPoint.y, pubkey->y);
    if(mpz_tstbit(TapPoint.y, 0) != 0) {
        Point_Negation(&TapPoint);
    }
    mpz_set_str(tap_t, Tagged_Hash(TapPoint.x), 16);
    Scalar_Multiplication(&AffinePoint, tap_t);
    Point_Addition(&TapPoint, &AffinePoint, &TapPoint);
    return mpz_get_str(NULL, 16, TapPoint.x);
}

const char * Point_To_Bech32m_P2TR_Address(const Affine_Point *pubkey) {
    gmp_snprintf(cpub, 68, "%064s", Point_To_Taproot_Tweaked_Pubkey(pubkey));
    hexs2bin(cpub, bin_publickey);
    segwit_addr_encode(bech32_output, "bc", 1, bin_publickey, 32);
    return bech32_output;
}

const char * Private_Key_To_WIF(mpz_t pk, bool compressed) {
    size_t wif_size = 54;
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

void secp256k1_context_init() {
    mpz_inits(EC.a, EC.b, EC.p, EC.n, NULL);
    mpz_set_str(EC.a, "0x0", 0); mpz_set_str(EC.b, "0x7", 0);
    mpz_set_str(EC.p, "0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 0);
    mpz_set_str(EC.n, "0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141", 0);
    init_set_AffinePoint(&Curve_Generator, "0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798", "0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8");
    init_AffinePoint(&A1Temp);
    init_AffinePoint(&A2Temp);
    init_AffinePoint(&AffinePoint);
    init_AffinePoint(&TapPoint);
    init_JacobianPoint(&JacobianPoint);
    init_JacobianPoint(&J1Temp);
    init_JacobianPoint(&J2Temp);
    mpz_inits(u1, v1, u, v, us2, vs2, vs3, usw2, vs2v2, _2vs2v2, a, vs3u2, NULL);
    mpz_inits(z2, x2, _3x2, w, s, s2, b, _8b, _8y2s2, y2, h, NULL);
    mpz_inits(temp, slope, tap_t, seckey, tap_s, NULL);
    mpz_inits(X, Y, NULL);
    return;
}

void secp256k1_context_clear() {
    mpz_clears(EC.a, EC.b, EC.p, EC.n, NULL);
    clear_AffinePoint(&Curve_Generator);
    clear_AffinePoint(&A1Temp);
    clear_AffinePoint(&A2Temp);
    clear_AffinePoint(&AffinePoint);
    clear_AffinePoint(&TapPoint);
    clear_JacobianPoint(&JacobianPoint);
    clear_JacobianPoint(&J1Temp);
    clear_JacobianPoint(&J2Temp);
    mpz_clears(u1, v1, u, v, us2, vs2, vs3, usw2, vs2v2, _2vs2v2, a, vs3u2, NULL);
    mpz_clears(z2, x2, _3x2, w, s, s2, b, _8b, _8y2s2, y2, h, NULL);
    mpz_clears(temp, slope, tap_t, seckey, tap_s, NULL);
    mpz_clears(X, Y, NULL);
    return;
}
