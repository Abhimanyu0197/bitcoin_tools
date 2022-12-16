#ifndef SECP256K1_H
#define SECP256K1_H

typedef struct { mpz_t a; mpz_t b; mpz_t p; mpz_t n; } Elliptic_Curve;
typedef struct { mpz_t x; mpz_t y; } Affine_Point;
typedef struct { mpz_t x; mpz_t y; mpz_t z; } Jacobian_Point;

void init_JacobianPoint(Jacobian_Point *J);
void clear_JacobianPoint(Jacobian_Point *J);
void init_AffinePoint(Affine_Point *A);
void clear_AffinePoint(Affine_Point *A);
void init_set_JacobianPoint(Jacobian_Point *J, char *x, char *y, char *z);
void init_set_AffinePoint(Affine_Point *A, char *x, char *y);
void print_JacobianPoint(const Jacobian_Point *J);
void print_AffinePoint(const Affine_Point *A);
void reduce_JacobianPoint(Jacobian_Point *J);
void JacobianPoint_To_AffinePoint(Jacobian_Point *J, Affine_Point *A);
void AffinePoint_To_JacobianPoint(const Affine_Point *A, Jacobian_Point *J);
void JacobianPoint_Doubling(Jacobian_Point *J, Jacobian_Point *J1);
void JacobianPoint_Addition(Jacobian_Point *J, Jacobian_Point *J1, Jacobian_Point *J2);
void Point_Doubling(Affine_Point *A, Affine_Point *A1);
void Point_Addition(Affine_Point *A, Affine_Point *A1, Affine_Point *A2);
void Point_Addition2(Affine_Point *A, Affine_Point *A1, Affine_Point *A2);
void Scalar_Multiplication(Affine_Point *A, const mpz_t m);
void Point_Multiplication(Affine_Point *A, Affine_Point *A1, const mpz_t m);
void Point_Negation(Affine_Point *A);
void Point_Subtraction(Affine_Point *A, Affine_Point *A1, Affine_Point *A2);
void Point_Division(Affine_Point *A, Affine_Point *A1, const mpz_t k);
bool Point_On_Curve(const Affine_Point *A);
const char * Point_To_Upub(const Affine_Point *A);
const char * Point_To_Cpub(const Affine_Point *A);
const char * Point_To_Hash160(const Affine_Point *pubkey, bool compressed);
const char * Point_To_Legacy_Address(const Affine_Point *pubkey, bool compressed);
const char * Point_To_P2SH_Address(const Affine_Point *pubkey);
const char * Point_To_Bech32_P2WPKH_Address(const Affine_Point *pubkey);
const char * Point_To_Bech32_P2WSH_Address(const Affine_Point *pubkey);
const char * Tagged_Hash(const mpz_t x);
const char * Taproot_Tweaked_PrivKey(const mpz_t k);
const char * Point_To_Taproot_Tweaked_Pubkey(const Affine_Point *pubkey);
const char * Point_To_Bech32m_P2TR_Address(const Affine_Point *pubkey);
const char * Private_Key_To_WIF(mpz_t pk, bool compressed);
void secp256k1_context_init();
void secp256k1_context_clear();

#endif
