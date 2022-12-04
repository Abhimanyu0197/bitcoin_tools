#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <gmp.h>

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
    mpz_init(temp);
    mpz_init(slope);
    	
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
    
    mpz_clear(temp);
    mpz_clear(slope);
}

void Point_Addition(struct Point P, struct Point Q, struct Point *R) {
	
    if(mpz_cmp_ui(P.x, 0) == 0 && mpz_cmp_ui(P.y, 0) == 0) {
	mpz_set(R->x, Q.x);
	mpz_set(R->y, Q.y);
	return;
    }
    
    if(mpz_cmp_ui(Q.x, 0) == 0 && mpz_cmp_ui(Q.y, 0) == 0) {
	mpz_set(R->x, P.x);
	mpz_set(R->y, P.y);
	return;
    }
    
    mpz_t temp;
    mpz_init(temp);
    
    if(mpz_cmp_ui(Q.y, 0) != 0) { 
	mpz_sub(temp, EC.p, Q.y);
	mpz_mod(temp, temp, EC.p);
    } else
	mpz_set_ui(temp, 0);
        
    if(mpz_cmp(P.y, temp) == 0 && mpz_cmp(P.x, Q.x) == 0) {
	mpz_set_ui(R->x, 0);
	mpz_set_ui(R->y, 0);
	mpz_clear(temp);
	return;
    }
	
    if(mpz_cmp(P.x, Q.x) == 0 && mpz_cmp(P.y, Q.y) == 0) {
	Point_Doubling(P, R);		
	mpz_clear(temp);
	return;		
    } else {
	mpz_t slope;
	mpz_init_set_ui(slope, 0);
        
	mpz_sub(temp, P.x, Q.x);
	mpz_mod(temp, temp, EC.p);
        
	mpz_invert(temp, temp, EC.p);
	mpz_sub(slope, P.y, Q.y);
	mpz_mul(slope, slope, temp);
	mpz_mod(slope, slope, EC.p);
        
	mpz_mul(R->x, slope, slope);
	mpz_sub(R->x, R->x, P.x);
	mpz_sub(R->x, R->x, Q.x);
	mpz_mod(R->x, R->x, EC.p);
        
	mpz_sub(temp, P.x, R->x);
	mpz_mul(R->y, slope, temp);
	mpz_sub(R->y, R->y, P.y);
	mpz_mod(R->y, R->y, EC.p);
        		
	mpz_clear(temp);
	mpz_clear(slope);
	return;
    }
}

void Scalar_Multiplication(struct Point *R, mpz_t m) {
	
    struct Point P;
    mpz_init(P.x); mpz_init(P.y);
    mpz_set(P.x, Curve_G.x); mpz_set(P.y, Curve_G.y);
    struct Point Q, T;
    mpz_init(Q.x); mpz_init(Q.y);
    mpz_init(T.x); mpz_init(T.y);
    long no_of_bits, loop;	
    no_of_bits = mpz_sizeinbase(m, 2);
    mpz_set_ui(R->x, 0);
    mpz_set_ui(R->y, 0);
    if(mpz_cmp_ui(m, 0) == 0)
	return;
        		
    mpz_set(Q.x, P.x);
    mpz_set(Q.y, P.y);
    
    if(mpz_tstbit(m, 0) == 1){
	mpz_set(R->x, P.x);
	mpz_set(R->y, P.y);
    }

    for(loop = 1; loop < no_of_bits; loop++) {
	mpz_set_ui(T.x, 0);
	mpz_set_ui(T.y, 0);
	Point_Doubling(Q, &T);
	mpz_set(Q.x, T.x);
	mpz_set(Q.y, T.y);
	mpz_set(T.x, R->x);
	mpz_set(T.y, R->y);
	if(mpz_tstbit(m, loop))
	    Point_Addition(T, Q, R);
    }
    
    mpz_clear(P.x); mpz_clear(P.y);
    mpz_clear(Q.x); mpz_clear(Q.y);
    mpz_clear(T.x); mpz_clear(T.y);    
}

void Point_Multiplication(struct Point P, struct Point *R, mpz_t m) {
	
    struct Point Q, T;
    mpz_init(Q.x); mpz_init(Q.y);
    mpz_init(T.x); mpz_init(T.y);
    long no_of_bits, loop;	
    no_of_bits = mpz_sizeinbase(m, 2);
    mpz_set_ui(R->x, 0);
    mpz_set_ui(R->y, 0);
    
    if(mpz_cmp_ui(m, 0) == 0)
	return;
		
    mpz_set(Q.x, P.x);
    mpz_set(Q.y, P.y);
    
    if(mpz_tstbit(m, 0) == 1){
	mpz_set(R->x, P.x);
	mpz_set(R->y, P.y);
    }

    for(loop = 1; loop < no_of_bits; loop++) {
	mpz_set_ui(T.x, 0);
	mpz_set_ui(T.y, 0);
	Point_Doubling(Q, &T);
	mpz_set(Q.x, T.x);
	mpz_set(Q.y, T.y);
	mpz_set(T.x, R->x);
	mpz_set(T.y, R->y);
	if(mpz_tstbit(m, loop))
	    Point_Addition(T, Q, R);
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
    }
    else {
        struct Point Q;
        Point_Negation(B);
        mpz_init(Q.x);
        mpz_init(Q.y);        
        mpz_set(Q.x, B->x);
	mpz_set(Q.y, B->y);
        Point_Addition(A, Q, R);
        mpz_clear(Q.x); mpz_clear(Q.y);
        return;      
   }
}

void Point_Division(struct Point D, struct Point *S, mpz_t k) {
    mpz_t m;
    mpz_init(m);
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
    if (mpz_cmp(X, Y) == 0) { return true; }
    else { return false; }
    
}

static char upub[132];
static char cpub[68];

const char * Point_To_Upub(struct Point A) {
    
    gmp_snprintf(upub, 132, "04%0.64Zx%0.64Zx", A.x, A.y);
    return upub;
    
}

const char * Point_To_Cpub(struct Point A) {
    
    if(mpz_tstbit(A.y, 0) == 0) { gmp_snprintf(cpub, 68, "02%0.64Zx", A.x); }
    else { gmp_snprintf(cpub, 68,"03%0.64Zx", A.x); }
    return cpub;
    
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
    
    struct Point R, Q;
    mpz_init_set_ui(R.x, 0); mpz_init_set_ui(R.y, 0);
    mpz_init_set_ui(Q.x, 0); mpz_init_set_ui(Q.y, 0);
	
    mpz_t m;
    mpz_init(m);
    
    puts("");
    printf("Point_Doubling secp256k1 generator point G -> 2G\n");
    Point_Doubling(Curve_G, &R); // doubling secp256k1 generator point G
    gmp_printf("X:%0.64Zx Y:%0.64Zx\n", R.x, R.y); // result -> 2G
    printf("Point_On_Curve: %s\n", Point_On_Curve(R)? "True" : "False");
    printf("%s\n", Point_To_Upub(R));
    printf("%s\n", Point_To_Cpub(R));	
    puts("");
    
    printf("Scalar_Multiplication -> 5*G = 5G\n");
    mpz_set_str(m, "0x5", 0);
    Scalar_Multiplication(&R, m); // scalar_multiplication -> 5G = 5*G 
    gmp_printf("X:%0.64Zx Y:%0.64Zx\n", R.x, R.y); // result -> 5G
    printf("Point_On_Curve: %s\n", Point_On_Curve(R)? "True" : "False");
    puts("");
    
    printf("Point_Multiplication -> 2*5G = 10G\n");
    mpz_set_str(m, "0x2", 0);
    mpz_set(Q.x, R.x); mpz_set(Q.y, R.y);    
    Point_Multiplication(Q, &R, m);// point_multiplication -> 10G = 2*5G
    gmp_printf("X:%0.64Zx Y:%0.64Zx\n", R.x, R.y);// result -> 10G
    printf("Point_On_Curve: %s\n", Point_On_Curve(R)? "True" : "False");
    puts("");
    
    printf("Point_Division (10G / 2) -> multiplicative_inverse_of(2)*10G = 5G\n");
    mpz_set(Q.x, R.x); mpz_set(Q.y, R.y);
    Point_Division(Q, &R, m); // point_division (10G / 2) -> 5G = multiplicative_inverse(2)*10G
    gmp_printf("X:%Zx Y:%Zx\n", R.x, R.y);// result -> 5G
    printf("Point_On_Curve: %s\n", Point_On_Curve(R)? "True" : "False");
    puts("");
    
    printf("Point_Negation(5G) = additive_inverse_of(5G) 5G+additive_inverse_of(5G) = P(0,0) = Point_At_infinity\n");
    Point_Negation(&R); //point_negation(5G) = additive_inverse_of(5G) 5G+additive_inverse_of(5G) = P(0,0) = Point_At_infinity   
    gmp_printf("X:%0.64Zx Y:%0.64Zx\n", R.x, R.y);// result -> 115792089237316195423570985008687907852837564279074904382605163141518161494332G
    printf("Point_On_Curve: %s\n", Point_On_Curve(R)? "True" : "False");
    puts("");
    
    printf("Point_Addition 115792089237316195423570985008687907852837564279074904382605163141518161494332G + 1G\n");
    mpz_set(Q.x, R.x); mpz_set(Q.y, R.y);
    Point_Addition(Q, Curve_G, &R); // point_addition 115792089237316195423570985008687907852837564279074904382605163141518161494332G + 1G
    gmp_printf("X:%0.64Zx Y:%0.64Zx\n", R.x, R.y);// result -> 115792089237316195423570985008687907852837564279074904382605163141518161494333G
    printf("Point_On_Curve: %s\n", Point_On_Curve(R)? "True" : "False");
    puts("");
    
    printf("Point_Subtraction 1G - 1G = P(0,0) = Point_At_infinity\n");
    Point_Subtraction(Curve_G, &Curve_G, &R);// point_subtraction 1G - 1G
    gmp_printf("X:%0.64Zx Y:%0.64Zx\n", R.x, R.y); // result -> P(0,0) = Point_At_infinity
    printf("Point_On_Curve: %s\n", Point_On_Curve(R)? "True" : "False");
    
    mpz_clear(EC.a); mpz_clear(EC.b); mpz_clear(EC.p); mpz_clear(EC.n);// free memory for mpz variables
    mpz_clear(Curve_G.x); mpz_clear(Curve_G.y);
    mpz_clear(R.x); mpz_clear(R.y);
    mpz_clear(Q.x); mpz_clear(Q.y);
    mpz_clear(m);	
    
    return 0;
}
