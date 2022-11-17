#include <stdio.h>
#include <stdlib.h>
#include <gmp.h> 

int main(int argc, char **argv)  {
    
    char *curve_N = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141";
    mpz_t A, B, C, secp256k1_N, n10, n16, arg_4;
    mpz_t inverse_multiplier;
    
    mpz_init_set_str(secp256k1_N, curve_N, 16);
    mpz_init_set_ui(A, 0);
    mpz_init_set_ui(B, 0);
    mpz_init_set_ui(C, 0);
    mpz_init_set_ui(n10, 10);
    mpz_init_set_ui(n16, 16);
    mpz_init_set_ui(arg_4, 0);	
    mpz_init(inverse_multiplier);
	
    if(argc < 5) {
	printf("Missing parameters\n");
	exit(0);
    }
		
    mpz_set_str(A, argv[1], 0);
    mpz_set_str(B, argv[3], 0);
    mpz_set_str(arg_4, argv[4], 0);
    
    switch(argv[2][0]) {
		case '+':
			mpz_add(C, A, B);	
			mpz_mod(C, C, secp256k1_N);
		break;
		case '-':
			mpz_sub(C, A, B);
			mpz_mod(C, C, secp256k1_N);
		break;
		case '/':
			mpz_invert(inverse_multiplier, B, secp256k1_N);
			mpz_mul(C, A, inverse_multiplier);
			mpz_mod(C, C, secp256k1_N);
		break;
		case 'x':
			mpz_mul(C, A, B);
			mpz_mod(C, C, secp256k1_N);
		break;		
	}
    
    if (mpz_cmp(arg_4, n10) == 0) gmp_printf("Result: %Zd\n", C);
    if (mpz_cmp(arg_4, n16) == 0) gmp_printf("Result: %Zx\n", C);
    
    mpz_clear(A);mpz_clear(B);mpz_clear(C); 
    mpz_clear(secp256k1_N);mpz_clear(n10);mpz_clear(n16);
    mpz_clear(inverse_multiplier);
    return 0;
}
