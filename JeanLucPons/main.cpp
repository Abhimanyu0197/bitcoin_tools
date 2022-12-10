#include <stdio.h>
#include <string>
#include <sys/time.h>
#include "secp256k1/SECP256k1.h"
using namespace std;

int main() {

    Secp256K1 *secp256k1 = new Secp256K1();
    secp256k1->Init();
    //Int privKey;
    //privKey.SetBase16("0000000000000000000000000000000000000000000000000000000000000001");
    //Point pub = secp256k1->ComputePublicKey(&privKey);
    //string calcAddress = secp->GetAddress(1, false, pub);
    //printf("Pub : %s \n", secp256k1->GetPublicKeyHex(false, pub));
    //Point S = (secp256k1->G);
    //Point Q;
    Point P = secp256k1->Double(S);
    //G2.Reduce();
    //printf("Pub : %s \n", secp256k1->GetPublicKeyHex(false, G2));
    //Point G3 = secp256k1->Add(G2, S);
    //G3.Reduce();
    //printf("Pub : %s \n", secp256k1->GetPublicKeyHex(true, G3));
    struct timeval  tv1, tv2;
    gettimeofday(&tv1, NULL);
    for (int i = 0; i < 1000000; i++) {
        //P.Reduce();
        //printf("Pub : %s \n", secp256k1->GetPublicKeyHex(true, P));
        //P = secp256k1->Add(P, secp256k1->G); // a little slower than Add2
        P = secp256k1->Add2(P, secp256k1->G); // fastest point addition
        //P = secp256k1->AddDirect(P, secp256k1->G); // slowest
    }
    puts("");
    gettimeofday(&tv2, NULL);
    printf ("Total time = %f seconds\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 + (double) (tv2.tv_sec - tv1.tv_sec));
    return 0;

}
