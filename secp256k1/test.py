import secp256k1

P  = secp256k1.scalar_multiplication(2)
print(secp256k1.point_to_upub(P))
print(secp256k1.point_to_cpub(P))    
G4 = secp256k1.point_multiplication(P, 2)
G6 = secp256k1.point_addition(P, G4)
G2 = secp256k1.point_subtraction(G6, G4)
G2 = secp256k1.point_division(G6, 3)
G12 = secp256k1.point_addition(G6, G6)
print(secp256k1.on_curve(G12))
print(secp256k1.point_negation(G6))

