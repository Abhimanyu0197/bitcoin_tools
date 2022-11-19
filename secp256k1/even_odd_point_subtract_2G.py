from secp256k1 import *

curve_2G = scalar_multiplication(2)
curve_3G = scalar_multiplication(3)
Q = scalar_multiplication(8897)
if Q == curve.g:
    raise ValueError('secp256k1 default Generator point is odd (has scalar 1)')
while True:
    if Q == curve_2G:
        print(f'Point is even (has Q=k(even scalar)*G)')
        break
    if Q == curve_3G:
        print(f'Point is odd (has Q=k(odd scalar)*G')
        break
    Q = point_subtraction(Q, curve_2G)
