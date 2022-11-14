from secp256k1 import *

Q = [89565891926547004231252920425935692360644145829622209833684329913297188986597,
     12158399299693830322967808612713398636155367887041628176798871954788371653930]
Q = scalar_multiplication(2)
if Q == curve.g:
    raise ValueError('secp256k1 default Generator point is odd (has scalar 1)')
A = point_subtraction(Q, curve.g)
B = curve.g
prev_A = prev_B = [0, 0]
while True:
    if A == B:
        print(f'Point is even (has Q=k(even scalar)*G)')
        break
    if A == prev_B and B == prev_A:
        print(f'Point is odd (has Q=k(odd scalar)*G')
        break
    prev_A = A
    prev_B = B
    A = point_subtraction(A, curve.g)
    B = point_addition(B, curve.g)
