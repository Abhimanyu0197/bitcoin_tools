from datetime import datetime
import secp256k1 as ice

G = bytes(bytearray.fromhex('0479be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8'))
dip = 10
scalar = 1171823477249417952904755003699827898
point = ice.scalar_multiplication(scalar).hex()
numb = 0
for i in range(dip):
    if scalar % 2 == 0: 
        numb = scalar // 2
    else:
        scalar = scalar - 1
        numb = scalar // 2
    scalar = numb
group = list()
group.append(point)
#final_group = list()
now = datetime.now()
time = now.strftime("%H:%M:%S")
print(f"[{time}] Starting out")
for i in range(dip):
    final_group = list()
    for p in group:
        point = bytes(bytearray.fromhex(p))
        final_group.append(ice.point_multiplication(57896044618658097711785492504343953926418782139537452191302581570759080747169, point).hex())
        point = bytes(bytearray.fromhex(ice.point_subtraction(point,G).hex()))        
        final_group.append(ice.point_multiplication(57896044618658097711785492504343953926418782139537452191302581570759080747169, point).hex())
    group = final_group
print(f'Number of points: {len(final_group)}')
for p in final_group:
    print(p)
index = 1
for p in final_group:
    if p == ice.scalar_multiplication(scalar).hex():
        print(f'Found: [{index}] {scalar} {p}')
    index += 1
now = datetime.now()
time = now.strftime("%H:%M:%S")
print(f"[{time}] Finished")
