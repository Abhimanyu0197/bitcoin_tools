from datetime import datetime
import secp256k1 as ice

G = bytes(bytearray.fromhex('0479be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8'))
point = '04ceb6cbbcdbdf5ef7150682150f4ce2c6f4807b349827dcdbdd1f2efa885a26302b195386bea3f5f002dc033b92cfc2c9e71b586302b09cfe535e1ff290b1b5ac' #120_puzzle_point
group = list()
group.append(point)
#final_group = list()
dip = 10
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
    with open("points_ice.txt", "a", encoding="utf-8") as f:
        f.write(f"{p} \n")
for p in final_group:
    print(p)
now = datetime.now()
time = now.strftime("%H:%M:%S")
print(f"[{time}] Finished")

