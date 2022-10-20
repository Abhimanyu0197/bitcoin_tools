from ecc import FieldElement, Point

prime = 115792089237316195423570985008687907853269984665640564039457584007908834671663
a,b = FieldElement(0, prime), FieldElement(7, prime)
x = FieldElement(55066263022277343669578718895168534326250603453777594175500187360389116729240, prime)
y = FieldElement(32670510020758816978083085130507043184471273380659243275938904335757337482424, prime)
G = Point(x, y, a, b)
x = FieldElement(21505829891763648114329055987619236494102133314575206970830385799158076338148, prime)
y = FieldElement(98003708678762621233683240503080860129026887322874138805529884920309963580118, prime)
G5 = Point(x, y, a, b)
x = FieldElement(86918276961810349294276103416548851884759982251107, prime)
y = FieldElement(28597260016173315074988046521176122746119865902901063272803125467328307387891, prime)
Gmid = Point(x, y, a, b)
#print(G + G)
#print(115792089237316195423570985008687907852837564279074904382605163141518161494336 * G)
#print((115792089237316195423570985008687907852837564279074904382605163141518161494336 * G)-(4*G))
x = FieldElement(72488970228380509287422715226575535698893157273063074627791787432852706183111, prime)
y = FieldElement(62070622898698443831883535403436258712770888294397026493185421712108624767191, prime)
G10 = Point(x, y, a, b)
result = 4*G10 / 2
#print(result)
for s in range(1,11):
    result = s*G
    print('{}*G=({},{})'.format(s,result.x.num,result.y.num))
print()
print(G10 - G10)
print(0*G)