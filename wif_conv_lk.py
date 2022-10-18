import hashlib
import base58
import binascii
from bit import Key 

wif = 'L21yJ6n3oRmeo2oQyAurBCsusiA4E3iFk3JYDUF1Mf6jtsWkYHVL'
first_encode = base58.b58decode(wif)
print(first_encode)
private_key_full = binascii.hexlify(first_encode)
private_key = private_key_full[2:-8]
print(private_key_full)
private_key_hex = private_key.decode("utf-8")
print(private_key_hex)
keyU = Key.from_hex(private_key_hex[0:64])
print(f"{keyU.to_wif()}")
#123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz - Base58
#********************Pattern WIF*********************   ********************Actual WIF*********************
