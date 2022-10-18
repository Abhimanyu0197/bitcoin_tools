import hashlib
import base58
import binascii
from bit import Key 

first_encode = base58.b58decode("5KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKSmnqY")
private_key_full = binascii.hexlify(first_encode)
private_key = private_key_full[2:-8]
private_key_hex = private_key.decode("utf-8")
print(private_key_hex)
keyU = Key.from_hex(private_key_hex)
keyU._public_key = keyU._pk.public_key.format(compressed=False)
print(f"{keyU.to_wif()}")
#123456789ABCDEFGHJKLMNPQRSTUVWXYZabcdefghijkmnopqrstuvwxyz - Base58
#********************Pattern WIF********************    ********************Actual WIF*********************
#5J1FeexV6bAHb8ybZjqQMjJrcCrHGW9sb6uFAddressRichRich -> 5J1FeexV6bAHb8ybZjqQMjJrcCrHGW9sb6uFAddressRieKzgx2(18e6b71103e7e34378acc33b7b654d0c57c35246799eb6b07dd7305540c99a6f)
#5JingLeBeLLsJingLeBeLLsJingLeALLTheWayWhatAFunRidet -> 5JingLeBeLLsJingLeBeLLsJingLeALLTheWayWhatAFuk51fAz(77362453f9c9f4835fc1f1bf838d1807c1d5932c4cb883d2572e5fa4ba866d84)
#5JudgementDayJudgementDayJudgementDayJudgementDay11
#5JBadAddressDoNotSendAssetsYouMightLoseitALLFor1111
#5JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
#5KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
#5KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKSmnqY  

#5JQw5LhP2WPob4NBZZLf3ECp8zuXg9pvu3fCujK81Uk8raNPJDj - 4eab43dc61cdf3b298d8b6721d634f30ea1ff6060def27c3f8527532178aaeb1
#5JTimesChanceLLoronbrinkofsecondbaiLoutforbanks1111
#5JSatoshiSatoshiSatoshiSatoshiSatoshiSatoshiSatoshi
#5KFuckingBitcoinFuckingBitcoinFuckingBitcoinBitcoin
