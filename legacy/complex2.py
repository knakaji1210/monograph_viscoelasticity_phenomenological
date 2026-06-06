a = 3
b = 5

k = a + b*(2j/2)
print(k.real)
print(k.imag)


k0 = 1
k = 10
c = 100
f = 1
#numer = complex(0, k*c*f)
numer = k*c*f*(2j/2)
#denom = complex(k, c*freq)
denom = k + c*f*(2j/2)
comMod = k0 + numer/denom
print(numer)
print(denom)
print(comMod)
print(comMod.real)
print(comMod.imag)