k0 = 1
k = 10
c = 100
freq = 1
numer = complex(0, k*c*freq)
denom = complex(k, c*freq)
comMod = k0 + numer/denom
print(numer)
print(denom)
print(comMod)
print(comMod.real)
print(comMod.imag)