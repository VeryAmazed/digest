from fractions import Fraction

# 1/k + 2/k + 3/k + ... + k/k
def lb0(k):
    assert k > 0
    return Fraction(k+1,2)

# 1/k + 2/(k+1) + 3/(k+2) + ... + k/(2k-1)
# linear with slope approaching 1-log(2)
def lb1(k):
    assert k > 0
    return sum(Fraction(i+1,k+i) for i in range(k))

# add
def lb2(k):
    a = lb1(k)


from itertools import count
# from math import log
for k in count(1,100):
    print(k, float(lb0(k) / lb1(k)))
