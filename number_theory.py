import random

MILLER_RABIN_ROUNDS = 6

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def gcd(a, b):
    return gcd(b, a % b) if a % b else b


def mod_inv(a, n):
    x0, x1 = 1, 0
    while n != 0:
        q, a, n = a // n, n, a % n
        x0, x1 = x1, x0 - q * x1
    return x0

def prime_mod_inv(a, prime):
    return pow(a, prime-2, prime)

def sieve_era(B): # returns primes <= b
    l = [True] * (B+1)
    l[0] = l[1] = False

    for i in range(2,int(B**0.5)+1):
        if l[i]:
            for a in range(i*i,B+1,i):
                l[a] = False

    return [i for i,j in enumerate(l) if j]

def sieve_in_range(small_primes, start, end):
    sieve = [True] * (end-start)

    # Sieve out multiples of the prime
    for p in small_primes:
        closest_above = closest_multiple_above(p, start)

        for i in range(closest_above, end, p):
            if i >= start:
                sieve[i-start] = False

    return sieve

def closest_multiple_above(prime, target):
    return ((target +  prime - 1) // prime) * prime

# Miller-Rabin primality test
def is_prime(n):
    n = int(n)

    if n in [0,1,4,6,8,9]:
        return False

    if n in [2,3,5,7]:
        return True

    if n % 2 == 0:
        return False

    a = 0
    m = n - 1
    k = 0

    while m % 2 == 0:
        k += 1
        m = m // 2

    m = int(m)

    for i in range(MILLER_RABIN_ROUNDS):
        a = random.randrange(2, n-1)
        x = pow(a, m, n)

        if x == 1 or x == n - 1:
           continue

        for j in range(1, k):
            x = pow(x,2,n)

            if x == 1:
                return False
            if x == n - 1:
                break
        else: # only executed if break above not executed
            return False

    return True
