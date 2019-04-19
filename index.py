"""
Usage: python index.py beta p alpha
where beta is some number such that beta \equiv alpha^x (mod p), where alpha is a primitive root mod p
"""
import logging
import sys
import math
import random


def precompute_prime_base(B):
    """
    :param B: bound for our index calculus
    :return: list of primes below B

    Basically just the sieve of Erastothenes
    """
    l = [True] * (B + 1)
    l[0] = l[1] = False
    for i in range(2, int(B ** 0.5) + 1):
        if l[i]:
            for a in range(i * i, B + 1, i):
                l[a] = False
    return [i for i, j in enumerate(l) if j]


def prime_factorization(n, primes):
    ret = dict()
    for p in primes:
        ret[p] = 0
        while n % p == 0:
            n //= p
            ret[p] += 1
    return ret if n == 1 else None


def index_calculus(beta, p, alpha):
    """
    :param beta:
    :param p:
    :param alpha:
    :return:
    """
    factor_base = precompute_prime_base(find_bound(p))
    unused = set(factor_base)
    factor_products = []
    checked = set()
    while len(unused) > 0 or len(factor_products) < len(factor_base):
        if len(checked) == p - 1:
            print("The bound B wasn't big enough!")
            exit(0)
        exp = random.randint(1, p - 1)
        if exp in checked:
            continue
        checked.add(exp)
        prime_factors = prime_factorization(pow(alpha, exp, p), factor_base)
        if prime_factors:
            for p in prime_factors:
                unused.remove(p)
            factor_products.append(prime_factors)





def find_bound(n):
    return int(math.exp((2**(-0.5)*(math.log(2*n)*math.log(math.log(2*n)))**0.5)))



def main():
    if len(sys.argv) < 4:
        print("Usage: python index.py beta p alpha, where beta = alpha^x (mod p)")
        return 1

    logging.basicConfig(filename='index.log', filemode='w', level=logging.DEBUG, format='%(message)s')

    beta = int(sys.argv[1])
    p = int(sys.argv[2])
    alpha = int(sys.argv[3])
    print("The value of x is",  index_calculus(beta, p, alpha))
    return 0


if __name__ == "__main__":
   main()