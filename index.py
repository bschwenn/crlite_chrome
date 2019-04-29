"""
Usage: python index.py beta p alpha
where beta is some number such that beta \equiv alpha^x (mod p), where alpha is a primitive root mod p
"""
import sys
import math
from sympy.ntheory import factorint
import random
import logging
import flint
import numpy as np
from collections import defaultdict
from number_theory import *
from functools import reduce


def prime_factorization(n, primes):
    ret = dict()
    nn = n
    for p in primes:
        ret[p] = 0
        while n % p == 0:
            n //= p
            ret[p] += 1
    result =  [ret[i] for i in ret] if n == 1 else None
    return result


def rref_crt(mat, p, fz, factor_base, g):
    congruences = [[] for i in range(len(factor_base))]

    logging.debug("Reducing equations in prime subgroups...")
    for f,exp in fz.items():
        if exp > 0:
            reduction = rref_prime(np.copy(mat), f)

            for i,r in enumerate(reduction):
                congruences[i].append(r[-1])

    sols = {}
    n = [k for k in fz.keys()]

    logging.debug("Applying Chinese remainder theorem...")
    for i,a in enumerate(congruences):
        sols[i] = chinese_remainder(n, [int(i) for i in a])

    transl = {k:idx for idx,k in enumerate(fz)}

    modulus = reduce(lambda a,b: a*b, n)

    logging.debug("Fixing possible errors in factor base logarithms...")
    for i,s in sols.items():
        if not confirm_log(g, p, s, factor_base[i]):
            sols[i] = fix_prime_power_log(g, p, s, factor_base[i], fz, congruences[i], modulus, transl)

    return sols


def fix_prime_power_log(g, p, s, base_factor, fz, congruences, modulus, index_map):
    possibilities = [s]

    for f,exp in fz.items():
        if exp > 1:
            mod_f = int(modulus/f)
            new_possibilities = []
            power = f**exp

            for pos in possibilities:
                for t in range(0, f**exp, f):
                    r = chinese_remainder([mod_f, f**exp], [pos % mod_f, t + int(congruences[index_map[f]])])
                    new_possibilities.append(r)

            possibilities = new_possibilities

    for t in possibilities:
        if confirm_log(g, p, t, base_factor):
            return t

    return s

def rref_prime(mat, p):
    res, _ = flint.nmod_mat(mat.tolist(), p).rref()
    return res.table()


def solve_calculus(beta, p, alpha, prime_discrete_logs, factor_base):
    ret = None
    # Only use the prime in the factor base for which we got a valid log
    # computation. This might not be the case in particular if p-1
    # contains prime powers in its factorization.
    fb =  [f for i,f in enumerate(factor_base) if prime_discrete_logs[i]]
    prime_discrete_logs = [f for i,f in prime_discrete_logs.items() if f]

    logging.debug("Solving calculus with {} small discrete logs solved...".format(len(prime_discrete_logs)))
    for exp in range(1, p-1):
        calc = (beta*pow(alpha, exp, p)) % p
        prime_factors = prime_factorization(calc, fb)

        if prime_factors:
            ret = (-exp) % (p-1)

            for idx, ele in enumerate(prime_factors):
                if ele > 0:
                    ret = (ret + (ele * prime_discrete_logs[idx])) % (p-1)

            return ret

    return ret

def find_maximally_ind_set(mat, p, fz):
    all_i = set(range(len(mat[0])))
    ind = set(range(len(mat[0])))

    logging.debug("Finding maximally independent set...")
    for f,exp in fz.items():
        if exp > 0:
            ind_f = set()
            reduced = rref_prime(np.copy(mat).transpose(), f)
            r, c = 0, 0

            while r < len(reduced) and c < len(reduced[r]):
                if reduced[r][c].__int__() == 1:
                    ind_f.add(c)
                    r += 1
                    c += 1
                else: c+= 1

            not_ind = all_i.difference(ind_f)
            ind = ind.difference(not_ind)

    return ind

def index_calculus(beta, p, alpha):
    """
    :param beta:
    :param p:
    :param alpha:
    :return:
    """
    B = find_bound(p)
    factor_base = sieve_era(B)
    unused = set(factor_base)
    factor_products = []
    checked = set()
    order_fz = factorint(p-1)

    logging.debug("Finding dependencies; need {}...".format(len(factor_base)))
    while len(checked) < p-1:
        while len(unused) > 0 or len(factor_products) < len(factor_base):
            exp = random.randint(2, p - 1)

            if exp in checked:
                continue

            checked.add(exp)
            power = pow(alpha, exp, p)

            if power == 1:
                continue

            prime_factors = prime_factorization(power, factor_base)

            if prime_factors:
                prime_factors.append(exp)

                for index, i in enumerate(prime_factors[:-1]):
                    if factor_base[index] in unused:
                        unused.remove(factor_base[index])

                factor_products.append(prime_factors)

        mat = np.array(factor_products)
        independent = find_maximally_ind_set(mat, p, order_fz)

        if len(independent) == len(factor_base):
            logging.debug("Found enough factorizations!")
            factor_products = [e for i,e in enumerate(factor_products) if i in independent]
            prime_discrete_logs = rref_crt(np.array(factor_products), p, order_fz, factor_base, alpha)
            prime_discrete_logs = check_crt_results(prime_discrete_logs, factor_base, alpha, p)
            ret = solve_calculus(beta, p, alpha, prime_discrete_logs, factor_base)

            if ret:
                return ret
            else:
                print("We were unsuccessful. Make sure that your base is a primitive root and try again!")
                exit(0)
        else:
            logging.debug("Not enough independent factorizations! Need {}, had {}.".format(len(factor_base), len(independent)))
            factor_products = [e for i,e in enumerate(factor_products) if i in independent]


def check_crt_results(discrete_log_map, factor_base, g, p):
    logging.debug("Checking CRT results...")
    successful = []
    new_map = dict()

    for k,v in discrete_log_map.items():
        new_map[k] = confirm_log(g, p, v, factor_base[k])

    return new_map

# Returns the log if it is the log, else none
def confirm_log(g, p, x, f):
    t = pow(g,x,p)
    neg_f = -f % p

    # If g is truely a primitive root, then if we end up getting x
    # with g^x = -s (mod p) for some prime s in the factor base,
    # then we know that s = g^(x - (p-1)/2) (mod p). This is true
    # because p-1 = g^((p-1)/2) (mod p) when g is a primitive root,
    # since -1 = p - 1 (mod p) and g^(p-1) = 1 (mod p).
    logging.debug("{} =?= {} =?= {}".format(f, t, neg_f))

    if t == f:
        return x
    elif t == neg_f:
        return (x - int((p-1)/2)) % (p-1)
    else:
        return None

# Optimal bounds as described in https://doi.org/10.1016/j.mcm.2011.02.022
def find_bound(p):
    return int(3.33*L_p(p, 0.5, 0.476))


def L_p(p,s,c):
    return math.exp(c*(math.log(p)**s)*(math.log(math.log(p)))**s)


def main():
    if len(sys.argv) < 4:
        print("Usage: python index.py beta p alpha, where beta = alpha^x (mod p)")
        return 1

    logging.basicConfig(filename='index.log', filemode='w', level=logging.DEBUG)

    beta = int(sys.argv[1])
    p = int(sys.argv[2])
    alpha = int(sys.argv[3])
    x = index_calculus(beta, p, alpha)
    print("The value of x is",  x)
    print("Verification:", beta, "==", pow(alpha, x, p))
    return 0


if __name__ == "__main__":
    main()
