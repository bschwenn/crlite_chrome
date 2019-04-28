"""
Usage: python index.py beta p alpha
where beta is some number such that beta \equiv alpha^x (mod p), where alpha is a primitive root mod p
"""
import sys
import math
from sympy import Matrix
import random
import flint
import numpy as np
from collections import defaultdict
from number_theory import *

EXCESS_EQUATIONS = 0

def precompute_prime_base(B):
    """
    :param B: bound for our index calculus
    :return: list of primes below B

    Basically just the sieve of Erastothenes
    """
    return sieve_era(B)


def prime_factorization(n, primes):
    #ret = []

    #for p in primes:
    #    k = 1
    #    while n % p**k == 0:
    #        k += 1

    #    ret.append(k-1)

    ##print(n)
    ##print(primes)
    ##print(ret)
    #prod = 1
    #for i,p in enumerate(primes):
    #    prod *= primes[i] ** ret[i]
    #return ret if prod == n else None
    ret = dict()
    nn = n
    for p in primes:
        ret[p] = 0
        while n % p == 0:
            n //= p
            ret[p] += 1
    result =  [ret[i] for i in ret] if n == 1 else None
    return result


def find_row_to_elim_nonzero(mat, r_idx, c_idx, n, p):
    pivot_val = mat[r_idx][c_idx]
    for r_iter in range(r_idx, n):
        if gcd(mat[r_iter][c_idx], p-1) == 1:
            return r_iter, None, None #can just swap the two rows
        d, x, y = egcd(pivot_val, mat[r_iter][c_idx])
        if d == 1:
            return r_iter, x, y
    return None, None, None


def subtract_rows(row_1, row_2, multiple,p):
    row_1 -= (multiple * row_2) % p
    row_1 %= p
    return row_1

def rref_crt(mat, p, fz):
    reductions = {}
    for f,exp in fz.items():
        if exp > 0:
            print("Factor:", f)
            reductions[f] = rref_prime(np.copy(mat), f).table()
            r = reductions[f]

            for i in range(len(r[0])-1):

                row = r[i]
                print(row[i], end='')
            print()

    vals_for_idx = defaultdict(lambda: ([],[]))
    for f,reduced in reductions.items():
        for i in range(len(reduced[0])-1):
            row = reduced[i]
            mods, vals = vals_for_idx[int(i)]
            mods.append(f)
            vals.append(int(row[-1]))
            vals_for_idx[int(i)] = mods,vals

    sols = {}
    for k,v in vals_for_idx.items():
        n, a = v
        print(k, n, a)
        sols[k] = chinese_remainder(n, a)
    return sols

def rref_prime(mat, p):
    res, rank = flint.nmod_mat(mat.tolist(), p).rref()
    print(p, ":", rank, "=", len(mat[0]), "=", len(mat))
    return res
    return rref(mat,p+1)
    r_index, c_index = 0, 0
    while r_index < mat.shape[0] and c_index < mat.shape[1]:
        val = int(mat[r_index][c_index])

        if val != 1:
            vinv = prime_mod_inv(val, p)
            mat[r_index] = (mat[r_index] * vinv) % p

        for r in range(0, mat.shape[0]):
            if r != r_index:
                below = mat[r][c_index]

                if below != 0:
                    mat[r] = (mat[r] - below * mat[r_index]) % p

        r_index += 1
        c_index += 1

    return mat

def rref(mat, p):
    n = len(mat)
    r_index, c_index = 0, 0
    bottom = n-1
    while r_index < mat.shape[0] and c_index < mat.shape[1]:
        if mat[r_index][c_index] == 0 or gcd(mat[r_index][c_index], p - 1) != 1:
            r_idx_to_use, x, y = find_row_to_elim_nonzero(mat, r_index, c_index, n, p)
            if r_idx_to_use is None:
                #return  # Failure
                r_index+=1
                c_index+=1
                continue

            if x is not None:
                mat[r_index] = (x * mat[r_index] + y * mat[r_idx_to_use]) % (p - 1)
                r_index += 1
                c_index += 1
                continue
            else:  # means we're at the first return case, i.e. just need to swap rows
                mat[[r_idx_to_use, r_index]] = mat[[r_index, r_idx_to_use]]
        div = mod_inv(mat[r_index][c_index], p - 1)
        for k in range(0, mat.shape[1]):
            if k != c_index:
                mat[r_index][k] = mulmod(mat[r_index][k], div, p-1)#(mat[r_index][k] * div) % (p - 1)
            mat[r_index][c_index] = 1
        for k in range(0, mat.shape[0]):
            if k != r_index:
                mat[k] = subtract_rows(mat[k], mat[r_index], mat[k][c_index], p - 1)
        r_index += 1
        c_index += 1
    return mat % (p-1)


def zero_row_absent(matrix):
    if matrix is None:
        return False
    for i in matrix[-1]:
        if i != 0:
            return True
    return False


def linearly_independent(prime_factors, factor_products, p):
    factor_products.append(prime_factors)
    matrix = np.array(factor_products)
    factor_products.pop()
    matrix = rref(matrix, p)
    return zero_row_absent(matrix)


def solve_calculus(beta, p, alpha, prime_discrete_logs, factor_base):
    ret = None
    for exp in range(1, p-1):
        calc = (beta*pow(alpha, exp, p)) % p
        prime_factors = prime_factorization(calc, factor_base)
        if prime_factors:
            ret = -exp % (p-1)
            for idx, ele in enumerate(prime_factors):
                if ele > 0:
                    #ret = (ret + (ele * prime_discrete_logs[idx][-1])) % (p-1)
                    ret = (ret + (ele * prime_discrete_logs[idx])) % (p-1)

            return ret

    return ret

def linearly_independent_mod_2(mat):
    pass

def find_maximally_ind_set(mat, p, fz):
    all_i = set(range(len(mat[0])))
    ind = set(range(len(mat[0])))

    for f,exp in fz.items():
        if exp > 0:
            ind_f = set()
            reduced = rref_prime(np.copy(mat).transpose(), f).table()
            r, c = 0, 0

            if  f == 2:
                print(np.array(reduced))

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
    factor_base = precompute_prime_base(find_bound(p))
    unused = set(factor_base)
    factor_products = []
    checked = set()
    order_fz = factor_by_division(p-1)
    print (order_fz)
    if not order_fz:
        print("Couldn't factor p-1!")
        exit(0)
    real_discrete_logs = dict()
    while len(checked) < p-1:
        while len(unused) > 0 or len(factor_products) < len(factor_base)+ EXCESS_EQUATIONS:
            exp = random.randint(1, p - 1)
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
            factor_products = [e for i,e in enumerate(factor_products) if i in independent]
            prime_discrete_logs = rref_crt(np.array(factor_products), p, order_fz)
            check_crt_results(prime_discrete_logs, factor_base, alpha, p)
            ret = solve_calculus(beta, p, alpha, prime_discrete_logs, factor_base)

            if ret:
                #print(np.array(factor_products))
                return ret
            else:
                print("We were unsuccessful. Try again!")
        else:
            factor_products = [e for i,e in enumerate(factor_products) if i in independent]


def check_crt_results(discrete_log_map, factor_base, g, p):
    successful = dict()

    for k,v in discrete_log_map.items():
        power = pow(g, v, p)
        print(power, "==", factor_base[k])

        if power == factor_base[k]:
            successful[factor_base[k]] = v

    return successful

def find_bound(p):
    return int(3.33*L_p(p, 0.5, 0.476))
    #return int(math.exp((2**(-0.5)*(math.log(2*n)*math.log(math.log(2*n)))**0.5)))

def L_p(p,s,c):
    return math.exp(c*(math.log(p)**s)*(math.log(math.log(p)))**s)


def main():
    if len(sys.argv) < 4:
        print("Usage: python index.py beta p alpha, where beta = alpha^x (mod p)")
        return 1

    beta = int(sys.argv[1])
    p = int(sys.argv[2])
    alpha = int(sys.argv[3])
    x = index_calculus(beta, p, alpha)
    # print(alpha)
    print("The value of x is",  x)
    print("Verification:", beta, "==", pow(alpha, x, p))
    return 0


if __name__ == "__main__":
   main()
