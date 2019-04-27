"""
Usage: python index.py beta p alpha
where beta is some number such that beta \equiv alpha^x (mod p), where alpha is a primitive root mod p
"""
import sys
import math
import random
import flint
import numpy as np
from collections import defaultdict
from number_theory import *

EXCESS_EQUATIONS = 5

def precompute_prime_base(B):
    """
    :param B: bound for our index calculus
    :return: list of primes below B

    Basically just the sieve of Erastothenes
    """
    return sieve_era(B)


def prime_factorization(n, primes):
    ret = dict()
    for p in primes:
        ret[p] = 0
        while n % p == 0:
            n //= p
            ret[p] += 1
    return [ret[i] for i in ret] if n == 1 else None


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
            reductions[f] = rref_prime(np.copy(mat), f)

    vals_for_idx = defaultdict(lambda: ([],[]))
    for f,reduced in reductions.items():
        for i, row in enumerate(reduced.table()):
            mods, vals = vals_for_idx[int(i)]
            mods.append(f)
            vals.append(int(row[-1]))
            vals_for_idx[int(i)] = mods,vals

    sols = {}
    for k,v in vals_for_idx.items():
        n, a = v
        sols[k] = chinese_remainder(n, a)

    return sols

def rref_prime(mat, p):
    res, rank = flint.nmod_mat(mat.tolist(), p).rref()
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

    #print(mat)
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

def find_maximally_ind_set(mat, p):
    trans = mat.transpose()
    reduced = rref(trans, p)
    ind = set()
    for r in range(len(reduced)):
        for c in range(len(reduced[r])):
            if reduced[r][c] == 1:
                ind.add(c)
            if reduced[r][c] != 0:
                break

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
    factors_pm1 = prime_factorization(p-1, factor_base)
    factors_pm1 = factor_by_division(p-1)
    if not factors_pm1:
        exit(0)
    pm1_map = factors_pm1
    #for i in range(len(factors_pm1)):
    #    pm1_map[factor_base[i]] = factors_pm1[i]
    while len(checked) < p-1:
        while len(unused) > 0 or len(factor_products) < len(factor_base):# + EXCESS_EQUATIONS:
            #if len(checked) == p - 1:
            #    print("The bound B wasn't big enough!")
            #    exit(0)
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
                if True:#linearly_independent(prime_factors, factor_products, p):
                    for index, i in enumerate(prime_factors[:-1]):
                        if factor_base[index] in unused:
                            unused.remove(factor_base[index])
                    factor_products.append(prime_factors)
        # print(factor_products)
        mat = np.array(factor_products)
        ind = find_maximally_ind_set(np.copy(mat), p)
        if len(ind) == len(mat):

            prime_discrete_logs = rref_crt(np.array(factor_products), p, pm1_map)
            ret = solve_calculus(beta, p, alpha, prime_discrete_logs, factor_base)
            if ret:
                return ret
            else:
                print("We were unsuccessful. Try again!")
        else:
            factor_products = [e for i,e in enumerate(factor_products) if i in ind]


def find_bound(n):
    return int(math.exp((2**(-0.5)*(math.log(2*n)*math.log(math.log(2*n)))**0.5)))


def main():
    if len(sys.argv) < 4:
        print("Usage: python index.py beta p alpha, where beta = alpha^x (mod p)")
        return 1

    #logging.basicConfig(filename='index.log', filemode='w', level=logging.DEBUG, format='%(message)s')

    beta = int(sys.argv[1])
    p = int(sys.argv[2])
    alpha = int(sys.argv[3])
    # print(alpha)
    print("The value of x is",  index_calculus(beta, p, alpha))
    return 0


if __name__ == "__main__":
   main()
