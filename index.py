"""
Usage: python index.py beta p alpha
where beta is some number such that beta \equiv alpha^x (mod p), where alpha is a primitive root mod p
"""
import sys
import math
from sympy.ntheory import factorint
import random
import flint
import numpy as np
from collections import defaultdict
from number_theory import *

EXCESS_EQUATIONS = 0

def precompute_prime_base(B):
    return sieve_era(B)


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
            r = reductions[f]

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
        sols[k] = chinese_remainder(n, a)

    return sols

def rref_prime(mat, p):
    if p == 2:
        return rref_mod_2(mat)
    res, _ = flint.nmod_mat(mat.tolist(), p).rref()
    return res.table()



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
            ret = (-exp) % (p-1)

            for idx, ele in enumerate(prime_factors):
                if ele > 0:
                    ret = (ret + (ele * prime_discrete_logs[idx])) % (p-1)

            return ret

    return ret

def rref_mod_2(mat):
    res = np.copy(mat) % 2
    r, c = 0, 0

    while r < len(res) and c < len(res[0]):
        if res[r][c] != 1:
            swap_idx = find_nonzero(res, c, r)

            if swap_idx:
                res[[r, swap_idx]] = res[[swap_idx, r]]

        if res[r][c] == 1:
            clear_col(res, r, c)

        r += 1
        c += 1

    return res


def find_nonzero(mat, col_idx, into_row_idx):
    for r in range(into_row_idx+1, len(mat)):
        if mat[r][col_idx] != 0:
            return r

def clear_col(mat, row_idx, col_idx):
    for r in range(0, len(mat)):
        if r != row_idx and mat[r][col_idx] == 1:
            mat[r] = mat[r] ^ mat[row_idx]

def find_maximally_ind_set(mat, p, fz):
    all_i = set(range(len(mat[0])))
    ind = set(range(len(mat[0])))

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
    factor_base = precompute_prime_base(B)
    unused = set(factor_base)
    factor_products = []
    checked = set()
    order_fz = factorint(p-1)

    while len(checked) < p-1:
        while len(unused) > 0 or len(factor_products) < len(factor_base)+ EXCESS_EQUATIONS:
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
            factor_products = [e for i,e in enumerate(factor_products) if i in independent]
            prime_discrete_logs = rref_crt(np.array(factor_products), p, order_fz)
            #good_indices = check_crt_results(prime_discrete_logs, factor_base, alpha, p)

            #if len(good_indices) < math.sqrt(B):
            #    print("We were unsuccessful. Make sure the base is actually a primitive root, and try again!")
            #    exit(0)

            #prime_discrete_logs = [p for i,p in prime_discrete_logs.items() if i in good_indices]
            prime_discrete_logs = check_crt_results(prime_discrete_logs, factor_base, alpha, p)

            #factor_base_to_use = [p for i,p in enumerate(factor_base) if i in good_indices]
            ret = solve_calculus(beta, p, alpha, prime_discrete_logs, factor_base)#_to_use)

            if ret:
                return ret
            else:
                print("We were unsuccessful. Try again!")
        else:
            factor_products = [e for i,e in enumerate(factor_products) if i in independent]


def check_crt_results(discrete_log_map, factor_base, g, p):
    successful = []
    new_map = dict()
    for k,v in discrete_log_map.items():
        power = pow(g, v, p)
        neg_prime = (-factor_base[k]) % p

        # If g is truely a primitive root, then if we end up getting x
        # with g^x = -s (mod p) for some prime s in the factor base,
        # then we know that s = g^(x - (p-1)/2) (mod p). This is true
        # because p-1 = g^((p-1)/2) (mod p) when g is a primitive root,
        # since -1 = p - 1 (mod p) and g^(p-1) = 1 (mod p).
        if power == neg_prime:
            new_map[k] = (v - int((p-1)/2)) % (p-1)
        else:
            new_map[k] = v

    return new_map


def find_bound(p):
    return int(3.33*L_p(p, 0.5, 0.476))


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
    print("The value of x is",  x)
    print("Verification:", beta, "==", pow(alpha, x, p))
    return 0


if __name__ == "__main__":
   main()
