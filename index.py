"""
Usage: python index.py beta p alpha
where beta is some number such that beta \equiv alpha^x (mod p), where alpha is a primitive root mod p
"""
import logging
import sys
import math
import random
import numpy as np


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
    return [ret[i] for i in ret] if n == 1 else None


def egcd(a, b):
    if a == 0:
        return b, 0, 1
    else:
        g, y, x = egcd(b % a, a)
        return g, x - (b // a) * y, y


def find_row_to_elim_nonzero(mat, r_idx, c_idx, n, p):
    pivot_val = mat[r_idx][c_idx]
    for r_iter in range(r_idx, n):
        if gcd(mat[r_iter][c_idx], p-1) == 1:
            return r_iter, None, None #can just swap the two rows
        d, x, y = egcd(pivot_val, mat[r_iter][c_idx])
        if d == 1:
            return r_iter, x, y
    return None, None, None

def solve_relations(matrix, p):
   return

def gcd(a, b):
    return gcd(b, a % b) if a % b else b


def mod_inv(a, n):
    x0, x1 = 1, 0
    while n != 0:
        q, a, n = a // n, n, a % n
        x0, x1 = x1, x0 - q * x1
    return x0


def subtract_rows(row_1, row_2, multiple,p):
    row_1 -= (multiple * row_2) % p
    row_1 %= p
    return row_1


def rref(mat, p):
    n = len(mat)
    r_index, c_index = 0, 0
    while r_index < mat.shape[0] and c_index < mat.shape[1]:
        if mat[r_index][c_index] == 0 or gcd(mat[r_index][c_index], p - 1) != 1:
            r_idx_to_use, x, y = find_row_to_elim_nonzero(mat, r_index, c_index, n, p)
            if r_idx_to_use is None:
                return  # Failure
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
                mat[r_index][k] *= div
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
    matrix = np.array(factor_products, dtype=np.float64)
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
                    ret = (ret + (ele * prime_discrete_logs[idx][-1])) % (p-1)
    return ret



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
        power = pow(alpha, exp, p)
        if power == 1:
            continue
        prime_factors = prime_factorization(power, factor_base)
        if prime_factors:
            prime_factors.append(exp)
            if linearly_independent(prime_factors, factor_products, p):
                for index, i in enumerate(prime_factors[:-1]):
                    if factor_base[index] in unused:
                        unused.remove(factor_base[index])
                factor_products.append(prime_factors)
    # print(factor_products)
    prime_discrete_logs = rref(np.array(factor_products), p)
    ret = solve_calculus(beta, p, alpha, prime_discrete_logs, factor_base)
    if ret:
        return ret
    else:
        print("We were unsuccessful. Try again!")


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
    # print(alpha)
    print("The value of x is",  index_calculus(beta, p, alpha))
    return 0


if __name__ == "__main__":
   main()
