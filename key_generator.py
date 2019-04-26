import argparse
import json
import math
import random
import sys
from number_theory import *

SEARCH_INTERVAL = 1000
SMALL_PRIMES_BOUND = 10000

def generate_and_write_keys(digits):
    prime = generate_prime(digits)
    generator = select_generator(prime)
    priv_key_a = generate_private_key(prime, generator)
    priv_key_b = generate_private_key(prime, generator)
    g_to_a = pow(generator, priv_key_a, prime)
    g_to_b = pow(generator, priv_key_b, prime)
    shared_key = pow(g_to_a, priv_key_b, prime)
    assert(shared_key == pow(g_to_b, priv_key_a, prime))

    pub_dict = { 'prime':prime, 'generator':generator }

    with open('key.pub', 'w') as outfile:
        json.dump(pub_dict, outfile)

    priv_a_f = open("a_key.priv", "w")
    priv_a_f.write(str(priv_key_a) + "\n")
    priv_a_f.close()

    priv_b_f = open("b_key.priv", "w")
    priv_b_f.write(str(priv_key_b) + "\n")
    priv_b_f.close()

def generate_prime(num_digits):
    # End of the range should be enough to ensure that there is a prime
    # between x and 10**num_digits (which has num_digits + 1 digits)
    x = random.randrange(10**(num_digits-1), 10**num_digits - num_digits**2)
    return next_prime(x)

def next_prime(num):
    start = num
    small_primes = sieve_era(SMALL_PRIMES_BOUND)

    while start < num**2: # arbitrary failure bound
        sieve = sieve_in_range(small_primes, start, start + SEARCH_INTERVAL)
        for candidate in range(start, start + SEARCH_INTERVAL):
            index = candidate - start

            if sieve[index] and is_prime(candidate):
                return candidate

        start += SEARCH_INTERVAL

    return -1


# Per https://crypto.stackexchange.com/questions/820/how-does-one-calculate
# -a-primitive-root-for-diffie-hellman, it seems that it is sufficient in
# practice to just take a random integer modulo p. TODO: confirm this.
def select_generator(p):
    return random.randrange(2, p)

def generate_private_key(prime, generator):
    a = random.randrange(2, prime-1)
    return pow(generator, a, prime)


# Command-line interface stuff
def main():
    parser = argparse.ArgumentParser(description='Diffie-Hellman key generator. Produces files key.pub, a_key.priv, b_key.priv.')
    parser.add_argument('--digits', type=int, help='Number of digits in the public key, when generating a public key. Default is 30.', default=30)

    args = parser.parse_args()
    generate_and_write_keys(args.digits)

    return 0

if __name__ == "__main__":
   main()
