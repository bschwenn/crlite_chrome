import json
from index import index_calculus

def main():
    pub_dict = dict()

    with open("key.pub") as df:
        pub_dict = json.load(df)

    g_to_a = int(pub_dict['g_to_a'])
    p = int(pub_dict['prime'])
    g = int(pub_dict['generator'])
    x = int(index_calculus(g_to_a, p, g))
    print("p", p, "g", g, "power", g_to_a)
    print("The value of x is",  x)
    print("Verification:", g_to_a, "==", pow(g, x, p))

if __name__ == "__main__":
   main()
