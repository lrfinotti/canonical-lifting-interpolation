import sys
import tracemalloc

load('witt/witt.sage')
load('witt/canonical_lifting.sage')
load('witt/canonical_lifting_interpolate.sage')

def print_usage_and_exit():
    print('Usage: profile.sage <profile_type> <lift_alg> <p> <prec>')
    print('  profile_type = speed | memory')
    print('      lift_alg = classic | interp')
    print('             p = prime greater than 3')
    print('          prec = number of components to calculate')
    exit(1)

def parse_args_or_exit():
    if len(sys.argv) < 5: print_usage_and_exit()
    
    profile_type = sys.argv[1]
    if profile_type not in ['speed', 'memory']: print_usage_and_exit()
    
    lift_alg = sys.argv[2]
    if lift_alg not in ['classic', 'interp']: print_usage_and_exit()
    
    p = sys.argv[3]
    try:
        p = int(p)
    except:
        print_usage_and_exit()
    if p < 5 or p not in Primes(): print_usage_and_exit()
    
    prec = sys.argv[4]
    try:
        prec = int(prec)
    except:
        print_usage_and_exit()
    
    return profile_type, lift_alg, p, prec

def magma_is_installed():
    return False

def main():
    profile_type, lift_alg, p, prec = parse_args_or_exit()
    
    if lift_alg == 'classic':
        command = f"F.<a,b> = Frac(GF({p})['a,b']); _ = canonical_lifting_odd_char(a, b, {p}, {prec})"
    else:
        command = f"_ = canonical_lifting_odd_char_interpolate_formulas({p}, {prec})"
    
    if profile_type == 'speed':
        print(f'Timing "{command}"')
        print(timeit(command, number=5))
    else:
        print(f'Analyzing Memory Usage of "{command}"')
        tracemalloc.start()
        
        if lift_alg == 'classic':
            F.<a,b> = Frac(GF(p)['a,b'])
            _ = canonical_lifting_odd_char(a, b, p, prec)
        else:
            _ = canonical_lifting_odd_char_interpolate_formulas(p, prec)
        
        current, peak = tracemalloc.get_traced_memory()
        print(f'Usage (in bytes): {current=}, {peak=}')
        tracemalloc.stop()

if __name__ == "__main__":
    main()
