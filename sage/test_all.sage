# load necessary files
load('witt/witt.sage')
load('witt/canonical_lifting.sage')
load('witt/canonical_lifting_interpolate.sage')

def test_all(p, prec, n_proc):
    """
    Test all procedures to compute formulas for the canonical lifting and
    elliptic Teichmuller lift.

    p: characteristic (prime)
    prec: precision (number of coordinates)
    n_proc: number of processes to be used in parallelized version
    """


    # symbolic field for formulas
    F.<a, b> = Frac(GF(p)['a,b'])

    # "usual" way
    print("Computing the \"usual way\": ")
    time res1 = canonical_lifting_odd_char(a, b, p, prec)
    print("Done!\n\n")

    # interpolation
    print("Computing with interpolation: ")
    time res2 = canonical_lifting_odd_char_interpolate_formulas(p, prec)
    print("Done!\n\n")

    # interpolation parallel
    print(f"Computing with interpolation in parallel ({n_proc} cores): ")
    time res3 = canonical_lifting_odd_char_interpolate_formulas_parallel(p, prec, n_proc=n_proc)
    print("Done!\n\n")

    P1 = res1[1][1][1].parent()
    P2 = res2[1][1][1].parent()
    P3 = res3[1][1][1].parent()


    phi = P1.hom([P2.0, P2.1, 0, 0, 0, 0])
    psi = P1.hom([P3.0, P3.1, 0, 0, 0, 0])

    print("Check if same results:")
    print((res1[0] == res2[0]) and (phi(res1[1][0][1]) == res2[1][0][1]) and (phi(res1[1][1][1]) == res2[1][1][1]) and
          (res1[0] == res3[0]) and (psi(res1[1][0][1]) == res3[1][0][1]) and (psi(res1[1][1][1]) == res3[1][1][1]))
