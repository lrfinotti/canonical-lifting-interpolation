def get_ab_powers(w):
    """
    Returns all pairs $(i, j)$ that are solutions to $4i + 6j = w$ with
    $i$ and $j$ non-negative.
    """
    if w%2 != 0: return []
    powers = []
    lb = int(math.ceil(w/6))
    ub = int(math.floor(w/4))
    for k in range(ub, lb-1, -1):
        powers.append((-w//2 + 3*k, w//2 - 2*k))
    return powers

def conjectured_h_bound(p, n):
    return n*p^(n-1) + (n-1)*p^(n-2)

# Only works for p ≥ 5
def canonical_lifting_odd_char_interpolate_formulas(p, prec, verbose=False):
    # Setup final base ring
    F.<a,b> = Frac(GF(p)['a,b'])
    S.<x,y> = F['x,y']
    vF = [x]; vH = [S(1)]; vy = [y]; va = [a]; vb = [b]

    # Calculate Hasse invariant formula
    h = EllipticCurve([a,b]).hasse_invariant()

    # Figure out the ring for j-invariants
    max_weight = 6*p^(prec-1) # weight of b_{prec-1}
    max_unknowns = max_weight//4 - ceil(max_weight/6) + 1
    r = ceil(log(max_unknowns + p^2, p)) # Add p^2 to account for supersingular j-invariants
    B = GF(p^r)

    # Setup ring to be passe to one_step function
    R = B['x0,y0,Fn,Hn,an,bn']
    x0, y0 = R.gens()[:2]

    # Compute x/y ahead of time, if necessary
    if prec > 2:
        x_vars = [f'x{i}p' for i in range(prec)]
        y_vars = [f'y{i}p' for i in range(prec)]
        P = Frac(PolynomialRing(GF(p), x_vars + y_vars))
        Wp = WittRing(P, prec=prec)
        x_over_y_witt = Wp(x_vars) / Wp(y_vars)
    else:
        x_over_y_witt = None

    for n in range(1, prec):
        h_power = conjectured_h_bound(p, n)
        denom_weight = (p-1)*h_power

        # Setup matrices and calculate exponents
        an_numerator_weight = 4*p^n + denom_weight
        an_numerator_powers = get_ab_powers(an_numerator_weight)
        an_unknown_count = len(an_numerator_powers)
        an_matr = MatrixSpace(B, 0, an_unknown_count).an_element()
        an_result = []

        bn_numerator_weight = 6*p^n + denom_weight
        bn_numerator_powers = get_ab_powers(bn_numerator_weight)
        bn_unknown_count = len(bn_numerator_powers)
        bn_matr = MatrixSpace(B, 0, bn_unknown_count).an_element()
        bn_result = []

        ci_numerator_powers = []
        ci_unknown_counts = []
        num_cis = F_n_degree_bound(n, p)//p + 1
        for i in range(num_cis):
            ci_numerator_weight = 2*p^n - 2*i*p + denom_weight
            ci_numerator_powers.append(get_ab_powers(ci_numerator_weight))
            ci_unknown_counts.append(len(ci_numerator_powers[-1]))
        ci_matrs = [MatrixSpace(B, 0, len(l)).an_element() for l in ci_numerator_powers]
        ci_results = [[] for _ in ci_numerator_powers]

        di_numerator_powers = []
        di_unknown_counts = []
        num_dis = H_n_degree_bound(n, p) + 1
        for i in range(num_dis):
            di_numerator_weight = 3*p^n - 3 - 2*i + denom_weight
            di_numerator_powers.append(get_ab_powers(di_numerator_weight))
            di_unknown_counts.append(len(di_numerator_powers[-1]))
        di_matrs = [MatrixSpace(B, 0, len(l)).an_element() for l in di_numerator_powers]
        di_results = [[] for _ in di_numerator_powers]

        # Interpolate
        for j in B:
            # Get a new elliptic curve
            E = EllipticCurve(j=j)
            h0 = E.hasse_invariant()
            if h0 == 0: continue
            a0, b0 = E.a4(), E.a6()
            f = x0^3 + a0*x0 + b0
            h0_to_power = h0^h_power

            # Calculate first n coordinates of canonical lifting
            va0 = [B(ai(a0, b0)) for ai in va]
            vb0 = [B(bi(a0, b0)) for bi in vb]
            vF0 = [sum(B(c(a0,b0))*x0^mon.degree() for c, mon in Fi) for Fi in vF]
            vH0 = [sum(B(c(a0,b0))*x0^mon.degree() for c, mon in Hi) for Hi in vH]
            vy0 = [y0*Hi for Hi in vH0]

            # Get the next step of the canonical lifting
            an, bn, Fn, Hn = canonical_lifting_odd_char_one_step(
                a0, b0, p, n,
                f, h0,
                va0, vb0, vF0, vH0, vy0,
                R, x_over_y_witt,
                lambda x,y: None, False
            )

            done = True

            # Add to system for an
            if an_matr.rank() < an_unknown_count:
                done = False
                row = vector(a0^i*b0^j/h0_to_power for i, j in an_numerator_powers)
                an_matr = an_matr.stack(row)
                an_result.append(B(an))

            # Add to system for bn
            if bn_matr.rank() < bn_unknown_count:
                done = False
                row = vector(a0^i*b0^j/h0_to_power for i, j in bn_numerator_powers)
                bn_matr = bn_matr.stack(row)
                bn_result.append(B(bn))

            # Add to system for c_i
            for k in range(num_cis):
                if ci_matrs[k].rank() < ci_unknown_counts[k]:
                    done = False
                    row = vector(a0^i*b0^j/h0_to_power for i, j in ci_numerator_powers[k])
                    ci_matrs[k] = ci_matrs[k].stack(row)
                    ci_results[k].append(B(Fn.monomial_coefficient(x0^(k*p))))

            # Add to system for d_i
            for k in range(num_dis):
                if di_matrs[k].rank() < di_unknown_counts[k]:
                    done = False
                    row = vector(a0^i*b0^j/h0_to_power for i, j in di_numerator_powers[k])
                    di_matrs[k] = di_matrs[k].stack(row)
                    di_results[k].append(B(Hn.monomial_coefficient(x0^k)))

            if done: break

        # Solve the systems
        soln = an_matr.solve_right(vector(an_result))
        an = sum(F(soln[k])*a^i*b^j for k, (i,j) in enumerate(an_numerator_powers)) / h^h_power
        va.append(an)

        soln = bn_matr.solve_right(vector(bn_result))
        bn = sum(F(soln[k])*a^i*b^j for k, (i,j) in enumerate(bn_numerator_powers)) / h^h_power
        vb.append(bn)

        Fn_deriv = teichmuller_x_coord_derivative(x^3 + a*x + b, n, vF, h)
        Fn = formal_integrate(Fn_deriv, x, F_n_degree_bound(n, p), False)
        for k in range(num_cis):
            soln = ci_matrs[k].solve_right(vector(ci_results[k]))
            ck = sum(F(soln[k])*a^i*b^j for k, (i,j) in enumerate(ci_numerator_powers[k])) / h^h_power
            Fn += ck*x^(p*k)
        vF.append(Fn)

        Hn = 0
        for k in range(num_dis):
            soln = di_matrs[k].solve_right(vector(di_results[k]))
            dk = sum(F(soln[k])*a^i*b^j for k, (i,j) in enumerate(di_numerator_powers[k])) / h^h_power
            Hn += dk*x^k
        vH.append(Hn)
        vy.append(y*Hn)

    return (va, vb), (vF, vy)


# Only works for p ≥ 5
# parallelized version, 8 cores/processors by default
def canonical_lifting_odd_char_interpolate_formulas_parallel(p, prec, verbose=False, n_proc=8):

    @parallel(n_proc)
    def solve_mat(n, M, v):
        return M.solve_right(v)

    @parallel(n_proc)
    def p_canonical_lifting_odd_char_one_step(
        a0, b0, p, n,
        f, h,
        va, vb, vF, vH, vy,
        R, x_over_y_witt,
        log=lambda x, y: None, verbose=False
    ):
        return canonical_lifting_odd_char_one_step(
        a0, b0, p, n,
        f, h,
        va, vb, vF, vH, vy,
        R, x_over_y_witt,
        log, verbose
        )

    # Setup final base ring
    F.<a,b> = Frac(GF(p)['a,b'])
    S.<x,y> = F['x,y']
    vF = [x]; vH = [S(1)]; vy = [y]; va = [a]; vb = [b]

    # Calculate Hasse invariant formula
    h = EllipticCurve([a,b]).hasse_invariant()

    # Figure out the ring for j-invariants
    max_weight = 6*p^(prec-1) # weight of b_{prec-1}
    max_unknowns = max_weight//4 - ceil(max_weight/6) + 1
    r = ceil(log(max_unknowns + p^2, p)) # Add p^2 to account for supersingular j-invariants
    B = GF(p^r)

    # Setup ring to be passe to one_step function
    R = B['x0,y0,Fn,Hn,an,bn']
    x0, y0 = R.gens()[:2]

    # Compute x/y ahead of time, if necessary
    if prec > 2:
        x_vars = [f'x{i}p' for i in range(prec)]
        y_vars = [f'y{i}p' for i in range(prec)]
        P = Frac(PolynomialRing(GF(p), x_vars + y_vars))
        Wp = WittRing(P, prec=prec)
        x_over_y_witt = Wp(x_vars) / Wp(y_vars)
    else:
        x_over_y_witt = None

    for n in range(1, prec):
        h_power = conjectured_h_bound(p, n)
        denom_weight = (p-1)*h_power

        # Setup matrices and calculate exponents
        an_numerator_weight = 4*p^n + denom_weight
        an_numerator_powers = get_ab_powers(an_numerator_weight)
        an_unknown_count = len(an_numerator_powers)
        an_matr = MatrixSpace(B, 0, an_unknown_count).an_element()
        an_result = []

        bn_numerator_weight = 6*p^n + denom_weight
        bn_numerator_powers = get_ab_powers(bn_numerator_weight)
        bn_unknown_count = len(bn_numerator_powers)
        bn_matr = MatrixSpace(B, 0, bn_unknown_count).an_element()
        bn_result = []

        ci_numerator_powers = []
        ci_unknown_counts = []
        num_cis = F_n_degree_bound(n, p)//p + 1
        for i in range(num_cis):
            ci_numerator_weight = 2*p^n - 2*i*p + denom_weight
            ci_numerator_powers.append(get_ab_powers(ci_numerator_weight))
            ci_unknown_counts.append(len(ci_numerator_powers[-1]))
        ci_matrs = [MatrixSpace(B, 0, len(l)).an_element() for l in ci_numerator_powers]
        ci_results = [[] for _ in ci_numerator_powers]

        di_numerator_powers = []
        di_unknown_counts = []
        num_dis = H_n_degree_bound(n, p) + 1
        for i in range(num_dis):
            di_numerator_weight = 3*p^n - 3 - 2*i + denom_weight
            di_numerator_powers.append(get_ab_powers(di_numerator_weight))
            di_unknown_counts.append(len(di_numerator_powers[-1]))
        di_matrs = [MatrixSpace(B, 0, len(l)).an_element() for l in di_numerator_powers]
        di_results = [[] for _ in di_numerator_powers]

        # built input list
        cis_unknown_count = max(ci_unknown_counts)
        dis_unknown_count = max(di_unknown_counts)
        total_unknown = max([an_unknown_count, bn_unknown_count, cis_unknown_count, dis_unknown_count])

        input_list = []
        count = 1
        # Interpolate
        for j in B:
            # Get a new elliptic curve
            E = EllipticCurve(j=j)
            h0 = E.hasse_invariant()
            if h0 == 0: continue
            a0, b0 = E.a4(), E.a6()
            f = x0^3 + a0*x0 + b0
            h0_to_power = h0^h_power

            # Calculate first n coordinates of canonical lifting
            va0 = [B(ai(a0, b0)) for ai in va]
            vb0 = [B(bi(a0, b0)) for bi in vb]
            vF0 = [sum(B(c(a0,b0))*x0^mon.degree() for c, mon in Fi) for Fi in vF]
            vH0 = [sum(B(c(a0,b0))*x0^mon.degree() for c, mon in Hi) for Hi in vH]
            vy0 = [y0*Hi for Hi in vH0]

            input_list.append((a0, b0, p, n,
                f, h0,
                va0, vb0, vF0, vH0, vy0,
                R, x_over_y_witt));

            if count >= total_unknown + 10:
                break

            count += 1

        # compute the canonical liftings in parallel
        print("computing canonical liftings!")
        # with Pool(processes=no_proc) as pool:
        #     results = pool.starmap(canonical_lifting_odd_char_one_step, input_list)
        results = p_canonical_lifting_odd_char_one_step(input_list)
        print("done")

        # print(list(results))
        # return list(results)

        for res in results:
            done = True

            coef = res[0][0]
            a0 = coef[0]
            b0 = coef[1]
            h0 = coef[5]
            h0_to_power = h0^h_power
            an, bn, Fn, Hn = res[-1]

            # Add to system for an
            if an_matr.rank() < an_unknown_count:
                done = False
                row = vector(a0^i*b0^j/h0_to_power for i, j in an_numerator_powers)
                an_matr = an_matr.stack(row)
                an_result.append(B(an))

            # Add to system for bn
            if bn_matr.rank() < bn_unknown_count:
                done = False
                row = vector(a0^i*b0^j/h0_to_power for i, j in bn_numerator_powers)
                bn_matr = bn_matr.stack(row)
                bn_result.append(B(bn))

            # Add to system for c_i
            for k in range(num_cis):
                if ci_matrs[k].rank() < ci_unknown_counts[k]:
                    done = False
                    row = vector(a0^i*b0^j/h0_to_power for i, j in ci_numerator_powers[k])
                    ci_matrs[k] = ci_matrs[k].stack(row)
                    ci_results[k].append(B(Fn.monomial_coefficient(x0^(k*p))))

            # Add to system for d_i
            for k in range(num_dis):
                if di_matrs[k].rank() < di_unknown_counts[k]:
                    done = False
                    row = vector(a0^i*b0^j/h0_to_power for i, j in di_numerator_powers[k])
                    di_matrs[k] = di_matrs[k].stack(row)
                    di_results[k].append(B(Hn.monomial_coefficient(x0^k)))

            if done: break


        # # Solve the systems
        # soln = an_matr.solve_right(vector(an_result))
        # an = sum(F(soln[k])*a^i*b^j for k, (i,j) in enumerate(an_numerator_powers)) / h^h_power
        # va.append(an)

        # soln = bn_matr.solve_right(vector(bn_result))
        # bn = sum(F(soln[k])*a^i*b^j for k, (i,j) in enumerate(bn_numerator_powers)) / h^h_power
        # vb.append(bn)

        # Fn_deriv = teichmuller_x_coord_derivative(x^3 + a*x + b, n, vF, h)
        # Fn = formal_integrate(Fn_deriv, x, F_n_degree_bound(n, p), False)
        # for k in range(num_cis):
        #     soln = ci_matrs[k].solve_right(vector(ci_results[k]))
        #     ck = sum(F(soln[k])*a^i*b^j for k, (i,j) in enumerate(ci_numerator_powers[k])) / h^h_power
        #     Fn += ck*x^(p*k)
        # vF.append(Fn)

        # Hn = 0
        # for k in range(num_dis):
        #     soln = di_matrs[k].solve_right(vector(di_results[k]))
        #     dk = sum(F(soln[k])*a^i*b^j for k, (i,j) in enumerate(di_numerator_powers[k])) / h^h_power
        #     Hn += dk*x^k
        # vH.append(Hn)
        # vy.append(y*Hn)

        input_list = [(0, an_matr, vector(an_result)), (1, bn_matr, vector(bn_result))]
        for k in range(num_cis):
            input_list.append((k + 2, ci_matrs[k], vector(ci_results[k])))
        for k in range(num_dis):
            input_list.append((k + 2 + num_cis, di_matrs[k], vector(di_results[k])))

        # print(len(input_list), input_list)

        print("solving the systems!")
        # with Pool(processes=no_proc) as pool:
        #     results = pool.map(solve_mat, input_list)
        results = list(solve_mat(input_list))

        results.sort()
        print("done!")

        # return results

        soln = results[0][-1]
        an = sum(F(soln[k])*a^i*b^j for k, (i,j) in enumerate(an_numerator_powers)) / h^h_power
        va.append(an)

        soln = results[1][-1]
        bn = sum(F(soln[k])*a^i*b^j for k, (i,j) in enumerate(bn_numerator_powers)) / h^h_power
        vb.append(bn)

        Fn_deriv = teichmuller_x_coord_derivative(x^3 + a*x + b, n, vF, h)
        Fn = formal_integrate(Fn_deriv, x, F_n_degree_bound(n, p), False)
        for kk in range(num_cis):
            soln = results[kk + 2][-1]
            ck = sum(F(soln[k])*a^i*b^j for k, (i,j) in enumerate(ci_numerator_powers[kk])) / h^h_power
            Fn += ck*x^(p*kk)
        vF.append(Fn)

        Hn = 0
        for kk in range(num_dis):
            soln = results[kk + 2 + num_cis][-1]
            dk = sum(F(soln[k])*a^i*b^j for k, (i,j) in enumerate(di_numerator_powers[kk])) / h^h_power
            Hn += dk*x^kk
        vH.append(Hn)
        vy.append(y*Hn)


    return (va, vb), (vF, vy)
