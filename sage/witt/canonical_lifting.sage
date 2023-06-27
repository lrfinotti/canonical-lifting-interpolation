from datetime import datetime

def x_n_order_bound(n, p):
    return -(n+2)*p^n + n*p^(n-1)

def y_n_order_bound(n, p):
    return -(n+3)*p^n + n*p^(n-1)

def F_n_degree_bound(n, p):
    return -x_n_order_bound(n, p)//2

def H_n_degree_bound(n, p):
    return (-y_n_order_bound(n, p) - 3)//2

def poly_order(poly):
    return min(mon_order(mon) for mon in poly.monomials())

def mon_order(mon):
    x0_degree, y0_degree = mon.degrees()[:2]
    return -2*x0_degree + -3*y0_degree

def formal_integrate(f, var, max_degree, include_pth_power_coefficients=True):
    R = f.parent()
    p = R.characteristic()
    F = R(0)
    g = f.polynomial(var)
    for mon in g.monomials():
        degree = mon.degree()
        coeff = g.monomial_coefficient(mon)
        F += coeff / (degree+1) * var^(degree+1)
    if include_pth_power_coefficients:
        num_unknowns = max_degree//p + 1
        new_var_names = ['c' + str(i) for i in range(num_unknowns)]
        S = PolynomialRing(R, new_var_names)
        unknowns = S.gens()
        for n in range(num_unknowns):
            F += unknowns[n] * var^(n*p)
        flatten = S.flattening_morphism()
        return flatten(F)
    else:
        return F

def teichmuller_x_coord_derivative(f, n, prev_x_coords, hasse_invariant=None):
    if len(prev_x_coords) < n:
        raise ValueError(f'We need to know the previous {n} x-coordinates ' +
                          'of the teichmuller to find the next one.')

    if n == 0:
        return f.parent().one()

    P = f.parent()
    p = P.characteristic()
    if p == 2:
        return P.zero()

    x0 = P.gens()[0]
    if hasse_invariant is None:
        hasse_invariant = (f^((p-1)//2)).monomial_coefficient(x0^(p-1))

    F = prev_x_coords # rename to match paper
    final_term = sum(F[i]^(p^(n-i) - 1) * F[i].derivative(x0) for i in range(1, n))
    return hasse_invariant^(-(p^n-1)//(p-1)) * f^((p^n-1)//2) - x0^(p^n-1) - final_term

# def magma_is_installed():
#     #return False
#     try:
#         _ = magma.version()
#     except RuntimeError:
#         return False
#     return True

# def magma_solve(A, b, base, compute_kernel_basis):
#     mag_A = magma(A.transpose())
#     mag_b = magma(b)

#     mag_soln, mag_kernel = magma.Solution(mag_A, mag_b, nvals=2)
#     soln = []
#     for i in range(A.ncols()):
#         soln.append(base(str(mag_soln[i+1])))

#     kernel_basis = None
#     if compute_kernel_basis:
#         V2 = VectorSpace(base, A.ncols())
#         mag_basis = magma.Basis(mag_kernel)
#         kernel_basis = []
#         for mag_basis_vec in mag_basis:
#             basis_vec = []
#             for i in range(A.ncols()):
#                 basis_vec.append(base(str(mag_basis_vec[i+1])))
#             kernel_basis.append(V2(basis_vec))

#     return soln, kernel_basis

def make_zero(coefficients, first_variable_index, compute_kernel_basis=True):
    P = coefficients[0].parent()
    B = P.base()
    relevant_gens = P.gens()[first_variable_index:]

    A = []
    b = []
    for coeff in coefficients:
        for mon in relevant_gens:
            A.append(P(coeff).monomial_coefficient(mon))
        b.append(-P(coeff).constant_coefficient())

    M = MatrixSpace(B, len(b), len(A)//len(b))
    V = VectorSpace(B, len(b))

    # We'll use Magma to solve the system, if it's available.
    # Magma is much faster at solving linear systems than Sage.
    # (Magma is much faster than Sage in general.)
    # if magma_is_installed():
    #     soln, kernel_basis = magma_solve(M(A), V(b), B, compute_kernel_basis)
    # else:
    #     soln = M(A).solve_right(V(b))

    soln = M(A).solve_right(V(b))
    kernel_basis = None
    if compute_kernel_basis:
        kernel_basis = M(A).right_kernel(V(b)).basis()

    return soln, kernel_basis

def canonical_lifting_odd_char(a0, b0, p, prec, verbose=False):
    def log(s, should_log):
        if should_log:
            now = datetime.now()
            now_str = now.strftime('%b %d, %H:%M:%S.%f')[:-3]
            print(f'[{now_str}] {s}')

    if 4*a0^3 + 27*b0^2 == 0:
        raise ValueError(f'The coefficients {a} and {b} define a singular curve.')

    # Setup
    log('Starting...', verbose)

    # Determine common parent of a0 and b0
    # We do this to allow a0=1 instead of e.g. a0=GF(5)(1)
    B1 = a0.parent()
    B2 = b0.parent()
    if a0 in B2:   B = B2
    elif b0 in B1: B = B2
    else:
        raise ValueError(f'The coefficients {a} and {b} have no common parent.')

    # Build the the base ring and starting vectors
    R = B['x0,y0,Fn,Hn,an,bn']
    x0, y0 = R.gens()[:2]
    vF = [x0]; vH = [R(1)]; vy = [y0]; va = [a0]; vb = [b0]

    if p == 3:
        f = x0^3 + a0*x0^2 + b0
        h = a0
    else:
        f = x0^3 + a0*x0 + b0
        h = (f^((p-1)//2)).monomial_coefficient(x0^(p-1))

    if h == 0:
        raise ValueError( 'The canonical lifting only exists for ' +
                          'ordinary elliptic curves, but ' +
                         f'E/F_{p} : y0^2 = {f} is supersingular.')

    f_prime = f.derivative(x0)

    # Compute x/y ahead of time, if necessary
    if prec > 2:
        x_vars = [f'x{i}p' for i in range(prec)]
        y_vars = [f'y{i}p' for i in range(prec)]
        P = Frac(PolynomialRing(GF(p), x_vars + y_vars))
        Wp = WittRing(P, prec=prec)
        x_over_y_witt = Wp(x_vars) / Wp(y_vars)
    else:
        x_over_y_witt = None

    log('Computed x_over_y.\n', verbose)

    for n in range(1, prec):
        an, bn, Fn, Hn = canonical_lifting_odd_char_one_step(
            a0, b0, p, n,
            f, h,
            va, vb, vF, vH, vy,
            R, x_over_y_witt,
            log, verbose
        )

        va.append(an)
        vb.append(bn)
        vF.append(Fn)
        vH.append(Hn)
        vy.append(y0*Hn)

        log('Appended.\n', verbose)

    return (va, vb), (vF, vy)

def canonical_lifting_odd_char_one_step(
    a0, b0, p, n,
    f, h,
    va, vb, vF, vH, vy,
    R, x_over_y_witt,
    log, verbose
):
    x0, y0, Fn, Hn, an, bn = R.gens()

    log(f'n = {n}', verbose)
    log('-----', verbose)

    ## Calculate LHS and RHS of Greenberg Transform
    W = WittRing(R, prec=n+1)

    gt_lhs = W(vy[:n] + [y0*Hn])^2
    lhs0 = gt_lhs.vec[n]
    log('Computed lhs.', verbose)

    if p == 3:
        gt_rhs = W(vF[:n] + [Fn])^3 \
               + W(va[:n] + [an])*W(vF[:n] + [Fn])^2 \
               + W(vb[:n] + [bn])
    else:
        gt_rhs = W(vF[:n] + [Fn])^3 \
               + W(va[:n] + [an])*W(vF[:n] + [Fn]) \
               + W(vb[:n] + [bn])

    rhs0 = gt_rhs.vec[n]
    log('Computed rhs.', verbose)

    # Build F_n including unknowns
    F_n_deriv = teichmuller_x_coord_derivative(f, n, vF, h)
    F_n_with_unknowns = formal_integrate(F_n_deriv, x0, F_n_degree_bound(n, p))
    num_ci_s = F_n_with_unknowns.parent().ngens() - 6
    log('Computed Fn_with_unknowns.', verbose)

    # Enforce τ*(x/y) having a zero at infinity
    if n > 1:
        Mprime = (3*p^(n-1) - 1)//2

        if False:
            # Use Finotti's condition rather than using τ*(x/y)
            x1_sq = vF[1]^2
            known_cis = []
            for i in range(Mprime+1, num_ci_s):
                cprime = x1_sq.monomial_coefficient(x0^(i+p))
                known_cis.append(3*cprime^p / 4)
        else:
            known_cis = canonical_lifting_odd_char_ensure_regularity(
                a0, b0, p, n,
                f, h,
                va, vb, vF, vH, vy,
                R, x_over_y_witt,
                F_n_with_unknowns,
                lhs0, rhs0,
                log, verbose
            )

        still_unknown = F_n_with_unknowns.parent().gens()[:Mprime + 7]
        subs = still_unknown + tuple(known_cis)
        F_n_with_unknowns = F_n_with_unknowns(subs)
        log('  Simplified F_n_with_unknowns.', verbose)

    # Build H_n including unknowns
    num_di_s = H_n_degree_bound(n, p)+1
    new_var_names = ['d' + str(i) for i in range(num_di_s)]
    P = PolynomialRing(F_n_with_unknowns.parent(), new_var_names)
    flatten = P.flattening_morphism()
    H_n_with_unknowns = flatten(sum(P(f'd{i}')*x0^(i) for i in range(num_di_s)))
    log('Computed H_n_with_unknowns.', verbose)

    # Expand all powers of y0 and substitute the expressions
    # for F_n and H_n (with unknowns). The paper guarantees
    # that all powers of y0 will be even.
    lhs1 = 0
    for coeff, mon in lhs0:
        x0_deg = mon.degrees()[0]
        y0_deg = mon.degrees()[1]
        Hn_deg = mon.degrees()[3]
        lhs1 += coeff * x0^x0_deg * _fcppow(f, y0_deg//2, p) * H_n_with_unknowns^Hn_deg

    rhs1 = rhs0(x0, y0, F_n_with_unknowns, H_n_with_unknowns, an, bn)

    # Calculate the entire expression that should be zero.
    eqn = lhs1 - rhs1
    should_be_zero = eqn.polynomial(eqn.parent().gens()[0])
    log('Computed what should be zero.', verbose)

    # Include condition that c_{p^(n-1)} = 0 / a_n = 0 (from paper).
    if p >= 5:
        special_coeff = should_be_zero.parent().base()(f'c{p**(n-1)}')
    else:
        special_coeff = should_be_zero.parent().base()('an')

    coeffs = should_be_zero.coefficients() + [special_coeff]
    soln, kernel_basis = make_zero(coeffs, 3, False)
    log('Solved the system.', verbose)

    subs = [x0, y0, Fn, Hn] + list(soln)

    an = soln[0]
    bn = soln[1]
    Fn = F_n_with_unknowns(subs[:6+num_ci_s])
    Hn = H_n_with_unknowns(subs)

    return an, bn, Fn, Hn

def canonical_lifting_odd_char_ensure_regularity(
    a0, b0, p, n,
    f, h,
    va, vb, vF, vH, vy,
    R, x_over_y_witt,
    F_n_with_unknowns,
    lhs0, rhs0,
    log, verbose
):
    x0, y0, Fn, Hn, an, bn = R.gens()

    log('Enforcing τ*(x/y) being regular at infinity.', verbose)
    x_over_y = x_over_y_witt.vec[n]

    # Compute pullback with current information
    prec = x_over_y.parent().ngens() // 2
    argument = vF[:n] + [Fn] + [0]*(prec-n-1) + vy[:n] + [y0*Hn] + [0]*(prec-n-1)
    pullback_num = x_over_y.numerator()(argument)
    pullback_den = x_over_y.denominator()(argument)
    log('  Computed pullback num and denom.', verbose)

    # Replace y_0^(p^n + 1) Hn by the expression given on RHS of GT
    temp = lhs0.polynomial(Hn)
    subst_expression = rhs0 - temp.constant_coefficient()
    Hn_term      = -x0^(p^n) * y0^((n-1)*p^n + 1) * Hn
    Hn_term_subs = -x0^(p^n) * y0^((n-2)*p^n) * subst_expression / 2
    pullback_num_without_Hn = pullback_num - Hn_term + Hn_term_subs
    log('  Substituted Hn.', verbose)

    # Split up numerator
    scriptF = 0
    scriptG = 0
    scriptH = 0
    for coeff, mon in pullback_num_without_Hn:
        x0_deg, y0_deg = mon.degrees()[:2]
        if mon.degrees()[2] != 0: # Fn degree
            scriptH += coeff * x0^x0_deg * _fcppow(f, y0_deg//2, p)

        if mon.degrees()[2] == 0 and mon.degrees()[3] == 0 and\
           mon.degrees()[4] == 0 and mon.degrees()[5] == 0:
            if y0_deg % 2 == 0 and n%2 == 0:
                scriptF += coeff * x0^x0_deg * _fcppow(f, y0_deg//2, p)
            elif y0_deg % 2 == 1 and n%2 == 1:
                scriptG += coeff * x0^x0_deg * _fcppow(f, y0_deg//2, p)
    log('  Computed scriptF/scriptG and scriptH.', verbose)

    # Throw out all terms with small degree (they don't give information).
    F_n_1 = sum(coeff*mon for coeff, mon in F_n_with_unknowns if mon.degree(x0) >= (3*p^n + 1)/2)

    P = F_n_1.parent()
    if n%2 == 0:
        poly = (scriptH*F_n_1 + scriptF).polynomial(P.gens()[0])
        degree_bound = 3*(n+1)*p^n/2
    else:
        poly = (scriptH*F_n_1 + scriptG).polynomial(P.gens()[0])
        degree_bound = (3*(n+1)*p^n - 3)/2

    # Only hold onto terms that have a large enough degree.
    should_be_zero = 0
    for mon in poly.monomials():
        coeff = poly.monomial_coefficient(mon)
        if mon.degree(x0) >= degree_bound:
            should_be_zero += coeff*mon
    log('  Computed what should be zero.', verbose)

    # Solve the system.
    soln, _ = make_zero(should_be_zero.coefficients(), 3, False)
    log('  Solved system.', verbose)

    # +3 because
    #   We included an and bn in the solve, which adds 2
    #   We solve for c_i for i > M', which add 1
    Mprime = (3*p^(n-1) - 1)//2
    return soln[Mprime + 3:]
