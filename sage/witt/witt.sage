from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.finite_rings.finite_field_base import is_FiniteField
import sage.rings.ring as ring

_CommutativeRings = CommutativeRings()
_Primes = Primes()

def _fast_char_p_power(x, n, p=None):
    r"""
    Raise x^n power in characteristic p

    If x is not an element of a ring of characteristic p, throw an error.
    If x is an element of GF(p^k), this is already fast.
    However, is x is a polynomial, this seems to be slow?
    """
    if n not in ZZ:
        raise ValueError(f'Exponent {n} is not an integer')
    if n == 0 or x == 1:
        return x.parent().one()
    if x.parent().characteristic() not in _Primes:
        raise ValueError(f'{x} is not in a ring of prime characteristic')

    x_is_Polynomial = isinstance(x, Polynomial)
    x_is_MPolynomial = isinstance(x, MPolynomial)

    if not (x_is_Polynomial or x_is_MPolynomial):
        return x**n
    if (x_is_Polynomial and x.is_gen()) or (x_is_MPolynomial and x.is_generator()):
        return x**n
    if n < 0:
        x = x**-1  # This may throw an error.
        n = -n

    P = x.parent()
    if p is None:
        p = P.characteristic()
    base_p_digits = ZZ(n).digits(base=p)

    x_to_the_n = 1

    for p_exp, digit in enumerate(base_p_digits):
        if digit == 0:
            continue
        inner_term = x**digit
        term_dict = {}
        for e_int_or_tuple, c in inner_term.dict().items():
            power = p**p_exp
            new_c = _fast_char_p_power(c, power)
            new_e_tuple = None
            if x_is_Polynomial:  # Then the dict keys are ints
                new_e_tuple = e_int_or_tuple * power
            elif x_is_MPolynomial:  # Then the dict keys are ETuples
                new_e_tuple = e_int_or_tuple.emul(power)
            term_dict[new_e_tuple] = new_c
        term = P(term_dict)
        x_to_the_n *= term

    return x_to_the_n


_fcppow = _fast_char_p_power

class WittVector_base(CommutativeRingElement):
    def __init__(self, parent, vec=None):
        self.prec = parent.precision()
        B = parent.base()
        if vec is not None:
            if len(vec) != self.prec:
                raise ValueError(f'{vec} is not the correct length. Expected length to be {self.prec}.')
            self.vec = tuple(B(x) for x in vec)
        else:
            self.vec = (B(0) for i in range(self.prec))
        CommutativeRingElement.__init__(self, parent)

    def __hash__(self):
        return hash(self.vec)

    def _richcmp_(self, other, op):
        from sage.structure.richcmp import op_EQ, op_NE
        if op == op_EQ:
            return self.vec == other.vec
        elif op == op_NE:
            return self.vec != other.vec
        else:
            return NotImplemented

    def _repr_(self):
        return '(' + ', '.join(map(str, self.vec)) + ')'

    def _add_(self, other):
        P = self.parent()
        C = self.__class__

        # As a slight optimization, we'll check for zero ahead of time.
        # This has the benefit of allowing us to create polynomials,
        # even if ``P._algorithm`` is 'none'.
        if other == P.zero():
            return self
        elif self == P.zero():
            return other

        alg = P._algorithm
        if alg == 'standard':
            s = P.sum_polynomials
            # note here this is tuple addition, i.e. concatenation
            sum_vec = tuple(s[i](*(self.vec + other.vec)) for i in range(self.prec))
            return C(P, vec=sum_vec)
        else:
            return NotImplemented

    def _mul_(self, other):
        P = self.parent()
        C = self.__class__

        # As a slight optimization, we'll check for zero or one ahead of time.
        # This has the benefit of allowing us to create polynomials,
        # even if ``P._algorithm`` is 'none'.
        if self == P.zero() or other == P.zero():
            return P.zero()
        elif other == P.one():
            return self
        elif self == P.one():
            return other

        alg = P._algorithm
        if alg == 'standard':
            p = P.prod_polynomials
            # note here this is tuple addition, i.e. concatenation
            prod_vec = tuple(p[i](*(self.vec + other.vec)) for i in range(self.prec))
            return C(P, vec=prod_vec)
        else:
            return NotImplemented

    def _neg_(self):
        P = self.parent()
        C = self.__class__
        # If p == 2, -1 == (-1, -1, -1, ...)
        # Otherwise, -1 == (-1, 0, 0, ...)
        if P.prime == 2:
            all_ones = P(tuple(-1 for _ in range(self.prec)))
            return all_ones*self
        neg_vec = tuple(-self.vec[i] for i in range(self.prec))
        return C(P, vec=neg_vec)

    def _div_(self, other):
        P = self.parent()
        # As a slight optimization, we'll check for one ahead of time.
        # This has the benefit of allowing us to create polynomials,
        # even if ``P._algorithm`` is 'none'.
        if other == P.one():
            return self
        elif self == P.one():
            return other._invert_()

        return self * other._invert_()

    def _invert_(self):
        if not self.vec[0].is_unit():
            raise ZeroDivisionError(f"Inverse of {self} does not exist.")
        P = self.parent()
        C = self.__class__

        if self == P.one():
            return self
        if self.prec == 1:
            return P((self.vec[0]**-1, ))

        # Strategy: Multiply ``self`` by ``(Y_0, Y_1, ...)``, set equal
        # to (1, 0, 0, ...), and solve.
        var_names = [f'Y{i}' for i in range(1, self.prec)]
        poly_ring = PolynomialRing(P.base(), var_names)
        inv_vec = list((self.vec[0]**-1, ) + poly_ring.gens()) # We'll fill this in one-by-one

        W = WittRing(poly_ring, p=P.prime, prec=P.prec)
        prod_vec = (W(self.vec)*W(inv_vec)).vec
        for i in range(1, self.prec):
            poly = prod_vec[i](inv_vec[1:])
            Y_i = poly.parent().gens()[i-1]
            inv_vec[i] = -poly.constant_coefficient() / poly.monomial_coefficient(Y_i)

        return C(P, vec=inv_vec)


class WittVector_p_typical(WittVector_base):
    def _add_(self, other):
        P = self.parent()
        C = self.__class__

        # As a slight optimization, we'll check for zero ahead of time.
        # This has the benefit of allowing us to create polynomials,
        # even if ``P._algorithm`` is 'none'.
        if other == P.zero():
            return self
        elif self == P.zero():
            return other

        alg = P._algorithm
        if alg == 'standard':
            s = P.sum_polynomials
            # note here this is tuple addition, i.e. concatenation
            sum_vec = tuple(s[i](*(self.vec + other.vec)) for i in range(self.prec))
            return C(P, vec=sum_vec)
        elif alg == 'finotti':
            x = self.vec
            y = other.vec
            prec = P.precision()
            p = P.prime
            bin_table = P.binomial_table

            G = []
            for n in range(0, prec):
                G_n = [x[n], y[n]]
                for i in range(0, n):
                    G_n.append(P.eta_bar(G[i], n-i))
                G.append(G_n)
            sum_vec = tuple(sum(G[i]) for i in range(prec))
            return C(P, vec=sum_vec)
        elif alg == 'Zq_isomorphism':
            x = P._vector_to_series(self.vec)
            y = P._vector_to_series(other.vec)
            sum_vec = P._series_to_vector(x+y)
            return C(P, vec=sum_vec)
        else:
            return NotImplemented

    def _mul_(self, other):
        P = self.parent()
        C = self.__class__

        # As a slight optimization, we'll check for zero or one ahead of time.
        # This has the benefit of allowing us to create polynomials,
        # even if ``P._algorithm`` is 'none'.
        if self == P.zero() or other == P.zero():
            return P.zero()
        elif other == P.one():
            return self
        elif self == P.one():
            return other

        alg = P._algorithm
        if alg == 'standard':
            p = P.prod_polynomials
            # note here this is tuple addition, i.e. concatenation
            prod_vec = tuple(p[i](*(self.vec + other.vec)) for i in range(self.prec))
            return C(P, vec=prod_vec)
        elif alg == 'finotti':
            x = self.vec
            y = other.vec
            prec = P.precision()
            p = P.prime
            bin_table = P.binomial_table

            G = [[x[0] * y[0]]]
            for n in range(1, prec):
                G_n = [_fcppow(x[0], p**n) * y[n], _fcppow(y[0], p**n) * x[n]]
                G_n.extend(_fcppow(x[i], p**(n-i)) * _fcppow(y[n-i], p**i) for i in range(1, n))
                for i in range(0, n):
                    G_n.append(P.eta_bar(G[i], n-i))
                G.append(G_n)
            prod_vec = tuple(sum(G[i]) for i in range(prec))
            return C(P, vec=prod_vec)
        elif alg == 'Zq_isomorphism':
            x = P._vector_to_series(self.vec)
            y = P._vector_to_series(other.vec)
            sum_vec = P._series_to_vector(x*y)
            return C(P, vec=sum_vec)
        else:
            return NotImplemented


class WittVector_non_p_typical(WittVector_base):

    def _add_(self, other):
        P = self.parent()
        C = self.__class__

        # As a slight optimization, we'll check for zero ahead of time.
        # This has the benefit of allowing us to create polynomials,
        # even if ``P._algorithm`` is 'none'.
        if other == P.zero():
            return self
        elif self == P.zero():
            return other

        alg = P._algorithm
        if alg == 'standard':
            s = P.sum_polynomials
            # note here this is tuple addition, i.e. concatenation
            sum_vec = tuple(s[i](*(self.vec + other.vec)) for i in range(self.prec))
            return C(P, vec=sum_vec)
        elif alg == 'standard_otf':
            p = P.prime # we know p is a unit in this case!
            x = self.vec
            y = other.vec

            sum_vec = [x[0] + y[0]]
            for n in range(1, self.prec):
                next_sum = x[n] + y[n] + \
                    sum((x[i]**(p**(n-i)) + y[i]**(p**(n-i)) - sum_vec[i]**(p**(n-i))) / p**(n-i)
                        for i in range(0, n))
                sum_vec.append(next_sum)

            return C(P, vec=sum_vec)
        elif alg == 'IntegerMod_isomorphism':
            a = P._vector_to_coefficients(self)
            b = P._vector_to_coefficients(other)
            sum_coeffs = tuple(a[i] + b[i] for i in range(self.prec))
            return C(P, vec=P._coefficients_to_vector(sum_coeffs))
        else:
            return NotImplemented

    def _mul_(self, other):
        P = self.parent()
        C = self.__class__

        # As a slight optimization, we'll check for zero or one ahead of time.
        # This has the benefit of allowing us to create polynomials,
        # even if ``P._algorithm`` is 'none'.
        if self == P.zero() or other == P.zero():
            return P.zero()
        elif other == P.one():
            return self
        elif self == P.one():
            return other

        alg = P._algorithm
        if alg == 'standard':
            p = P.prod_polynomials
            # note here this is tuple addition, i.e. concatenation
            prod_vec = tuple(p[i](*(self.vec + other.vec)) for i in range(self.prec))
            return C(P, vec=prod_vec)
        elif alg == 'standard_otf':
            p = P.prime # we know p is a unit in this case!
            x = self.vec
            y = other.vec

            prod_vec = [x[0] * y[0]]
            for n in range(1, self.prec):
                next_prod = (
                    sum(p**i * x[i]**(p**(n-i)) for i in range(0, n+1)) *
                    sum(p**i * y[i]**(p**(n-i)) for i in range(0, n+1)) -
                    sum(p**i * prod_vec[i]**(p**(n-i)) for i in range(0, n))
                ) / p**n
                prod_vec.append(next_prod)

            return C(P, vec=prod_vec)
        elif alg == 'IntegerMod_isomorphism':
            p = P.prime

            a = P._vector_to_coefficients(self)
            b = P._vector_to_coefficients(other)

            prod_coeffs = [a[0]*b[0]]
            for i in range(1, self.prec):
                c_i = a[0]*b[i] + a[i]*b[0] + p**i * a[i]*b[i]
                prod_coeffs.append(c_i)

            return C(P, vec=P._coefficients_to_vector(prod_coeffs))
        else:
            return NotImplemented

class WittRing_base(CommutativeRing, UniqueRepresentation):

    Element = WittVector_base

    def __init__(self, base_ring, prec, prime, algorithm='none', category=None):
        self.prec = prec
        self.prime = prime

        self._algorithm = algorithm
        self.sum_polynomials = None
        self.prod_polynomials = None

        if algorithm == 'standard':
            self._generate_sum_and_product_polynomials(base_ring)

        if category is None:
            category = CommutativeRings()
        CommutativeRing.__init__(self, base_ring, category=category)

    def __iter__(self):
        from itertools import product, repeat
        for t in product(self.base(), repeat=self.prec):
            yield self(t)

    def _repr_(self):
        return f"Ring of Witt Vectors of length {self.prec} over {self.base()}"

    def _coerce_map_from_(self, S):
        # Question: do we return True is S == self.base()?
        # We have the teichmuller lift, but I don't think that's
        # a "coercion" map, per se.
        return (S is ZZ)

    def _element_constructor_(self, x):
        if x in ZZ:
            return self.element_class(self, self._int_to_vector(x))
        elif isinstance(x, tuple) or isinstance(x, list):
            return self.element_class(self, x)
        else:
            return NotImplemented

    def _int_to_vector(self, k):
        p = self.prime

        char = self.characteristic()
        if char != 0:
            k = k % char

        vec_k = [k]
        for n in range(1, self.prec):
            total = k - k**(p**n) - sum(p**(n-i) * vec_k[n-i]**(p**i) for i in range(1, n))
            total //= p**n
            vec_k.append(total)

        return vec_k

    def _generate_sum_and_product_polynomials(self, base):
        p = self.prime
        prec = self.prec
        x_var_names = ['X{}'.format(i) for i in range(prec)]
        y_var_names = ['Y{}'.format(i) for i in range(prec)]
        var_names = x_var_names + y_var_names

        # Okay, what's going on here? Sage, by default, relies on
        # Singular for Multivariate Polynomial Rings, but Singular uses
        # only SIXTEEN bits (unsigned) to store its exponents. So if we
        # want exponents larger than 2^16 - 1, we have to use the
        # generic implementation. However, after some experimentation,
        # it seems like the generic implementation is faster?
        #
        # After trying to compute S_4 for p=5, it looks like generic is
        # faster for  very small polys, and MUCH slower for large polys.
        # So we'll default to singular unless we can't use it.
        #
        # Remark: Since when is SIXTEEN bits sufficient for anyone???
        #
        if p**(prec-1) >= 2**16:
            implementation = 'generic'
        else:
            implementation = 'singular'

        # We first generate the "universal" polynomials and then project
        # to the base ring.
        R = PolynomialRing(ZZ, var_names, implementation=implementation)
        x_y_vars = R.gens()
        x_vars = x_y_vars[:prec]
        y_vars = x_y_vars[prec:]

        self.sum_polynomials = [0]*(self.prec)
        for n in range(prec):
            s_n = x_vars[n] + y_vars[n]
            for i in range(n):
                s_n += (x_vars[i]**(p**(n-i)) + y_vars[i]**(p**(n-i)) - self.sum_polynomials[i]**(p**(n-i))) / p**(n-i)
            self.sum_polynomials[n] = R(s_n)

        self.prod_polynomials = [x_vars[0] * y_vars[0]] + [0]*(self.prec)
        for n in range(1, prec):
            x_poly = sum([p**i * x_vars[i]**(p**(n-i)) for i in range(n+1)])
            y_poly = sum([p**i * y_vars[i]**(p**(n-i)) for i in range(n+1)])
            p_poly = sum([p**i * self.prod_polynomials[i]**(p**(n-i)) for i in range(n)])
            p_n = (x_poly*y_poly - p_poly) // p**n
            self.prod_polynomials[n] = p_n

        # We have to use generic here, because Singular doesn't support
        # Polynomial Rings over Polynomial Rings. For example,
        # ``PolynomialRing(GF(5)['x'], ['X', 'Y'], implementation='singular')``
        # will fail.
        S = PolynomialRing(base, x_y_vars, implementation='generic')
        for n in range(prec):
            self.sum_polynomials[n] = S(self.sum_polynomials[n])
            self.prod_polynomials[n] = S(self.prod_polynomials[n])

    def characteristic(self):
        p = self.prime
        if self.base()(p).is_unit():
            # If p is invertible, W_n(R) is isomorphic to R^n.
            return self.base().characteristic()
        else:
            # This is a conjecture. It's known for char(R) == p.
            return p**(self.prec-1) * self.base().characteristic()

    def precision(self):
        return self.prec

    def random_element(self):
        return self.element_class(self, tuple(self.base().random_element() for _ in range(self.prec)))

    def teichmuller_lift(self, x):
        if x not in self.base():
            raise Exception(f'{x} not in {self.base()}')
        else:
            return self.element_class(self, (x,) + tuple(0 for _ in range(self.prec-1)))

    def is_finite(self):
        return self.base().is_finite()

    def cardinality(self):
        return self.base().cardinality()**(self.prec)


class WittRing_p_typical(WittRing_base):

    Element = WittVector_p_typical

    def __init__(self, base_ring, prec, prime, algorithm=None, category=None):
        WittRing_base.__init__(self, base_ring, prec, prime,
                               algorithm=algorithm, category=category)

        if algorithm == 'finotti':
            self.generate_binomial_table()

    def _int_to_vector(self, k):
        p = self.prime
        R = Zp(p, prec=self.prec+1, type='fixed-mod')
        F = GF(p)
        B = self.base()

        series = R(k)
        witt_vector = []
        for _ in range(self.prec):
            # Probably slightly faster to do "series % p," but this way, temp is in F_p
            temp = F(series)
            witt_vector.append(B(temp)) # make sure elements of vector are in base
            series = (series - R.teichmuller(temp)) // p
        return witt_vector

    def generate_binomial_table(self):
        import numpy as np
        p = self.prime
        R = Zp(p, prec=self.prec+1, type='fixed-mod')
        v_p = ZZ.valuation(p)
        table = [[0]]
        for k in range(1, self.prec+1):
            row = np.empty(p**k, dtype=int)
            row[0] = 0
            prev_bin = 1
            for i in range(1, p**k // 2 + 1):
                val = v_p(i)
                # Instead of calling binomial each time, we compute the coefficients
                # recursively. This is MUCH faster.
                next_bin = prev_bin * (p**k - (i-1)) // i
                prev_bin = next_bin
                series = R(-next_bin // p**(k-val))
                for _ in range(val):
                    temp = series % p
                    series = (series - R.teichmuller(temp)) // p
                row[i] = ZZ(series % p)
                row[p**k - i] = row[i] # binomial coefficients are symmetric
            table.append(row)
        self.binomial_table = table

    def eta_bar(self, vec, eta_index):
        vec = tuple(x for x in vec if x != 0) # strip zeroes

        # special cases
        if len(vec) <= 1:
            return 0
        if eta_index == 0:
            return sum(vec)

        # renaming to match notation in paper
        k = eta_index
        p = self.prime
        # if vec = (x,y), we know what to do: Theorem 8.6
        if len(vec) == 2:
            # Here we have to check if we've pre-computed already
            x, y = vec
            scriptN = [[None] for _ in range(k+1)] # each list starts with None, so that indexing matches paper
            # calculate first N_t scriptN's
            for t in range(1, k+1):
                for i in range(1, p**t):
                    scriptN[t].append(self.binomial_table[t][i] * _fcppow(x, i) * _fcppow(y, p**t - i))
            indexN = [p**i - 1 for i in range(k+1)]
            for t in range(2, k+1):
                for l in range(1, t):
                    # append scriptN_{t, N_t+l}
                    next_scriptN = self.eta_bar(scriptN[t-l][1:indexN[t-l]+t-l], l)
                    scriptN[t].append(next_scriptN)
            return sum(scriptN[k][1:])

        # if vec is longer, we split and recurse: Proposition 5.4
        # This is where we need to using multiprocessing.
        else:
            m = len(vec) // 2
            v_1 = vec[:m]
            v_2 = vec[m:]
            s_1 = sum(v_1)
            s_2 = sum(v_2)
            total = 0
            scriptM = [[] for _ in range(k+1)]
            for t in range(1, k+1):
                scriptM[t].append(self.eta_bar(v_1,        t))
                scriptM[t].append(self.eta_bar(v_2,        t))
                scriptM[t].append(self.eta_bar((s_1, s_2), t))
            for t in range(2, k+1):
                for s in range(1, t):
                    result = self.eta_bar(scriptM[t-s], s)
                    scriptM[t].append(result)
            return sum(scriptM[k])


class WittRing_finite_field(WittRing_p_typical):
    def __init__(self, base_ring, prec, prime, category=None):
        WittRing_p_typical.__init__(self, base_ring, prec, prime,
                                    algorithm='Zq_isomorphism',
                                    category=category)

    def _series_to_vector(self, series):
        F = self.base() # known to be finite
        R = Zq(F.cardinality(), prec=self.prec, type='fixed-mod', modulus=F.polynomial(), names=['z'])
        K = R.residue_field()
        p = self.prime

        series = R(series)
        witt_vector = []
        for i in range(self.prec):
            temp = K(series)
            elem = temp.polynomial()(F.gen()) # hack to convert to F (K != F for some reason)
            witt_vector.append(elem**(p**i))
            series = (series - R.teichmuller(temp)) // p
        return witt_vector

    def _vector_to_series(self, vec):
        F = self.base()
        R = Zq(F.cardinality(), prec=self.prec, type='fixed-mod', modulus=F.polynomial(), names=['z'])
        K = R.residue_field()
        p = self.prime

        series = R(0)
        for i in range(0, self.prec):
            temp = vec[i].nth_root(p**i)
            elem = temp.polynomial()(K.gen()) # hack to convert to K (F != K for some reason)
            series += p**i * R.teichmuller(elem)
        return series


class WittRing_non_p_typical(WittRing_base):

    Element = WittVector_non_p_typical

    def __init__(self, base_ring, prec, prime, algorithm=None, category=None):
        WittRing_base.__init__(self, base_ring, prec, prime,
                               algorithm=algorithm, category=category)

    def _repr_(self):
        return f"Ring of {self.prime}-Witt Vectors of length {self.prec} over {self.base()}"


class WittRing_integers_mod_power_of_p(WittRing_non_p_typical):
    def __init__(self, base_ring, prec, prime, category=None):
        self.alpha = ZZ.valuation(prime)(base_ring.characteristic())
        WittRing_base.__init__(self, base_ring, prec, prime,
                               algorithm='IntegerMod_isomorphism',
                               category=category)

    def _coefficients_to_vector(self, c):
        p = self.prime
        n = self.prec
        alpha = self.alpha

        B = self.base()
        R = Zmod(p**(alpha+n-1))

        v = []
        for i in range(n-1):
            # It may appear that we can simplify the computations below
            # by canceling common factors in the multiplication and
            # division. However, since the first computations are done in
            # R and then we switch to ZZ, this cancellation doesn't work.
            v_i = R(c[0]) + p**(i+1)*R(c[i+1])
            v_i -= sum(p**j * R(v[j])**(p**(i-j)) for j in range(i))
            v_i *= p**(n-i-1)
            v_i = ZZ(v_i) // p**(n-1)
            v.append(B(v_i))

        last_v = R(c[0]) - sum(p**j * R(v[j])**(p**(n-j-1)) for j in range(n-1))
        last_v = ZZ(last_v) // p**(n-1)
        v.append(B(last_v))

        return tuple(v)

    def _vector_to_coefficients(self, v):
        p = self.prime
        n = self.prec
        alpha = self.alpha

        R = Zmod(p**(alpha+n-1))
        S = Zmod(p**(alpha-1))

        c0 = sum(p**j * R(v.vec[j])**(p**(n-j-1)) for j in range(n))
        c = [c0]
        for i in range(1, n):
            # It may appear that we can simplify the computations below
            # by canceling common factors in the multiplication and
            # division. However, since the first computations are done in
            # R and then we switch to ZZ, this cancellation doesn't work.
            c_i = sum(p**j * R(v.vec[j])**(p**(i-j-1)) for j in range(i)) - c0
            c_i *= p**(n-i)
            c_i = ZZ(c_i) // p**n
            c.append(S(c_i))

        return tuple(c)


class WittRing_p_invertible(WittRing_non_p_typical):
    def __init__(self, base_ring, prec, prime, category=None):
        WittRing_non_p_typical.__init__(self, base_ring, prec, prime,
                                        algorithm='standard_otf',
                                        category=category)


def WittRing(base_ring, prec=1, p=None, algorithm='auto'):

    if not ring.is_Ring(base_ring):
        raise TypeError(f'Base ring {base_ring} must be a ring.')

    if base_ring not in _CommutativeRings:
        raise TypeError(f'Cannot create Ring of Witt Vectors over {base_ring}: {base_ring} is not a commutative ring.')

    char = base_ring.characteristic()
    prime = None
    if p is None:
        if char not in _Primes:
            raise ValueError(f'Cannot create Ring of Witt Vectors over {base_ring}: {base_ring} has non-prime characteristic and no prime was supplied.')
        else:
            prime = char
    elif p in _Primes:
        prime = p
    else:
        raise ValueError(f'Cannot create Ring of {p}-Witt Vectors: {p} is not a prime.')

    if algorithm is None:
        algorithm = 'none'
    elif algorithm not in ['none', 'auto', 'standard', 'finotti']:
        raise ValueError(f"'{algorithm}' is not a valid algorithm. It must be one of 'none', 'auto', 'standard', or 'finotti'.")

    if prime == char: # p-typical
        if base_ring.is_field() and base_ring.is_finite():
            if not is_FiniteField(base_ring):
                # This means we have something like Zmod(p)
                base_ring = base_ring.field()

            # TODO: document that this ignores the choice of algorithm
            return WittRing_finite_field(base_ring, prec, prime,
                                         category=_CommutativeRings)
        else:
            if algorithm == 'auto':
                algorithm = 'finotti'
            return WittRing_p_typical(base_ring, prec, prime,
                                      algorithm=algorithm,
                                      category=_CommutativeRings)
    else: # non-p-typical
        if algorithm == 'finotti':
            raise ValueError(f"The 'finotti' algorithm only works for p-typical Witt Rings.")
        if base_ring(prime).is_unit():
            # TODO: document that this ignores the choice of algorithm
            return WittRing_p_invertible(base_ring, prec, prime,
                                         category=_CommutativeRings)
        elif (
            isinstance(base_ring, IntegerModRing) and \
            p**(ZZ.valuation(p)(char)) == char
        ):
            return WittRing_integers_mod_power_of_p(
                base_ring, prec, prime,
                category=_CommutativeRings
            )
        else:
            if algorithm == 'auto':
                algorithm = 'standard'
            return WittRing_non_p_typical(base_ring, prec, prime,
                                          algorithm=algorithm,
                                          category=_CommutativeRings)
