# Canonical Lifting by Interpolation

*The code here was forked from Jacob Dennerlein's [Canonical Lifting Comparison](https://github.com/lrfinotti/canonical-lifting-interpolation)*.*

We have versions that computes the formulas for canonical liftings and corresponding elliptic Teichmüller lift for a given prime $p \geq 5$.  The second coordinates are proved, but from the third and later coordinates it depends on a conjecture on the powers that appear on denominators (made explicit in Jacob's Ph.D. thesis).

Code is provided for both [Sage](https://www.sagemath.org/) and [Magma](http://magma.maths.usyd.edu.au/magma/).  The former is free and open source, and thus widely available, but the latter is considerably faster.

The interpolation code was implemented by Jacob Dennerlein.  On the other hand, the main advantage of doing the computation through interpolation is that possibility of doing it in parallel, so we provide parallel versions of the code here as well.

*The parallelization of the Magma code could not have been done without the help of Geoff Bailey from the University of Sidney.*


## Sage Version

This code was done from scratch by Jacob Dennerlein.  Hopefully the code included in the `sage/witt/witt.sage` file (for computations with Witt vectors) will be included in a future version of Sage.

The file `sage/witt/canonical_lifting.sage` provides the function `canonical_lifting_odd_char` to compute canonical liftings.  Jacob's original code would detect if Magma is installed and used to solve systems, speeding up the code considerably.  I had to disable this feature, as it would not work in parallel.  (If you do have Magma installed, you should use Magma's native code below.)  This function can be used to obtain formulas, and it is used in the interpolation version, to compute several liftings over finite fields.

It also provides the function `canonical_lifting_odd_char_interpolate_formulas_parallel(p, prec, verbose=False, nproc=8)`, which does the interpolation in parallel (fastest method).

The file `sage/test_all.sage` provides the function `test_all` that tests the three methods ("usual" (symbolic computation), interpolation, and interpolation in parallel).

**TO DO:** The code can be improved by echlonizing matrices as we construct the systems in the interpolation.  (Already done in Magma.)



## Magma Version

This requires the code from [Witt Vectors and Canonical Liftings](https://github.com/lrfinotti/witt), which provides code for faster computation with Witt vectors.

**The parallel version only works for the second coordinate!**

Running code in parallel in Magma is a lot harder than in Sage.  (See the [manual entry for parallelization](https://magma.maths.usyd.edu.au/magma/handbook/text/70).)

Due to a bug in Magma version V2.28-1, we could not use the function `canonical_lifting_odd_char_one_step` from `magma/canonical_lifting.magma` in the parallelized version.  So we used `lift` from `lift.m` from  [Witt Vectors and Canonical Liftings](https://github.com/lrfinotti/witt).

**The path to the code from  [Witt Vectors and Canonical Liftings](https://github.com/lrfinotti/witt) must be entered at the top of `magma/canonical_lifting_interpolate.magma`.**  (It currently points to `~\code\witt`.)

The problem with parallelization in Magma is that I could not make it into a function.  I've made it a [BASH](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) script, to spawn the workers to do the parallel work.  The script is `magma/interpoate_parallel.sh`.  Make is executable (e.g. `chmod 755 magma/interpolate_parellel.sh`) and call it with two arguments: the prime and the number of processes.  By default, it writes the result in the file `results_p_2.txt`, where `p` is the prime that was passed to the script.

The file `magma/test_all.sh` is another BASH script that compares the times and results of all three methods.  Pass a prime and the number of cores for the parallelized interpolation.


## Examples

### Sage Examples

We assume here you start Sage from `sage/` folder in the code folder.

Here is how we can compute two coordinates the canonical lifting and elliptic Teichüller lift of the curve:

$$E/\mathbb{F}_{11}: y^2 = x^2 + 5x + 2.$$

```
sage: load("witt/witt.sage")
sage: load("witt/canonical_lifting.sage")
sage: F = GF(11)
sage: canonical_lifting_odd_char(F(5), F(2), 11, 2)
(([5, 4], [2, 2]),
 ([x0,
   -x0^16 - 5*x0^14 - 3*x0^13 + 4*x0^12 + 4*x0^10 - x0^9 - 3*x0^8 + 3*x0^7 - 2*x0^6 + x0^5 - 3*x0^4 - 4*x0^3 + x0^2 + 5*x0 + 4],
  [y0,
   4*x0^20*y0 - x0^18*y0 - 3*x0^17*y0 - x0^16*y0 + 2*x0^15*y0 - 3*x0^14*y0 - x0^13*y0 + 2*x0^12*y0 - x0^10*y0 + 2*x0^9*y0 - 2*x0^8*y0 + 2*x0^7*y0 - x0^6*y0 + 3*x0^5*y0 - 5*x0^4*y0 + 2*x0^3*y0 - 5*x0*y0]))
```

To obtain the formulas for two coordinates in characteristic 5, where the curve has Weierstrass equation $y^2 = x^2 + ax + b$, the "usual" way:

```
sage: load("witt/witt.sage")
sage: load("witt/canonical_lifting.sage")
sage: R.<a, b> = PolynomialRing(GF(5), 2)
sage: F = R.fraction_field()
sage: canonical_lifting_odd_char(F(a), F(b), 5, 2)
(([a, (2*a^3*b^2 + 2*b^4)/(2*a)], [b, -a^6*b + a^3*b^3 + b^5]),
 ([x0, -1/a*x0^7 + (-b)/a*x0^4 + a*x0^3 + (-2*b)*x0^2 + (-2*b^2)/a*x0 + (a*b)],
  [y0,
   1/a*x0^8*y0 + 2*x0^6*y0 + (-2*b)/a*x0^5*y0 + (2*a)*x0^4*y0 + (-2*b)*x0^3*y0 + (a^2)*x0^2*y0 + (-2*a*b)*x0*y0 + (-2*a^3 - 2*b^2)*y0]))
```

Or, we can do it with interpolation:

```
sage: load("witt/witt.sage")
sage: load("witt/canonical_lifting_interpolate.sage")
sage: canonical_lifting_odd_char_interpolate_formulas(5, 2)
(([a, (2*a^3*b^2 + 2*b^4)/(2*a)], [b, -a^6*b + a^3*b^3 + b^5]),
 ([x, -1/a*x^7 + (-b)/a*x^4 + a*x^3 + (-2*b)*x^2 + (-2*b^2)/a*x + (a*b)],
  [y,
   1/a*x^8*y + 2*x^6*y + (-2*b)/a*x^5*y + (2*a)*x^4*y + (-2*b)*x^3*y + (a^2)*x^2*y + (-2*a*b)*x*y + (-2*a^3 - 2*b^2)*y]))
```

To obtain the formulas by doing the interpolation in parallel, for $p=37$ and using 20 cores:

```
sage: load("witt/witt.sage")
sage: load("witt/canonical_lifting_interpolate.sage")
sage: canonical_lifting_odd_char_interpolate_formulas_parallel(37, 2, n_proc=20)
computing canonical liftings!
done
solving the systems!
done!
(([a,
   (16*a^43*b^2 - 3*a^40*b^4 + 11*a^37*b^6 + 14*a^34*b^8 - 13*a^31*b^10 + 4*a^28*b^12 - 10*a^25*b^14 + a^22*b^16 + 2*a^19*b^18 + 2*a^16*b^20 - 6*a^13*b^22 - 11*a^10*b^24 + 4*a^7*b^26 - 13*a^4*b^28 + 11*a*b^30)/(2*a^9 + 6*a^6*b^2 - 8*a^3*b^4 - 10*b^6)],
  [b,
   (-18*a^63*b + 4*a^60*b^3 - 15*a^57*b^5 + 16*a^54*b^7 + 16*a^51*b^9 + 4*a^48*b^11 - 17*a^45*b^13 - 16*a^42*b^15 - 3*a^39*b^17 + 2*a^36*b^19 - 18*a^33*b^21 - 16*a^30*b^23 + 3*a^27*b^25 - 9*a^24*b^27 - 15*a^21*b^29 + 6*a^18*b^31 - 9*a^15*b^33 - 17*a^12*b^35 - 6*a^9*b^37 + 18*a^6*b^39 - a^3*b^41 + 17*b^43)/(2*a^9 + 6*a^6*b^2 - 8*a^3*b^4 - 10*b^6)]),
 ([x,
   -1/(a^9 + 3*a^6*b^2 - 4*a^3*b^4 - 5*b^6)*x^55 + (-11*a)/(a^9 + 3*a^6*b^2 - 4*a^3*b^4 - 5*b^6)*x^53 + (8*b)/(a^9 + 3*a^6*b^2 - 4*a^3*b^4 - 5*b^6)*x^52 + (-17*a^2)/(a^9 + 3*a^6*b^2 - 4*a^3*b^4 - 5*b^6)*x^51 + (-11*a*b)/(a^9 + 3*a^6*b^2 - 4*a^3*b^4 - 5*b^6)*x^50 + (-3*a^3 + 11*b^2)/(a^9 + 3*a^6*b^2 - 4*a^3*b^4 - 5*b^6)*x^49 + (7*a^2*b)/(a^9 + 3*a

[output snipped]
```

To run a comparison between the three methods, verifying we obtain the same result, say for $p=37$ and using 20 cores for the parallelized version:

```
sage: test_all(37, 2, 20)
Computing the "usual way":
Time: CPU 24.52 s, Wall: 24.55 s
Done!


Computing with interpolation:
Time: CPU 4.00 s, Wall: 4.01 s
Done!


Computing with interpolation in parallel (20 cores):
computing canonical liftings!
done
solving the systems!
done!
Time: CPU 0.78 s, Wall: 1.22 s
Done!


Check if same results:
True
```

We can also compute the formulas for more than two coordinates, but there is a possibility that it might fail, as uses some conjectures (for the degrees in the denominators of the formulas) as assumption:

```
sage: canonical_lifting_odd_char_interpolate_formulas_parallel(11, 3, false, 20)
computing canonical liftings!
done
solving the systems!
done!
computing canonical liftings!
done
solving the systems!
done!
(([a,
   (2*a^12 + a^9*b^2 - 2*a^6*b^4 - 2*a^3*b^6 - 5*b^8)/(-2*a),
   (4*a^177 + 2*a^162*b^10 + 3*a^159*b^12 - 5*a^156*b^14 - 4*a^153*b^16 + 3*a^150*b^18 + 4*a^147*b^20 + 5*a^144*b^22 - 2*a^141*b^24 + 2*a^138*b^26 - 5*a^135*b^28 - 2*a^132*b^30 - a^129*b^32 + a^126*b^34 + 5*a^120*b^38 + 2*a^117*b^40 - 5*a^111*b^44 - 4*a^108*b^46 - 2*a^105*b^48 + 5*a^102*b^50 - 3*a^99*b^52 + 2*a^96*b^54 - 4*a^93*b^56 + 5*a^90*b^58 + 2*a^87*b^60 - 2*a^84*b^62 + 3*a^81*b^64 + 4*a^78*b^66 - 3*a^75*b^68 - 4*a^72*b^70 - 4*a^69*b^72 - 4*a^66*b^74 - 3*a^63*b^76 + 4*a^60*b^78 + 3*a^57*b^80 - 2*a^54*b^82 + 2*a^51*b^84 - 2*a^48*b^86 - 3*a^45*b^88 + 2*a^36*b^94 + a^33*b^96 - 5*a^30*b^98 + 2*a^27*b^100 - 2*a^24*b^102 - 4*a^21*b^104 - a^18*b^106 - a^15*b^108 - a^12*b^110 + a^9*b^112 - 2*a^6*b^114 - 2*a^3*b^116 - 5*b^118)/(3*a^23*b^22)],
  [b,
   (-5*a^18 + 2*a^15*b^2 + a^12*b^4 + 2*a^9*b^6 - a^6*b^8 - 3*a^3*b^10 + 2*b^12)/(-2*b),
   (-5*a^225 + 2*a^222*b^2 + a^219*b^4 + 2*a^216*b^6 - a^213*b^8 - 3*a^210*b^10 - a^207*b^12 - 4*a^204*b^14 + 4*a^201*b^16 + 5*a^198*b^18 - 4*a^195*b^20 - 4*a^189*b^24 - 2*a^186*b^26 - a^183*b^28 - a^180*b^30 + 5*a^177*b^32 + a^174*b^34 - 3*a^171*b^36 + 2*a^168*b^38 - a^165*b^40 + 5*a^162*b^42 + 4*a^159*b^44 - 2*a^156*b^46 + a^153*b^48 - 4*a^150*b^50 + 2*a^147*b^52 - a^144*b^54 + 4*a^141*b^56 + 5*a^138*b^58 - a^135*b^60 + 2*a^132*b^62 - 3*a^126*b^66 + 2*a^123*b^68 - 2*a^120*b^70 + 5*a^117*b^72 + 2*a^114*b^74 + 4*a^111*b^76 + 2*a^108*b^78 - 2*a^105*b^80 + 5*a^102*b^82 + 4*a^99*b^84 - 4*a^96*b^86 - 4*a^93*b^88 + 5*a^90*b^90 + 5*a^87*b^92 - 2*a^81*b^96 + 3*a^78*b^98 + 3*a^75*b^100 - 3*a^72*b^102 + 2*a^69*b^104 - 5*a^63*b^108 - 5*a^60*b^110 + 2*a^57*b^112 + a^54*b^114 + 5*a^51*b^116 - 5*a^48*b^118 - 2*a^45*b^120 + 2*a^42*b^122 - 4*a^36*b^126 + 2*a^33*b^128 - 2*a^30*b^130 - 4*a^27*b^132 - 2*a^24*b^134 + 2*a^21*b^136 - 2*a^18*b^138 - 3*a^15*b^140 - 5*a^12*b^142 - 3*a^6*b^146 - 4*a^3*b^148 - b^150)/(3*a^9*b^23)]),
 ([x,
   1/(a*b)*x^16 + 1/b*x^14 - 4/a*x^13 + (-5*a)/b*x^12 + (5*a^3 + 5*b^2)/(a*b)*x^10 + (2*a)*x^9 + (-a^3 + 5*b^2)/b*x^8 + (-3*a^3 + 4*b^2)/a*x^7 + (-a^4 + 3*a*b^2)/b*x^6 + (5*a^3 - 2*b^2)*x^5 + (-4*a^3*b - 2*b^3)/a*x^4 + (2*a*b^2)*x^3 + (-4*b^3)*x^2 + (5*b^4)/a*x + (3*a^7 + a^4*b^2 - 5*a*b^4)/b,

[output snipped]
```


### Magma Examples

We assume here you start Magma from the `magma` folder in the code folder.  We also assume the code from [Witt Vectors and Canonical Liftings](https://github.com/lrfinotti/witt) is in `~/code/witt/`.

Here is how we can compute two coordinates the canonical lifting and elliptic Teichüller lift of the curve:

$$E/\mathbb{F}_{11}: y^2 = x^2 + 5x + 2.$$

```
> load "canonical_lifting.magma";
Loading "canonical_lifting.magma"
Loading "/home/finotti/code/witt/gt.m"
Loading "/home/finotti/code/witt/witt.m"
Loading "/home/finotti/code/witt/etas.m"
> F := GF(11);
> canonical_lifting_odd_char(F!5, F!2, 2, false);
[
    5,
    4
]
[
    2,
    2
]
[
    x0,
    10*x0^16 + 6*x0^14 + 8*x0^13 + 4*x0^12 + 4*x0^10 + 10*x0^9 + 8*x0^8 + 3*x0^7 + 9*x0^6 + x0^5 + 8*x0^4 + 7*x0^3 +
        x0^2 + 5*x0 + 4
]
[
    y0,
    4*x0^20*y0 + 10*x0^18*y0 + 8*x0^17*y0 + 10*x0^16*y0 + 2*x0^15*y0 + 8*x0^14*y0 + 10*x0^13*y0 + 2*x0^12*y0 +
        10*x0^10*y0 + 2*x0^9*y0 + 9*x0^8*y0 + 2*x0^7*y0 + 10*x0^6*y0 + 3*x0^5*y0 + 6*x0^4*y0 + 2*x0^3*y0 + 6*x0*y0
]
```


To obtain the formulas for two coordinates in characteristic 5, where the curve has Weierstrass equation $y^2 = x^2 + ax + b$, the "usual" way:

```
> load "canonical_lifting.magma";
Loading "canonical_lifting.magma"
Loading "/home/finotti/code/witt/gt.m"
Loading "/home/finotti/code/witt/witt.m"
Loading "/home/finotti/code/witt/etas.m"
> F<a, b> := RationalFunctionField(GF(5), 2);
> canonical_lifting_odd_char(a, b, 2, false);
[
    a,
    (a^3*b^2 + b^4)/a
]
[
    b,
    4*a^6*b + a^3*b^3 + b^5
]
[
    x0,
    4/a*x0^7 + 4*b/a*x0^4 + a*x0^3 + 3*b*x0^2 + 3*b^2/a*x0 + a*b
]
[
    y0,
    1/a*x0^8*y0 + 2*x0^6*y0 + 3*b/a*x0^5*y0 + 2*a*x0^4*y0 + 3*b*x0^3*y0 + a^2*x0^2*y0 + 3*a*b*x0*y0 + (3*a^3 +
        3*b^2)*y0
]
```


We could also use the code from [Witt Vectors and Canonical Liftings](https://github.com/lrfinotti/witt):

```
> SetPath(":~/code/witt");
> load "lift.m";
Loading "/home/finotti/code/witt/lift.m"
Loading "/home/finotti/code/witt/gt.m"
Loading "/home/finotti/code/witt/witt.m"
Loading "/home/finotti/code/witt/etas.m"
> F<a, b> := RationalFunctionField(GF(5), 2);
> lift(a, b, 1);
[
    a,
    (a^3*b^2 + b^4)/a
]
[
    b,
    4*a^6*b + a^3*b^3 + b^5
]
[
    x0,
    4/a*x0^7 + 4*b/a*x0^4 + a*x0^3 + 3*b*x0^2 + 3*b^2/a*x0 + a*b
]
[
    1,
    1/a*x0^8 + 2*x0^6 + 3*b/a*x0^5 + 2*a*x0^4 + 3*b*x0^3 + a^2*x0^2 + 3*a*b*x0 + 3*a^3 + 3*b^2
]
```

(Note we use `1` instead of `2` for two coordinates with `lift`!)

Or, we can do it with interpolation:

```
> load "canonical_lifting_interpolate.magma";
Loading "canonical_lifting_interpolate.magma"
Loading "canonical_lifting.magma"
Loading "/home/finotti/code/witt/gt.m"
Loading "/home/finotti/code/witt/witt.m"
Loading "/home/finotti/code/witt/etas.m"
Loading "/home/finotti/code/witt/lift.m"
Loading "/home/finotti/code/witt/gt.m"
Loading "/home/finotti/code/witt/witt.m"
Loading "/home/finotti/code/witt/etas.m"
> canonical_lifting_odd_char_interpolate_formulas(5, 2);
[
    a,
    (a^3*b^2 + b^4)/a
]
[
    b,
    4*a^6*b + a^3*b^3 + b^5
]
[
    x,
    4/a*x^7 + 4*b/a*x^4 + a*x^3 + 3*b*x^2 + 3*b^2/a*x + a*b
]
[
    y,
    1/a*x^8*y + 2*x^6*y + 3*b/a*x^5*y + 2*a*x^4*y + 3*b*x^3*y + a^2*x^2*y + 3*a*b*x*y + (3*a^3 + 3*b^2)*y
]
```

To obtain the formulas by doing the interpolation in parallel, for $p=101$ and using 20 cores, we can call, *from the shell* (and not from Magma):

```
$ interpolation_parallel.sh 101 20
Loading "canonical_lifting_interpolate.magma"
Loading "canonical_lifting.magma"
Loading "/home/finotti/code/witt/gt.m"
Loading "/home/finotti/code/witt/witt.m"
Loading "/home/finotti/code/witt/etas.m"
Loading "/home/finotti/code/witt/lift.m"
Loading "/home/finotti/code/witt/gt.m"
Loading "/home/finotti/code/witt/witt.m"
Loading "/home/finotti/code/witt/etas.m"
Get interpolation values...
Done!  Time:  0.030

Lift (in parallel)...
Done!  Time:  0.380

Solve systems (in parallel)...
Remove lock file...

Find formulas...
Done!  Time:  0.090


#################################################
Total time:  1.370
#################################################
```

 By default, the result will be saved in `results_p_2.txt`, where `p` represents the prime:

 ```
 $ cat results_101_2.txt
[*
    [
        a,
        (86*a^123*b^2 + 53*a^120*b^4 + 52*a^117*b^6 + 48*a^114*b^8 +
            42*a^111*b^10 + 40*a^108*b^12 + 25*a^105*b^14 + 70*a^102*b^16 +
            76*a^99*b^18 + 57*a^96*b^20 + 23*a^93*b^22 + 92*a^90*b^24 +
            81*a^87*b^26 + 84*a^84*b^28 + 66*a^81*b^30 + 75*a^78*b^32 +
            67*a^75*b^34 + 2*a^72*b^36 + 97*a^69*b^38 + 19*a^66*b^40 +
            79*a^63*b^42 + 47*a^60*b^44 + 97*a^57*b^46 + 75*a^54*b^48 +
            a^51*b^50 + a^48*b^52 + 98*a^45*b^54 + 13*a^42*b^56 + 97*a^39*b^58 +
            42*a^36*b^60 + 7*a^33*b^62 + 45*a^30*b^64 + 15*a^27*b^66 +
            48*a^24*b^68 + 73*a^21*b^70 + 51*a^18*b^72 + 3*a^15*b^74 +
            81*a^12*b^76 + 89*a^9*b^78 + 83*a^6*b^80 + 84*a^3*b^82 + b^84)/(a^25
            + 9*a^22*b^2 + 21*a^19*b^4 + 9*a^16*b^6 + 43*a^13*b^8 + 39*a^10*b^10
            + 35*a^7*b^12 + 69*a^4*b^14 + 25*a*b^16)
    ],
    [
        b,
        (76*a^174*b + 57*a^171*b^3 + 78*a^168*b^5 + 83*a^165*b^7 + 90*a^162*b^9
            + 64*a^159*b^11 + 7*a^156*b^13 + 96*a^153*b^15 + 100*a^150*b^17 +
            2*a^147*b^19 + 73*a^144*b^21 + 55*a^141*b^23 + 54*a^138*b^25 +
            64*a^135*b^27 + 46*a^132*b^29 + 75*a^129*b^31 + 99*a^126*b^33 +
            89*a^123*b^35 + 71*a^120*b^37 + 93*a^117*b^39 + 70*a^114*b^41 +
            54*a^111*b^43 + 92*a^108*b^45 + 96*a^105*b^47 + 80*a^102*b^49 +
            a^99*b^51 + 59*a^96*b^53 + 59*a^93*b^55 + 66*a^90*b^57 +
            98*a^87*b^59 + 43*a^84*b^61 + 90*a^81*b^63 + 54*a^78*b^65 +
            90*a^75*b^67 + 14*a^72*b^69 + 93*a^69*b^71 + 94*a^66*b^73 +
            52*a^63*b^75 + 95*a^60*b^77 + 28*a^57*b^79 + 66*a^54*b^81 +
            56*a^51*b^83 + 79*a^48*b^85 + 12*a^45*b^87 + 4*a^42*b^89 +
            16*a^39*b^91 + 71*a^36*b^93 + 84*a^33*b^95 + 23*a^30*b^97 +
            11*a^27*b^99 + 27*a^24*b^101 + 90*a^21*b^103 + 33*a^18*b^105 +
            9*a^15*b^107 + 89*a^12*b^109 + 42*a^9*b^111 + 72*a^6*b^113 +
            66*a^3*b^115 + 55*b^117)/(a^24 + 9*a^21*b^2 + 21*a^18*b^4 +
            9*a^15*b^6 + 43*a^12*b^8 + 39*a^9*b^10 + 35*a^6*b^12 + 69*a^3*b^14 +
            25*b^16)
    ],
    [
        x,
        100/(a^25 + 9*a^22*b^2 + 21*a^19*b^4 + 9*a^16*b^6 + 43*a^13*b^8 +
            39*a^10*b^10 + 35*a^7*b^12 + 69*a^4*b^14 + 25*a*b^16)*x^151 +
            91/(a^24 + 9*a^21*b^2 + 21*a^18*b^4 + 9*a^15*b^6 + 43*a^12*b^8 +
            39*a^9*b^10 + 35*a^6*b^12 + 69*a^3*b^14 + 25*b^16)*x^149 +
            65*b/(a^25 + 9*a^22*b^2 + 21*a^19*b^4 + 9*a^16*b^6 + 43*a^13*b^8 +
            39*a^10*b^10 + 35*a^7*b^12 + 69*a^4*b^14 + 25*a*b^16)*x^148 +
            21*a/(a^24 + 9*a^21*b^2 + 21*a^18*b^4 + 9*a^15*b^6 + 43*a^12*b^8 +

[output snipped -- VERY long]
 ```


To run a comparison between the three methods, verifying we obtain the same result, say for we run *from the shell* (and not from Magma), say for $p=101$ again, and using 20 cores for the parallelized version:


```
$ test_all.sh 101 20
Computing the usual way...

real    0m5.457s
user    0m5.419s
sys     0m0.017s
Done!


Computing using interpolation...

real    0m49.155s
user    0m48.583s
sys     0m0.353s
Done!


Computing using PARALLEL interpolation...
Loading "canonical_lifting_interpolate.magma"
Loading "canonical_lifting.magma"
Loading "/home/finotti/code/witt/gt.m"
Loading "/home/finotti/code/witt/witt.m"
Loading "/home/finotti/code/witt/etas.m"
Loading "/home/finotti/code/witt/lift.m"
Loading "/home/finotti/code/witt/gt.m"
Loading "/home/finotti/code/witt/witt.m"
Loading "/home/finotti/code/witt/etas.m"
Get interpolation values...
Done!  Time:  0.030

Lift (in parallel)...
Done!  Time:  0.380

Solve systems (in parallel)...
Remove lock file...
Done!
Done!  Time:  0.900

Find formulas...
Done!  Time:  0.090


#################################################
Total time:  1.420
#################################################
Are results equal?
true
```
