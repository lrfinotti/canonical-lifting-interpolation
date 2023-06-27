# Canonical Lifting by Interpolation

*The code here was forked from Jacob Dennerlein's [Canonical Lifting Comparison](https://github.com/lrfinotti/canonical-lifting-interpolation)*.*

We have versions that computes the formulas for canonical liftings and corresponding elliptic Teichmüller lift for a given prime $p \geq 5$.  The second coordinates are proved, but from the third and later coordinates it depends on a conjecture on the powers that appear on denominators (made explicit in Jacob's Ph.D. thesis).

Code is provided for both [Sage](https://www.sagemath.org/) and [Magma](http://magma.maths.usyd.edu.au/magma/).  The former is free and open source, and thus widely available, but the latter is much faster.

The interpolation code was implemented by Jacob Dennerlein.  On the other hand, the main advantage of doing the computation through interpolation is that possibility of doing it in parallel, so we provide parallel versions of the code here as well.

*The parallelization of the Magma code could not have been done without the help of Geoff Bailey from the University of Sidney.*


## Sage Version

This code was done from scratch by Jacob Dennerlein.  Hopefully the code included in the `sage/witt/witt.sage` file (for computations with Witt vectors) will be included in a future version of Sage.

The file `sage/witt/canonical_lifting.sage` provides the function `canonical_lifting_odd_char` to compute canonical liftings.  Jacob's original code would detect if Magma is installed and used to solve systems, speeding up the code considerably.  I had to disable this feature, as it would not work in parallel.  (If you do have Magma installed, you should use Magma's native code below.)  This function can be used to obtain formulas, and it is used in the interpolation version, to compute several liftings over finite fields.

It also provides the function `canonical_lifting_odd_char_interpolate_formulas_parallel(p, prec, verbose=False, nproc=8)`, which does the interpolation in parallel (fastest method).

The file `sage/test_all.sage` tests the three methods ("usual" (symbolic computation), interpolation, and interpolation in parallel).  Edit the value of the prime, precision, and number of processes as you wish.

**To do:** The code can be improved by echlonizing matrices as we construct the systems in the interpolation.  (Already done in Magma.)


## Magma Version

This requires the code from [Witt Vectors and Canonical Liftings](https://github.com/lrfinotti/witt), which provides code for faster computation with Witt vectors.

**The parallel version only works for the second coordinate!**

Running code in parallel in Magma is a lot harder than in Sage.  (See the [manual entry for parallelization](https://magma.maths.usyd.edu.au/magma/handbook/text/70).)

Due to a bug in Magma version V2.28-1, we could not use the function `canonical_lifting_odd_char_one_step` from `magma/canonical_lifting.magma` in the parallelized version.  So we used `lift` from `lift.m` from  [Witt Vectors and Canonical Liftings](https://github.com/lrfinotti/witt).

**The path to the code from  [Witt Vectors and Canonical Liftings](https://github.com/lrfinotti/witt) must be entered at the top of `magma/canonical_lifting_interpolate.magma`.**  (It currently points to `~\code\witt`.)

The problem with parallelization in Magma is that I could not make it into a function.  I've made it a [BASH](https://en.wikipedia.org/wiki/Bash_(Unix_shell)) script, to spawn the workers to do the parallel work.  The script is `magma/interpoate_parallel.sh`.  Make is executable (e.g. `chmod 755 magma/interpolate_parellel.sh`) and call it with two arguments: the prime and the number of processes.  By default, it writes the result in the file `results_p_2.txt`, where `p` is the prime that was passed to the script.

The file `magma/test_all.sh` is another BASH script that compares the times and results of all three methods.  Pass a prime and the number of cores for the parallelized interpolation.
