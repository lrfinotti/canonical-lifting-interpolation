#! /bin/bash

# USAGE: test_all.sh prime no_cores

# prime
p="$1"

# number of processes
NO_CORES="$2"


TIMESTAMP=$(date +%s)

# usual
cat <<EOF >  "usual_${TIMESTAMP}.magma"
SetPath(":~/code/witt");
load "canonical_lifting.magma";

p := ${p};
prec := 2;  // for parallel, only second coordinate

F<a,b> := RationalFunctionField(GF(p), 2);
va, vb, vF, vG := canonical_lifting_odd_char(a, b, prec, false);
Write("results_${TIMESTAMP}_usual_${p}_2.txt", [* va, vb, vF, vG *] : Overwrite:=true);
// print va, vb, vF, vG;
EOF

echo "Computing the usual way..."
time magma -b <  "usual_${TIMESTAMP}.magma" 1> /dev/null
echo "Done!"

echo ""
echo ""

# convert x0->x and y0->y
sed -i -e 's/x0/x/g' -e 's/y0/y/g' "results_${TIMESTAMP}_usual_${p}_2.txt"


rm -f "usual_${TIMESTAMP}.magma"

# interpolation
cat <<EOF >  "interpol_${TIMESTAMP}.magma"
SetPath(":~/code/witt");
load "canonical_lifting_interpolate.magma";

p := ${p};
prec := 2;  // for parallel, only second coordinate

va, vb, vF, vG := canonical_lifting_odd_char_interpolate_formulas(p, prec);
Write("results_${TIMESTAMP}_interpol_${p}_2.txt", [* va, vb, vF, vG *] : Overwrite:=true);
// print va, vb, vF, vG;
EOF

echo "Computing using interpolation..."
time magma -b <  "interpol_${TIMESTAMP}.magma" 1> /dev/null
echo "Done!"

echo ""
echo ""

rm -f "interpol_${TIMESTAMP}.magma"

echo "Computing using PARALLEL interpolation..."
./interpolation_parallel.sh $p $NO_CORES

sleep 1

# COMPARE results!

cat <<EOF > "compare_${TIMESTAMP}.magma"
p := ${p};
F<a,b> := RationalFunctionField(GF(p), 2);
PR<x, y> := PolynomialRing(F, 2);

res1 := Read("results_${TIMESTAMP}_usual_${p}_2.txt");
res2 := Read("results_${TIMESTAMP}_interpol_${p}_2.txt");
res3 := Read("results_${p}_2.txt");

print (res1[1] eq res2[1]) and (res1[2] eq res2[2]) and (res1[3] eq res2[3]) and (res1[4] eq res2[4]) and
      (res1[1] eq res3[1]) and (res1[2] eq res3[2]) and (res1[3] eq res3[3]) and (res1[4] eq res3[4]);

EOF

echo "Are results equal?"

magma -b < "compare_${TIMESTAMP}.magma"

rm -f "compare_${TIMESTAMP}.magma" "results_${TIMESTAMP}_"*
