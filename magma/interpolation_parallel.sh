#! /bin/bash

# USAGE: interpolation_parallel.sh prime no_cores

# prime
p="$1"

# number of cores
NO_CORES="$2"


# avoid file conflicts using time stamp
TIMESTAMP=$(date +%s)

# main computaion file (manager)
MAIN="main_${TIMESTAMP}.magma"

# worker for lifting in parallel
WORKER_L="worker_lift_${TIMESTAMP}.magma"

# worker for solving systems in parallel
WORKER_S="worker_system_${TIMESTAMP}.magma"


# write main file
cat <<EOF > $MAIN
Write("lock_file", "");
t := Realtime();
HOST := "localhost";
PORT := 12345;
S := Socket(: LocalHost := HOST, LocalPort := PORT);

load "canonical_lifting_interpolate.magma";

print "Get interpolation values...";
it := Realtime();
v := get_interp_values($p);
print "Done!  Time: ", Realtime(it);
print "";


print "Lift (in parallel)...";
it := Realtime();
lift_data := DistributedManager(S, v : initial_results := [* *]);
print "Done!  Time: ", Realtime(it);
print "";

prec := 2;
n := 1;
p := ${p};
h_power := Integers()!(n*p^(n-1) + (n-1)*p^(n-2));
denom_weight := (p-1)*h_power;


// generate input
// an:
an_numerator_weight := 4*p^n + denom_weight;
an_numerator_powers := get_ab_powers(an_numerator_weight);
v := [* [* lift_data, an_numerator_powers, "a", 0 *] *];

// bn:
bn_numerator_weight := 6*p^n + denom_weight;
bn_numerator_powers := get_ab_powers(bn_numerator_weight);
Append(~v,  [* lift_data, bn_numerator_powers, "b", 0 *]);

// ci's:
Fn_degree_bound := ((n+2)*p^n - n*p^(n-1)) div 2;
num_cis := Fn_degree_bound div p + 1;
for i in [1..num_cis] do
    ci_numerator_weight := 2*p^n - 2*(i-1)*p + denom_weight;
    powers := get_ab_powers(ci_numerator_weight);
    Append(~v, [* lift_data, powers, "c", i *]);

end for;

// di's:
num_dis := ((n+3)*p^n - n*p^(n-1) - 3) div 2 + 1;
for i in [1..num_dis] do
    di_numerator_weight := 3*p^n - 3 - 2*(i-1) + denom_weight;
    powers := get_ab_powers(di_numerator_weight);
    Append(~v, [* lift_data, powers, "d", i *]);
end for;

// solve in parallel
PORT2 := 12346;
S2 := Socket(: LocalHost := HOST, LocalPort := PORT2);

print "Solve systems (in parallel)...";
it := Realtime();

print "Remove lock file...";
%! rm -f lock_file
print "Done!";

solutions := DistributedManager(S2, v : initial_results := [* *]);

print "Done!  Time: ", Realtime(it);
print "";

// Write("solutions.txt", solutions);


// solutions := [* *];
// for s in systems do Append(~solutions, solve_system(s)); end for;

print "Find formulas...";
it := Realtime();
// _ := find_formulas(${p}, solutions);
va, vb, vF, vG := find_formulas(${p}, solutions);  Write("results_${p}_2.txt", [* va, vb, vF, vG *] : Overwrite:=true);
// print find_formulas(${p}, solutions);
print "Done!  Time: ", Realtime(it);
print "";


print "";
print "#################################################";
print "Total time: ", Realtime(t);
print "#################################################";

EOF


# write worked file for lifting
cat <<EOF > $WORKER_L
HOST := "localhost";
PORT := 12345;

load "canonical_lifting_interpolate.magma";

DistributedWorker(HOST, PORT, lift_each);
// DistributedWorker(HOST, PORT, f);
quit;
EOF

magma -b < $MAIN &

# avoid running too soon
sleep 0.2

# run the correct number of workers
for (( i=1; i<="$NO_CORES"; i++ ))
do
    magma -b < $WORKER_L 1> /dev/null &
    # magma -b < $WORKER_L &
done

# main/manager will remove the lock file when ready to solve systems
until [ `ls -1 lock_file 2>/dev/null | wc -l ` -eq 0 ]; do sleep 0.1; done

# avoid running too soon
sleep 0.2

# worker to solve the systems
cat <<EOF > $WORKER_S
HOST := "localhost";
PORT := 12346;

load "canonical_lifting_interpolate.magma";

DistributedWorker(HOST, PORT, solve_systems_parallel);
// DistributedWorker(HOST, PORT, f);
quit;
EOF

# run the correct number of workers
for (( i=1; i<="$NO_CORES"; i++ ))
do
    magma -b < $WORKER_S 1> /dev/null &
    # magma -b < $WORKER_S &
done

sleep 0.5

# clean up files
rm -f $MAIN $WORKER_L $WORKER_S
