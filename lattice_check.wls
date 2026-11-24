(* Independent four-site check. Run with wolframscript -file lattice_check.wls *)
n = 4;
op[a_, j_] := KroneckerProduct @@ Table[If[k == j, PauliMatrix[a], IdentityMatrix[2]], {k, n}];
id = IdentityMatrix[2^n];
z = Table[op[3, j], {j, n}];
q = Table[(z[[j]] + (-1)^j id)/2, {j, n}];
flux = Accumulate[q][[1 ;; n - 1]];
density = Total[Table[(id + (-1)^j z[[j]])/(2 n), {j, n}]];
h = Sum[(op[1,j].op[1,j+1] + op[2,j].op[2,j+1])/2, {j,n-1}] +
    Sum[(-1)^j z[[j]]/4, {j,n}] + Total[(#.#)& /@ flux];
vac = UnitVector[2^n, 6]; (* |0101>, site 1 most significant *)
comm[a_, b_] := a.b - b.a;
curvature = -vac.comm[h, comm[h, density]].vac;
If[curvature =!= 3, Exit[1]];
file = FileNameJoin[{DirectoryName[$InputFileName], "exact_dynamics.csv"}];
rows = Rest[Import[file,"CSV"]];
error = Max[Table[
    psi = MatrixExp[-I row[[1]] N[h]].vac;
    Max[Abs[(Re[Conjugate[psi].#.psi]& /@ Prepend[flux,density]) - Rest[row]]],
    {row,rows}]];
Print[{"density curvature",curvature,"maximum discrepancy",error}];
If[!TrueQ[error < 10^-10], Exit[1]];
