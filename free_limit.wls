(* Continuum free-fermion check. Requires FeynCalc in the Wolfram kernel. *)
Needs["FeynCalc`"];
(* Normalized Dirac trace, with tr(1)=2 for a two-component spinor. *)
kinetic = DiracSimplify[DiracTrace[GSD[p].GSD[p], TraceOfOne -> 2]]/2;
invariant = Simplify[FCI[kinetic] - FCI[SPD[p,p]]];
If[!TrueQ[invariant == 0], Print[invariant]; Exit[1]];
(* The squared staggered free band contains sin(k a)^2/a^2. *)
Print[Series[m^2 + Sin[k a]^2/a^2, {a,0,2}]];
