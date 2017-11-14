(* ::Package:: *)

BeginPackage["Effmodel`"]


toTensor2D::usage::"Takes a flattened vector to a square array";

cross2D::usage::"The 2D cross product used to find the orthogonal vector";

calc$m::usage::"Finds the new unit vector direction of M";

calc$s::usage::"Finds the new unit vector direction orthogonal of M";


lambda$M::usage::"compute the stretch along M";

lambda$S::usage::"compute the stretch along S";

phi::usage::"computes the shear angle";

gamma$1::usage::"compute the invariant";

gamma$2::usage::"compute the invariant";

gamma$3::usage::"compute the invariant";


W1$fromdata::usage::"computes the first response function from some given data";

W2$fromdata::usage::"computes the second response function from some given data";

W3$fromdata::usage::"computes the third response function from some given data";


Begin["`Private`"]


toTensor2D = Function[{f},{{f[[1]],f[[2]]},{f[[3]],f[[4]]}}];

dot2D = Function[{F,M},
					{F[[1]]*M[[1]]+F[[2]]*M[[2]], \
					F[[3]]*M[[1]]+F[[4]]*M[[2]]}];

cross2D = Function[{M},{-M[[2]],M[[1]]}];

calc$m = Function[{F,M}, 
			Module[{tF, vec, m}, \
				vec = {F[[1]]*M[[1]]+F[[2]]*M[[2]], \
					F[[3]]*M[[1]]+F[[4]]*M[[2]]};
				vec/(vec.vec)
			]];

calc$s = Function[{F,M}, 
			Module[{m}, 
				m = calc$m[F,M];
				cross2D[m]
				]
		];


lambda$M = Function[{F,M},
				Module[{m},
					m = calc$m[F,M];
					m.dot2D[F,M]
					]];
					
lambda$S = Function[{F,M},
			Module[{s,S},
				S = cross2D[M];
				s = calc$s[F,M];
				s.dot2D[F,S]
			]];
			
phi = Function[{F,M},
		Module[{m,S},
			S = cross2D[M];
			m = calc$m[F,M];
			m.dot2D[F,S]
			]];
			
gamma$1 = Function[{F,M}, 
			Log[lambda$M[F,M]]];
			
gamma$2 = Function[{F,M}, 
			Log[lambda$S[F,M]]];
			
gamma$3 = Function[{F,M}, 
			phi[F,M]];


End[]


EndPackage[]
