(* ::Package:: *)

BeginPackage["Effmodel`"]


toTensor2D::usage::"Takes a flattened vector to a square array";

cross2D::usage::"The 2D cross product used to find the orthogonal vector";

dot2D::usage::"The 2D dot product used to find the vector";

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


responsefunctions::usage::"computes the response functions and puts them in an array";

strainenergy::usage::"computes the strain energy given the parameters, and strain";


(* ::Section:: *)
(*Begin Function Definitions*)


Begin["`Private`"];


rtOPTs = {"CatchMachineOverflow"->False, "CatchMachineUnderflow"->False,
	"CatchMachineIntegerOverflow"->False, "CompareWithTolerance"->True,
	"EvaluateSymbolically"->False,"RuntimeErrorHandler"->Evaluate,"WarningMessages"->True};


cpOPTs = {CompilationOptions->{"ExpressionOptimization"->True, "InlineExternalDefinitions"->True,"InlineCompiledFunctions"->True},
	RuntimeAttributes->{Listable}, 
	Parallelization->True,
	RuntimeOptions->rtOPTs};


(* ::Section:: *)
(*These are the Utility functions*)


toTensor2D = Compile[{{f,_Real,1}},{{f[[1]],f[[2]]},{f[[3]],f[[4]]}},Evaluate@cpOPTs];

dot2D = Compile[{{F,_Real,1},{M,_Real,1}},
					{F[[1]]*M[[1]]+F[[2]]*M[[2]], \
					F[[3]]*M[[1]]+F[[4]]*M[[2]]},Evaluate@cpOPTs];

cross2D = Compile[{{M,_Real,1}},{-M[[2]],M[[1]]},Evaluate@cpOPTs];

calc$m = Compile[{{F,_Real,1},{M,_Real,1}}, 
			Module[{tF, vec, m}, \
				vec = {F[[1]]*M[[1]]+F[[2]]*M[[2]], \
					F[[3]]*M[[1]]+F[[4]]*M[[2]]};
				vec/Sqrt[vec.vec]
			],Evaluate@cpOPTs];

calc$s = Compile[{{F,_Real,1},{M,_Real,1}}, 
			Module[{m}, 
				m = calc$m[F,M];
				cross2D[m]
				]
		,Evaluate@cpOPTs];


(* ::Subsubsection:: *)
(*This defines the invariants*)


lambda$M = Compile[{{F,_Real,1},{M,_Real,1}},
				Module[{m},
					m = calc$m[F,M];
					m.dot2D[F,M]
					],Evaluate@cpOPTs];
					
lambda$S = Compile[{{F,_Real,1},{M,_Real,1}},
			Module[{s,S},
				S = cross2D[M];
				s = calc$s[F,M];
				s.dot2D[F,S]
			],Evaluate@cpOPTs];
			
phi = Compile[{{F,_Real,1},{M,_Real,1}},
		Module[{m,S},
			S = cross2D[M];
			m = calc$m[F,M];
			m.dot2D[F,S] / m.dot2D[F,M]
			],Evaluate@cpOPTs];
			
gamma$1 = Compile[{{F,_Real,1},{M,_Real,1}}, 
			Log[lambda$M[F,M]],Evaluate@cpOPTs];
			
gamma$2 = Compile[{{F,_Real,1},{M,_Real,1}}, 
			Log[lambda$S[F,M]],Evaluate@cpOPTs];
			
gamma$3 = Compile[{{F,_Real,1},{M,_Real,1}}, 
			phi[F,M],Evaluate@cpOPTs];
			
$gamma = Compile[{{F,_Real,1},{M,_Real,1}},
			Block[{S, m, s, FS, lambdaM},
				S = cross2D[M];
				m = calc$m[F,M];
				s = cross2D[m];
				FS = dot2D[F,S];
				lambdaM = m.dot2D[F,M];
				{Log[lambdaM], Log[s.FS], m.FS/lambdaM}
			]
			,Evaluate@cpOPTs];


(* ::Subsubsection:: *)
(*This defines Q and its derivatives*)


$Q = Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y3s},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			{c[[1]]*y1s, c[[2]]*y2s, c[[3]]*y3s, \
				c[[4]]*y1s*y2s, c[[5]]*y1s*y3s, c[[6]]*y2s*y3s, \
				c[[7]]*y1s*y2s*y3s}
		],Evaluate@cpOPTs];
		
$Q1 = Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y3s},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			2.0*{c[[1]]*y[[1]], 0.0, 0.0, \
				c[[4]]*y[[1]]*y2s, c[[5]]*y[[1]]*y3s, 0.0, \
				c[[7]]*y[[1]]*y2s*y3s}
		],Evaluate@cpOPTs];
		
$Q2 = Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y3s},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			2.0*{0.0, c[[2]]*y[[2]], 0.0, \
				c[[4]]*y1s*y[[2]], 0.0, c[[6]]*y[[2]]*y3s, \
				c[[7]]*y1s*y[[2]]*y3s}
		],Evaluate@cpOPTs];
		
$Q3 = Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y3s},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			2.0*{0.0, 0.0, c[[3]]*y[[3]], \
				0.0, c[[5]]*y1s*y[[3]], c[[6]]*y2s*y[[3]], \
				c[[7]]*y1s*y2s*y[[3]]}
		],Evaluate@cpOPTs];


(* ::Section:: *)
(*The model responses are below*)


responsefunctions = Compile[{{c,_Real,1}, {F,_Real,1}, {M,_Real,1}, {ymax,_Real,1}},
	Module[{y,Q,Q1,Q2,Q3,exp},
		y = $gamma[F,M];
		Q = $Q[c[[2;;8]],y] - $Q[c[[2;;8]],ymax];
		Q1 = $Q1[c[[2;;8]],y];
		Q2 = $Q2[c[[2;;8]],y];
		Q3 = $Q3[c[[2;;8]],y];
		exp = Exp[Q];
		c[[1]]*{Total[Q1*exp], Total[Q2*exp], Total[Q3*exp]}
	], Evaluate@cpOPTs
];


strainenergy = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {F,_Real,1}, {M,_Real,1}, {ymax,_Real,1}},
		Module[{y,Q,exp},
			y = $gamma[F,M];
			Q = $Q[c[[2;;8]],y] - $Q[c[[2;;8]],ymax];
			exp = Exp[Q];
			c[[1]]*Total[(Exp[Q] - 1)]
		],cpOPTs
		]]


(* ::Subsection:: *)
(*The data processing functions are defined below*)


W1$fromdata = With[{cpOPTs = cpOPTs},
	Compile[{{T,_Real,1},{F,_Real,1},{M,_Real,1}},
		Module[{m},
			m = calc$m[F];
			m.dot2D[T,m]
		], cpOPTs
	]]

W2$fromdata = With[{cpOPTs = cpOPTs},
	Compile[{{T,_Real,1},{F,_Real,1},{M,_Real,1}},
		Module[{s},
			s = calc$s[F];
			s.dot2D[T,s]
		], cpOPTs
	]]
	
W3$fromdata = With[{cpOPTs = cpOPTs},
	Compile[{{T,_Real,1},{F,_Real,1},{M,_Real,1}},
		Module[{m,s},
			m = calc$m[F,M];
			s = calc$s[F,M];
			m.dot2D[T,s] * lambda$M[F,M] / lambda$S[F,M]
		], cpOPTs
	]]


End[]


EndPackage[]
