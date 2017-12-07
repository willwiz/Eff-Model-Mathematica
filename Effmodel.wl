(* ::Package:: *)

BeginPackage["Effmodel`"]


toTensor2D::usage::"Takes a flattened vector to a square array"

cross2D::usage::"The 2D cross product used to find the orthogonal vector"

dot2D::usage::"The 2D dot product used to find the vector"

mult2D::usage::"matrix multiplication between 2 matrix in vector form"

trp2D::usage::"matrix transpose in vector form"

calc$m::usage::"Finds the new unit vector direction of M"

calc$s::usage::"Finds the new unit vector direction orthogonal of M"


lambda$M::usage::"compute the stretch along M"

lambda$S::usage::"compute the stretch along S"

phi::usage::"computes the shear angle"

gamma$1::usage::"compute the invariant"

gamma$2::usage::"compute the invariant"

gamma$3::usage::"compute the invariant"


W1$fromdata::usage::"computes the first response function from some given data"

W2$fromdata::usage::"computes the second response function from some given data"

W3$fromdata::usage::"computes the third response function from some given data"


rf$gamma::usage::"the response function as a function of gamma"

responsefunctions::usage::"computes the response functions and puts them in an array"

se$gamma::usage::"The strain energy as a function of gamma"

strainenergy::usage::"computes the strain energy given the parameters, and strain"


obj$W$gamma::usage::"compiled code for the objective function as a function of gamma"

objW$c::usage::"the obj fun using the compile code"

objW::usage::"the obj fun only partially using the compiled code"


eps$W$y::usage::"Calculate the error^2 at one data point"

grad::usage::"Calculate the gradient of the function at a point"

jac$W$gamma::usage::"Calculates the jacobian of the obj function obj$W$gamma"


(* ::Section:: *)
(*Begin Function Definitions*)


Begin["`Private`"];


rtOPTs = {"CatchMachineOverflow"->False, "CatchMachineUnderflow"->False,
	"CatchMachineIntegerOverflow"->False, "CompareWithTolerance"->True,
	"EvaluateSymbolically"->False,"RuntimeErrorHandler"->Evaluate, 
	"WarningMessages"->True};


cpOPTs = {CompilationOptions->{"ExpressionOptimization"->True, 
	"InlineExternalDefinitions"->True,"InlineCompiledFunctions"->True},
	RuntimeAttributes->{Listable}, 
	Parallelization->True,
	RuntimeOptions->rtOPTs};


(* ::Section:: *)
(*These are the Utility functions*)


toTensor2D = Compile[{{f,_Real,1}},{{f[[1]],f[[2]]},{f[[3]],f[[4]]}},Evaluate@cpOPTs];

dot2D = Compile[{{F,_Real,1},{M,_Real,1}},
					{F[[1]]*M[[1]]+F[[2]]*M[[2]], \
					F[[3]]*M[[1]]+F[[4]]*M[[2]]},Evaluate@cpOPTs];

mult2D = With[{cpOPTs = cpOPTs},
	Compile[{{F,_Real,1},{G,_Real,1}},
		{F[[1]]*G[[1]]+F[[2]]*G[[3]], F[[1]]*G[[2]]+F[[2]]*G[[4]],
			F[[3]]*G[[1]]+F[[4]]*G[[3]],F[[3]]*G[[2]]+F[[4]]*G[[4]]}
		,cpOPTs]];
					
trp2D = With[{cpOPTs = cpOPTs},
	Compile[{{F, _Real,1}}, 
		F[[{1,3,2,4}]]
	,cpOPTs]];

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


(* ::Subsection:: *)
(*This defines Q and its derivatives*)


(* ::Subsubsection::Closed:: *)
(*Traditional form*)


(*$Q = Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y3s,y3c},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			y3c = y3s*y[[3]];
			y3q = y3s*y3s;
			{c[[1]]*y1s, c[[2]]*y2s, c[[3]]*y3s,
				c[[4]]*y[[1]]*y[[2]], c[[5]]*y[[1]]*y3s,
				c[[6]]*y[[2]]*y3s, c[[7]]*y[[1]]*y[[2]]*y3s
				}
		],Evaluate@cpOPTs];
		
$Q1 = Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y3s,y3c},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			y3c = y3s*y[[3]];
			{2.0*c[[1]]*y[[1]], 0.0, 0.0,
				c[[4]]*y[[2]], c[[5]]*y3s,
				0.0, c[[7]]*y[[2]]*y3s
				}
		],Evaluate@cpOPTs];
		
$Q2 = Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y3s,y3c},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			{0.0, 2.0*c[[2]]*y[[2]], 0.0,
				c[[4]]*y[[1]], 0.0,
				c[[6]]*y3s, c[[7]]*y[[1]]*y3s
				}
		],Evaluate@cpOPTs];
		
$Q3 = Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y3s},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			2.0*{0.0, 0.0, c[[3]]*y[[3]],
				0.0, c[[5]]*y[[1]]*y[[3]],
				c[[6]]*y[[2]]*y[[3]], c[[7]]*y[[1]]*y[[2]]*y[[3]]
				}
		],Evaluate@cpOPTs];*)


(* ::Subsubsection:: *)
(*No Shear component*)


$Q = Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y1q,y2q,y3s,y3q},
			y1s = y[[1]]*y[[1]];
			y1q = y1s*y1s;
			y2s = y[[2]]*y[[2]];
			y2q = y2s*y2s;
			y3s = y[[3]]*y[[3]];
			y3q = y3s*y3s;
			{c[[1]]*y1s, c[[2]]*y2s, 2.0*c[[3]]*y[[1]]*y[[2]], c[[4]]*y1q, 
			c[[5]]*y2q, c[[6]]*y1s*y[[1]]*y[[2]], c[[7]]*y2s*y[[1]]*y[[2]], 
			c[[8]]*y1s*y3s, c[[9]]*y2s*y3s, c[[10]]*y3q}
		],Evaluate@cpOPTs];
		
$Q1 = Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y1q,y2q,y3s,y3q},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			y3q = y3s*y3s;
			{2.0*c[[1]]*y[[1]], 0.0, 2.0*c[[3]]*y[[2]], 4.0*c[[4]]*y1s*y[[1]], 
			0.0, 3.0*c[[7]]*y1s*y[[2]], c[[8]]*y2s*y[[2]], 
			2.0*c[[9]]*y[[1]]*y3s, 0.0, 0.0}
		],Evaluate@cpOPTs];
		
$Q2 = Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y1q,y2q,y3s,y3q},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			y3q = y3s*y3s;
			{0.0, 2.0*c[[2]]*y[[2]], 2.0*c[[3]]*y[[1]], 0.0, 
			4.0*c[[5]]*y2s*y[[2]], c[[6]]*y1s*y[[1]], 3.0*c[[7]]*y2s*y[[1]], 
			0.0, 2.0*c[[9]]*y[[2]]*y3s, 0.0}
		],Evaluate@cpOPTs];
		
$Q3 = Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y1q,y2q,y3s,y3q},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			y3q = y3s*y3s;
			{0.0, 0.0, 0.0, 0.0, 
			0.0, 0.0, 0.0, 
			2.0*c[[8]]*y1s*y[[3]], 2.0*c[[9]]*y2s*y[[3]], 4.0*c[[10]]*y3s*y[[3]]}
		],Evaluate@cpOPTs];


(* ::Subsubsection:: *)
(*This defines the jacobians*)


$Q$grad = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y1q,y2q,y3s,y3q},
			y1s = y[[1]]*y[[1]];
			y1q = y1s*y1s;
			y2s = y[[2]]*y[[2]];
			y2q = y2s*y2s;
			y3s = y[[3]]*y[[3]];
			y3q = y3s*y3s;
			{0.0, y1s, y2s, 2.0*y[[1]]*y[[2]], y1q, 
			y2q, y1s*y[[1]]*y[[2]], y2s*y[[1]]*y[[2]], 
			y1s*y3s, y2s*y3s, y3q}
		],cpOPTs]]
		
$Q1$grad = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y1q,y2q,y3s,y3q},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			y3q = y3s*y3s;
			{0.0, 2.0*y[[1]], 0.0, 2.0*y[[2]], 4.0*y1s*y[[1]], 
			0.0, 3.0*y1s*y[[2]], y2s*y[[2]], 
			2.0*y[[1]]*y3s, 0.0, 0.0}
		],cpOPTs]];
		
$Q2$grad = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y1q,y2q,y3s,y3q},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			y3q = y3s*y3s;
			{0.0, 0.0, 2.0*y[[2]], 2.0*y[[1]], 0.0, 
			4.0*y2s*y[[2]], y1s*y[[1]], 3.0*y2s*y[[1]], 
			0.0, 2.0*y[[2]]*y3s, 0.0}
		],cpOPTs]];
		
$Q3$grad = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1},{y,_Real,1}},
		Module[{y1s,y2s,y1q,y2q,y3s,y3q},
			y1s = y[[1]]*y[[1]];
			y2s = y[[2]]*y[[2]];
			y3s = y[[3]]*y[[3]];
			y3q = y3s*y3s;
			{0.0, 0.0, 0.0, 0.0, 0.0, 
			0.0, 0.0, 0.0, 
			2.0*y1s*y[[3]], 2.0*y2s*y[[3]], 4.0*y3s*y[[3]]}
		],cpOPTs]];


(* ::Section:: *)
(*The model responses are below*)


rf$gamma = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {y,_Real,1}, {ymax,_Real,1}},
		Module[{Q,Q1,Q2,Q3,exp, $c=c[[2;;11]]},
			Q = $Q[$c,y] - $Q[$c,ymax];
			Q1 = $Q1[$c,y];
			Q2 = $Q2[$c,y];
			Q3 = $Q3[$c,y];
			exp = Exp[Total[Q]];
			c[[1]]*{Total[Q1]*exp, Total[Q2]*exp, Total[Q3]*exp}
	],cpOPTs]]
	

(*rf$gamma = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {y,_Real,1}, {ymax,_Real,1}},
		Module[{Q,Q1,Q2,Q3,exp},
			Q = $Q[c[[2;;8]],y] - $Q[c[[2;;8]],ymax];
			Q1 = $Q1[c[[2;;8]],y];
			Q2 = $Q2[c[[2;;8]],y];
			Q3 = $Q3[c[[2;;8]],y];
			exp = Exp[Q];
			c[[1]]*{Total[Q1*exp], Total[Q2*exp], Total[Q3*exp]}
	],cpOPTs]]*)
	
	
(*rf$gamma = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {y,_Real,1}, {ymax,_Real,1}},
		Module[{Q,Q1,Q2,Q3,exp},
			Q = $Q[c[[9;;16]],y] - $Q[c[[9;;16]],ymax];
			Q1 = $Q1[c[[9;;16]],y];
			Q2 = $Q2[c[[9;;16]],y];
			Q3 = $Q3[c[[9;;16]],y];
			exp = Exp[Q];
			{Total[c[[1;;8]]*Q1*exp], Total[c[[1;;8]]*Q2*exp], Total[c[[1;;8]]*Q3*exp]}
	],cpOPTs]]*)


(* ::Subsection:: *)
(*Jacobian*)


jac$gamma = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {y,_Real,1}, {ymax,_Real,1}},
		Module[{Q,Q1,Q2,Q3,exp, $c=c[[2;;11]]},
			Q = $Q[$c,y] - $Q[$c,ymax];
			Q1 = $Q1[$c,y];
			Q2 = $Q2[$c,y];
			Q3 = $Q3[$c,y];
			exp = Exp[Total[Q]];
			c[[1]]*{Total[Q1]*exp, Total[Q2]*exp, Total[Q3]*exp}
	],cpOPTs]]


(*rf$grad = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {y,_Real,1}, {ymax,_Real,1}},
		Q = $Q[c[[2;;8]],y] - $Q[c[[2;;8]],ymax];
		Qgrad = $Q$grad[c,y] - $Q$grad[c,y];
		
	,cpOPTs]]*)


responsefunctions = Compile[{{c,_Real,1}, {F,_Real,1}, {M,_Real,1}, {ymax,_Real,1}},
	Module[{y},
		y = $gamma[F,M];
		rf$gamma[c,y,ymax]
	], Evaluate@cpOPTs
]


(*se$gamma = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {y,_Real,1}, {ymax,_Real,1}},
		Module[{Q,scaling},
			Q = $Q[c[[2;;8]],y];
			scaling = Exp[-$Q[c[[2;;8]],ymax]];
			c[[1]]*Total[scaling*(Exp[Q] - 1)]
		]
	,cpOPTs]]*)


(*strainenergy = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {F,_Real,1}, {M,_Real,1}, {ymax,_Real,1}},
		Module[{y},
			y = $gamma[F,M];
			se$gamma[c,y,ymax]
		],cpOPTs
		]]*)


(* ::Section:: *)
(*The Obj function are below*)


obj$W$gamma = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {ymax,_Real,1},{ydata,_Real,2},{Wdata,_Real,2}},
		Module[{W,err,sse,size},
			W = Map[rf$gamma[c,#,ymax]&,ydata];
			(*W = rf$gamma[c,ydata,ymax];*)
			err = (W - Wdata);
			Total[Map[#.#&,Transpose[err]]]
		]
	,cpOPTs]]


objW$c[c_?(VectorQ[#,NumericQ]&),ymax_?(VectorQ[#,NumericQ]&),
	ydata_?(ArrayQ[#,2,NumericQ]&),Wdata_?(ArrayQ[#,2,NumericQ]&)]:=
		obj$W$gamma[c,ymax,ydata,Wdata]


objW[c_?(VectorQ[#,NumericQ]&),ymax_?(VectorQ[#,NumericQ]&),
	ydata_?(ArrayQ[#,2,NumericQ]&),Wdata_?(ArrayQ[#,2,NumericQ]&)]:=
		Block[{W,err},
			W = rf$gamma[c,ydata,ymax];
			err = (W - Wdata);
			Total[err *err,2]
		]


(* ::Subsection:: *)
(*Calculating the covariance*)


Needs["NumericalCalculus`"]


c$eps$W$y = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {ymax,_Real,1},{ydata,_Real,1},{Wdata,_Real,1}},
		Module[{W,err,sse,size},
			err = (rf$gamma[c,ydata,ymax] - Wdata);
			err.err
		]
	,cpOPTs]]
	
eps$W$y[c_?(VectorQ[#,NumericQ]&),ymax_?(VectorQ[#,NumericQ]&),
	ydata_?(VectorQ[#,NumericQ]&),Wdata_?(VectorQ[#,NumericQ]&)]:=
		c$eps$W$y[c,ymax,ydata,Wdata]


SetAttributes[hold, HoldAll]

hold[f_,r_]:=f/.r

grad[optres_, ymax_, ydata_, Wdata_]:= 
	Module[{n=Length[optres],x,par}, 
		par=Table[x[j],{j,n}];
		Table[ND[hold[eps$W$y[par, ymax, ydata, Wdata],Drop[Thread[par->optres],{j}]],
			par[[j]], optres[[j]]],{j,Length[par]}]]


jac$W$gamma[optres_, ymax_, ydata_, Wdata_]:=
	Module[{m},
		m = Length[ydata];
		Table[grad[optres, ymax, ydata[[i]], Wdata[[i]]],{i,m}]
	]


(* ::Section:: *)
(*The data processing functions are defined below*)


W1$fromdata = With[{cpOPTs = cpOPTs},
	Compile[{{T,_Real,1},{F,_Real,1},{M,_Real,1}},
		Module[{m},
			m = calc$m[F,M];
			m.dot2D[T,m]
		], cpOPTs
	]]

W2$fromdata = With[{cpOPTs = cpOPTs},
	Compile[{{T,_Real,1},{F,_Real,1},{M,_Real,1}},
		Module[{s},
			s = calc$s[F,M];
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
