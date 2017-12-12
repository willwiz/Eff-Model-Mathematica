(* ::Package:: *)

BeginPackage["Effmodel`"]


toTensor2D::usage::"Takes a flattened vector to a square array";

cross2D::usage::"The 2D cross product used to find the orthogonal vector";

dot2D::usage::"The 2D dot product used to find the vector";

mult2D::usage::"matrix multiplication between 2 matrix in vector form";

trp2D::usage::"matrix transpose in vector form";

calc$m::usage::"Finds the new unit vector direction of M";

calc$s::usage::"Finds the new unit vector direction orthogonal of M";


lambda$M::usage::"compute the stretch along M";

lambda$S::usage::"compute the stretch along S";

phi::usage::"computes the shear angle";

gamma$1::usage::"compute the invariant";

gamma$2::usage::"compute the invariant";

gamma$3::usage::"compute the invariant";


$Qp;
$Q1p;
$Q2p;
$Q3p;


W1$fromdata::usage::"computes the first response function from some given data";

W2$fromdata::usage::"computes the second response function from some given data";

W3$fromdata::usage::"computes the third response function from some given data";


rf$gamma::usage::"the response function as a function of gamma";

responsefunctions::usage::"computes the response functions and puts them in an array";

se$gamma::usage::"The strain energy as a function of gamma";

strainenergy::usage::"computes the strain energy given the parameters, and strain";


obj$W$gamma::usage::"compiled code for the objective function as a function of gamma";

objW$c::usage::"the obj fun using the compile code";

obj$grad$W::usage::"the gradient of the obj fun using the compile code";

grad$W$c::usage::"the gradient of the obj fun using the compile code";

objW::usage::"the obj fun only partially using the compiled code";


eps$W$y::usage::"Calculate the error^2 at one data point";

rf$grad1::usage::"Calculate the gradient of the function at a point";
rf$grad2::usage::"Calculate the gradient of the function at a point";
rf$grad3::usage::"Calculate the gradient of the function at a point";

grad::usage::"Calculate the gradient of the function at a point";

grad$W$gamma::usage::"Calculates the gradient of the obj function obj$W$gamma at one point";

jac$W$gamma::usage::"Calculates the jacobian of the obj function obj$W$gamma";

cjac$W$gamma::usage::"Calculates the jacobian of the obj function obj$W$gamma using the analytical solution";


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


(* ::Subsubsection:: *)
(*No Shear component*)


$Qp = With[{cpOPTs = cpOPTs},
	Compile[{{y,_Real,1},{ys,_Real,1}},
		Module[{},
			{0.0, ys[[1]], ys[[2]], 2.0*y[[1]]*y[[2]], 
			ys[[1]]*ys[[1]], ys[[2]]*ys[[2]], ys[[1]]*y[[1]]*y[[2]], ys[[2]]*y[[1]]*y[[2]], 
			ys[[1]]*ys[[3]], ys[[2]]*ys[[3]], ys[[3]]*ys[[3]]}
		],cpOPTs]];
		
$Q1p = With[{cpOPTs = cpOPTs},
	Compile[{{y,_Real,1},{ys,_Real,1}},
		Module[{},
			{0.0, 2.0*y[[1]], 0.0, 2.0*y[[2]], 
			4.0*y[[1]]*ys[[1]], 0.0, 3.0*ys[[1]]*y[[2]], ys[[2]]*y[[2]], 
			2.0*y[[1]]*ys[[3]], 0.0, 0.0}
		],cpOPTs]];
		
$Q2p = With[{cpOPTs = cpOPTs},
	Compile[{{y,_Real,1},{ys,_Real,1}},
		Module[{},
			{0.0, 0.0, 2.0*y[[2]], 2.0*y[[1]], 
			0.0, 4.0*y[[2]]*ys[[2]], ys[[1]]*y[[1]], 3.0*ys[[2]]*y[[1]], 
			0.0, 2.0*y[[2]]*ys[[3]], 0.0}
		],cpOPTs]];
		
$Q3p = With[{cpOPTs = cpOPTs},
	Compile[{{y,_Real,1},{ys,_Real,1}},
		Module[{},
			{0.0, 0.0, 0.0, 0.0, 
			0.0, 0.0, 0.0, 0.0, 
			2.0*ys[[1]]*y[[3]], 2.0*ys[[2]]*y[[3]], 4.0*y[[3]]*ys[[3]]}
		],cpOPTs]];


(* ::Subsubsection:: *)
(*This defines the jacobians*)


(* ::Section:: *)
(*The model responses are below*)


rf$gamma = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {y,_Real,1}, {ymax,_Real,1}},
		Module[{Q,Q1,Q2,Q3,exp},
			With[{ys = y*y, ymaxs = ymax*ymax},
			Q = Total[c*($Qp[y,ys] - $Qp[ymax, ymaxs])];
			Q1 = Total[c*$Q1p[y,ys]];
			Q2 = Total[c*$Q2p[y,ys]];
			Q3 = Total[c*$Q3p[y,ys]];];
			exp = Exp[Q];
			c[[1]]*{Q1*exp, Q2*exp, Q3*exp}]
	,{{$Qp,_Real,1}}
	,cpOPTs]];


rf$grad1 = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {y,_Real,1}, {ymax,_Real,1}},
		Module[{ys = y*y,ymaxs = ymax*ymax,Qp,Q1p,Q1,exp,res(* = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}*)},
			Qp = ($Qp[y,ys] - $Qp[ymax,ymaxs]);
			Q1p = $Q1p[y,ys];
			exp = Exp[Total[c*Qp]];
			Q1 = Total[c*Q1p];
			(*Do[res = c[[1]]*(Q1p*exp + Q1*exp*Qp),{11}]*)
			res = c[[1]]*(Q1p*exp + Q1*exp*Qp);
			res[[1]] = Q1*exp;
			res
	],{{$Qp,_Real,1},{$Q1p,_Real,1}},cpOPTs]];


rf$grad2 = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {y,_Real,1}, {ymax,_Real,1}},
		Module[{ys = y*y,ymaxs = ymax*ymax,Qp,Q2p,Q2,exp,res(* = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}*)},
			Qp = ($Qp[y,ys] - $Qp[ymax,ymaxs]);
			Q2p = $Q2p[y,ys];
			exp = Exp[Total[c*Qp]];
			Q2 = Total[c*Q2p];
			(*Do[res = c[[1]]*(Q1p*exp + Q1*exp*Qp),{11}]*)
			res = c[[1]]*(Q2p*exp + Q2*exp*Qp);
			res[[1]] = Q2*exp;
			res
	],{{$Qp,_Real,1},{$Q2p,_Real,1}},cpOPTs]];


rf$grad3 = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {y,_Real,1}, {ymax,_Real,1}},
		Module[{ys = y*y,ymaxs = ymax*ymax,Qp,Q3p,Q3,exp,res(* = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}*)},
			Qp = ($Qp[y,ys] - $Qp[ymax,ymaxs]);
			Q3p = $Q3p[y,ys];
			exp = Exp[Total[c*Qp]];
			Q3 = Total[c*Q3p];
			(*Do[res = c[[1]]*(Q1p*exp + Q1*exp*Qp),{11}]*)
			res = c[[1]]*(Q3p*exp + Q3*exp*Qp);
			res[[1]] = Q3*exp;
			res
	],{{$Qp,_Real,1},{$Q3p,_Real,1}},cpOPTs]];


(* ::Subsection:: *)
(*Jacobian*)


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
			err = (W - Wdata);
			Total[Map[#.#&,Transpose[err]]]
		]
	,cpOPTs]]


objW$c[c_?(VectorQ[#,NumericQ]&),ymax_?(VectorQ[#,NumericQ]&),
	ydata_?(ArrayQ[#,2,NumericQ]&),Wdata_?(ArrayQ[#,2,NumericQ]&)]:=
		obj$W$gamma[c,ymax,ydata,Wdata]


grad$W$gamma = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {ymax,_Real,1},{ydata,_Real,1},{Wdata,_Real,1}},
		Module[{W,res = c, sigma},
			W = rf$gamma[c,ydata,ymax];
			sigma = (W - Wdata);
			res = 2.0*(sigma[[1]]*rf$grad1[c,ydata,ymax]+sigma[[2]]*rf$grad2[c,ydata,ymax]+
			sigma[[3]]*rf$grad3[c,ydata,ymax]);
			res
		]
	,cpOPTs]]


cjac$W$gamma = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {ymax,_Real,1},{ydata,_Real,2},{Wdata,_Real,2}},
		Table[grad$W$gamma[c,ymax,ydata[[i]],Wdata[[i]]],{i,Length[ydata]}],{{grad$W$gamma,_Real,1}}
	,cpOPTs]]


obj$grad$W = With[{cpOPTs = cpOPTs},
	Compile[{{c,_Real,1}, {ymax,_Real,1},{ydata,_Real,2},{Wdata,_Real,2}},
		Total[Table[grad$W$gamma[c,ymax,ydata[[i]],Wdata[[i]]],{i,Length[ydata]}]],{{grad$W$gamma,_Real,1}}
	,cpOPTs]]


grad$W$c[c_?(VectorQ[#,NumericQ]&),ymax_?(VectorQ[#,NumericQ]&),
	ydata_?(ArrayQ[#,2,NumericQ]&),Wdata_?(ArrayQ[#,2,NumericQ]&)]:=
		obj$grad$W[c,ymax,ydata,Wdata]


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
