(* ::Package:: *)

BeginPackage["Tools`"]


SSEtotal::usage::"Computes the SSE of an given array"


SSEt$vec::usage::"Computes the SSE of an given array"


epsTotal::usage::"Computes the difference of each data point from the main"


Begin["`Private`"]


rtOPTs = {"CatchMachineOverflow"->False, "CatchMachineUnderflow"->False,
	"CatchMachineIntegerOverflow"->False, "CompareWithTolerance"->True,
	"EvaluateSymbolically"->False,"RuntimeErrorHandler"->Evaluate,"WarningMessages"->True};


cpOPTs = {CompilationOptions->{"ExpressionOptimization"->True, "InlineExternalDefinitions"->True,"InlineCompiledFunctions"->True},
	RuntimeAttributes->{Listable}, 
	Parallelization->True,
	RuntimeOptions->rtOPTs};


epsTotal = With[{cpOPTs = cpOPTs},
	Compile[{{x,_Real,1}},
		With[{xmean = Mean[x]},Map[#-xmean&,x]]
		,cpOPTs]];


SSEtotal = With[{cpOPTs = cpOPTs},
	Compile[{{x,_Real,1}},
		Module[{eps},
			eps = epsTotal[x];
			eps.eps]
		,cpOPTs]]


SSEt$vec = With[{cpOPTs = cpOPTs},
	Compile[{{x,_Real,2}},
		Module[{eps},
			eps = Map[epsTotal,Transpose[x]];
			Total[Map[#.#&,eps]]]
		,cpOPTs]]


End[]


EndPackage[]
