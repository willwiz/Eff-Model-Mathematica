(* ::Package:: *)

BeginPackage["plotting`"]


pointPlot3D::usage::"ListPointPlot3D with changed default options"

listPlot3D::usage::"ListPlot3D with changed default options"

surfPlot3D::usage::"Plot3D with changed default options"

label::usage::"Given text with the styles of Black, no single letter italics, font 16 arial"


Begin["`Private`"]


(* ::Subsubsection:: *)
(*options*)


opts$Plt3D = {BaseStyle->Directive[14,Black], AxesStyle->Directive[Black, 16],
	PlotRange->All, BoxRatios->{1,1,1/GoldenRatio}, 
	FaceGrids->{Back,Right,Bottom}, Boxed->False,ViewPoint->10*{-1,-1.5,.8},
	Background->None, Prolog->{{EdgeForm[],Texture[{{{0,0,0,0}}}],
	Polygon[#,VertexTextureCoordinates->#]&[{{0,0},{1,0},{1,1}}]}}}

opts$PPlt3D=Join[{PlotStyle->{Red,PointSize[Medium]}}, opts$Plt3D];

opts$SPlt3D=Join[opts$Plt3D,{PlotStyle->None, ClippingStyle->None, MeshStyle->{Black,Thick}}];


opts$lbl = {Black,SingleLetterItalics->False,FontSize->16,FontFamily->"Arial"}


(* ::Subsubsection:: *)
(*Functions*)


pointPlot3D[dat_,opts:OptionsPattern[]]:=ListPointPlot3D[dat, opts, opts$PPlt3D]


listPlot3D[dat_,opts:OptionsPattern[]]:=ListPlot3D[dat,opts, opts$SPlt3D]


SetAttributes[surfPlot3D, HoldAll]
surfPlot3D[f_, {x_,xmin_,xmax_}, {y_, ymin_, ymax_}, opts:OptionsPattern[]]:=
	Plot3D[f, {x,xmin,xmax}, {y, ymin,ymax}, opts, 
		Evaluate[FilterRules[{opts$SPlt3D},Options[Plot3D]]]]


label[text_,opts:OptionsPattern[]]:=Thread[Style[text,opts,Sequence@@opts$lbl]]


End[]


EndPackage[]
