(* ::Package:: *)

BeginPackage["KappaCurves`"]

kCurveClosed::usage="kCurveClosed[{point1, point2, ..., pointN}] gives the Bezier control point tuples. N mustt be at least three.";
kCurveOpen::usage="kCurveOpen[{point1, point2, ..., pointN}] gives the Bezier control point tuples. N mustt be at least three.";

Begin["Private`"]
(*parameter t at point p when p is local maximum curvature point*)
maxParam[{c0_,p_,c2_}]:=NSolve[(c2-c0).(c2-c0) t^3+3(c2-c0).(c0-p) t^2+(3c0-2p-c2).(c0-p)t-(c0-p).(c0-p)==0&&0<=t<=1,t,Reals][[1,1,2]]

area[{a_,b_,c_}]:=Abs[Det[{b-a,c-b}]/2]

(*the interpolation coefficient of the joint points*)
lambda[{c0_,c1_,c3_,c4_}]:=Sqrt[area[{c0,c1,c3}]]/(Sqrt[area[{c0,c1,c3}]]+Sqrt[area[{c1,c3,c4}]])//N


(*closed \[Kappa]Curves from given input points. pts must be at least three points*)
kCurveClosed[pts_, iterations_:50]:=Module[{n=Length[pts],ctrls,t,\[Lambda],mat},
	(*control points of the bezier control polygon*)
	ctrls=Table[{(pts[[i]]+pts[[Mod[i-1,n,1]]])/2,pts[[i]],(pts[[i]]+pts[[Mod[i+1,n,1]]])/2},{i,n}];
	
	Do[
		t=Table[maxParam[{ctrls[[i,1]],pts[[i]],ctrls[[i,3]]}],{i,n}];
		\[Lambda]=Table[lambda[Join[ctrls[[i,{1,2}]],ctrls[[Mod[i+1,n,1],{2,3}]]]],{i,n}];
		
		mat=Table[0,n,n];
		Do[mat[[i,{Mod[i-1,n,1],i,Mod[i+1,n,1]}]]={(1-\[Lambda][[Mod[i-1,n,1]]])(1-t[[i]])^2,\[Lambda][[Mod[i-1,n,1]]](1-t[[i]])^2+2(1-t[[i]])t[[i]]+(1-\[Lambda][[i]])t[[i]]^2,\[Lambda][[i]]t[[i]]^2},{i,n}];
		
		(*off-curve Bezier control points*)
		ctrls[[All,2]]=Inverse[mat].pts;
		(*joint points between quadratic segments*)
		Do[ctrls[[Mod[i+1,n,1],1]]=ctrls[[i,3]]=(1-\[Lambda][[i]])ctrls[[i,2]]+\[Lambda][[i]]ctrls[[Mod[i+1,n,1],2]],{i,n}],
		
		iterations(*50 by default*)
	];
	ctrls
]


(* open \[Kappa]Curve from given input points. pts must be at least three points*)
kCurveOpen[pts_, iterations_:50]:=Module[{n=Length[pts]-2,ctrls,t,\[Lambda],matA,matB},
	(*control points of the bezier control polygon*)
	ctrls=Table[{(pts[[i]]+pts[[i-1]])/2,pts[[i]],(pts[[i]]+pts[[i+1]])/2},{i,2,Length[pts]-1}];
	(*the two end points will be used as constants*)
	ctrls[[1,1]]=pts[[1]];
	ctrls[[-1,-1]]=pts[[-1]];

	Do[
		t=Table[maxParam[{ctrls[[i,1]],pts[[i+1]],ctrls[[i,3]]}],{i,n}];
		\[Lambda]=Table[lambda[Join[ctrls[[i,{1,2}]],ctrls[[i+1,{2,3}]]]],{i,n-1}];
		
		(*tridiagonal matrix A*)
		matA=DiagonalMatrix[Table[2(1-t[[i]])t[[i]],{i,n}]];
		Do[matA[[i,{i,i+1}]]+={1-\[Lambda][[i]],\[Lambda][[i]]}t[[i]]^2,{i,n-1}];
		Do[matA[[i,{i-1,i}]]+={1-\[Lambda][[i-1]],\[Lambda][[i-1]]}(1-t[[i]])^2,{i,2,n}];
		(*two end points are constants*)
		matB=pts[[2;;-2]];
		matB[[1]]-=(1-t[[1]])^2 pts[[1]];
		matB[[-1]]-=t[[-1]]^2 pts[[-1]];
		
		(*off-curve Bezier control points*)
		ctrls[[All,2]]=Inverse[matA].matB;
		(*joint points between quadratic segments*)
		Do[ctrls[[i+1,1]]=ctrls[[i,3]]=(1-\[Lambda][[i]])ctrls[[i,2]]+\[Lambda][[i]]ctrls[[i+1,2]],{i,n-1}],
		
		iterations(*50 by default*)
	];
	ctrls
]

End[]
EndPackage[]
