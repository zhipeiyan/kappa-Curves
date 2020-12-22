# &kappa;-Curves: Interpolation at Local Maximum Curvature
This is the original Wolfram Mathematica implementation of the Siggraph 2017 paper

http://people.tamu.edu/~yanzp/projects/kCurves/kCurves.pdf

## Usage
Put the .wl package in the same path of your notebook file and load the package
```
In[1]:= Get[NotebookDirectory[] <> "kCurves.wl"]
```
Compute an open curve interpolating the input points
```
In[2]:= kCurveOpen[{{-1,0},{0,1},{1,0},{2,1}}]
Out[2]:= {{{-1, 0}, {-0.08187, 1.71183}, {0.5, 0.5}}, {{0.5, 0.5}, {1.08187, -0.711831}, {2, 1}}}
```
Draw the curve
```
In[3]:= Graphics[BezierCurve[%]]
Out[3]:= 
```
![](http://people.tamu.edu/~yanzp/projects/kCurves/four.png)

## Live Demo
Manipulate the curves in real time.
```
In[4]:= Manipulate[Graphics[BezierCurve[kCurveClosed[pts]], PlotRange -> {{-5, 5}, {-5, 5}}], {{pts, {{-1, 0}, {0, 1}, {1, 0}, {2, 1}}}, Locator}]
Out[4]:= 
```
![](http://people.tamu.edu/~yanzp/projects/kCurves/closed.png)

## Code format
Plain Mathematica code format may be not reading friendly. Here's the screen shot of the code in the .nb format for some one doesnot have Wolfram Mathematica installed.

https://raw.githubusercontent.com/zhipeiyan/kappa-Curves/main/reading.svg
