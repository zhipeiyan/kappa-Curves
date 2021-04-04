# &kappa;-Curves: Interpolation at Local Maximum Curvature
This is the original Wolfram Mathematica implementation of the Siggraph 2017 paper

http://people.tamu.edu/~yanzp/projects/kCurves/kCurves.pdf

## Usage
### Wolfram Mathematica
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

### C++
Download the [Eigen](https://eigen.tuxfamily.org/index.php) library.

Include the header file `kCurves.h`

The declaration of the two computing functions are
```
std::vector<std::vector<Eigen::Vector2d>> kCurveClosed(std::vector<Eigen::Vector2d>);
std::vector<std::vector<Eigen::Vector2d>> kCurveOpen(std::vector<Eigen::Vector2d>);
```

The input format `std::vector<Eigen::Vector2d>` is a list of your control points in 2D.

The output format is a list of triple points. Each triple represents the three control points of a quadratic Bezier curve segment according to your input points. The output vector has the same length as your input points for closed curves, and two shorter than your input points for open curves.

## Live Demo
Manipulate the curves in real time.
```
In[4]:= Manipulate[Graphics[BezierCurve[kCurveClosed[pts]], PlotRange -> {{-5, 5}, {-5, 5}}], {{pts, {{-1, 0}, {0, 1}, {1, 0}, {2, 1}}}, Locator}]
Out[4]:= 
```
![](http://people.tamu.edu/~yanzp/projects/kCurves/closed.png)

## Standalone version
Check https://github.com/zhipeiyan/kappa-Curves/releases

No package files needed.

## Code format
Plain Mathematica code format may be not reading friendly. Here's the screen shot of the code in the .nb format for some one doesnot have Wolfram Mathematica installed.

https://raw.githubusercontent.com/zhipeiyan/kappa-Curves/main/reading.svg
