/*
Author: Zhipei Yan
Date: 04/03/2021
Description: This is a simplified C++ implementation of the kappa-Curves paper.
http://people.tamu.edu/~yanzp/projects/kCurves/kCurves.pdf
The code is mainly converted from the Mathematica version.
For practical engineering projects, further details e.g. numerical error checking,
colinear cases may be required. The complete release version of the C++ code is
still in progress.
Thanks.
*/

#pragma once

#include <vector>
#include <Eigen/Dense>

namespace KCURVE
{
	double maxParam(Eigen::Vector2d c0, Eigen::Vector2d cp, Eigen::Vector2d c2)
	{
		double a = (c2 - c0).dot(c2 - c0);
		double b = 3 * (c2 - c0).dot(c0 - cp);
		double c = (3 * c0 - 2 * cp - c2).dot(c0 - cp);
		double d = -(c0 - cp).dot(c0 - cp);
		// solve a t^3 + b t^2 + c t + d = 0 in [0, 1]
		// TODO: if c0==c2, solve a quadratic or linear equation.
		// https://en.wikipedia.org/wiki/Cubic_equation
		// u = t + b / 3a
		double p = (3 * a * c - b * b) / 3 / a / a;
		double q = (2 * b * b * b - 9 * a * b * c + 27 * a * a * d) / 27 / a / a / a;
		// solve u^3 + p u + q = 0
		if (4 * p * p * p + 27 * q * q >= 0)
		{
			// single real root
			return cbrt(-q / 2 + sqrt(q * q / 4 + p * p * p / 27)) + cbrt(-q / 2 - sqrt(q * q / 4 + p * p * p / 27)) - b / 3 / a;
		}
		else
		{
			// three real roots
			for (int k = 0; k < 3; ++k)
			{
				double t = 2 * sqrt(-p / 3) * cos(1. / 3 * acos(3 * q / 2. / p * sqrt(-3 / p)) - 2 * EIGEN_PI * k / 3.) - b / 3 / a;
				if (0 <= t && t <= 1)
					return t;
			}
		}
		// error
		return -1;
	}

	double area(Eigen::Vector2d A, Eigen::Vector2d B, Eigen::Vector2d C)
	{
		Eigen::Matrix2d mat;
		mat << B - A, C - B;
		return fabs(mat.determinant() / 2.);
	}

	double lambda(Eigen::Vector2d c0, Eigen::Vector2d c1, Eigen::Vector2d c3, Eigen::Vector2d c4)
	{
		return sqrt(area(c0, c1, c3)) / (sqrt(area(c0, c1, c3)) + sqrt(area(c1, c3, c4)));
	}

	// when i = -1, 0, 1, ..., n-1, n, n+1
	inline int mod(int i, int n)
	{
		return ((i + n) % n);
	}

	std::vector<std::vector<Eigen::Vector2d>> kCurveClosed(std::vector<Eigen::Vector2d> pts)
	{
		int n = pts.size();

		// return bezier control points, a vector of triple of 2d points, {{q0,q1,q2}, {q0,q1,q2}, ...}
		std::vector<std::vector<Eigen::Vector2d>> ctrls(n, std::vector<Eigen::Vector2d>(3, Eigen::Vector2d::Zero()));
		for (int i = 0; i < n; ++i)
		{
			ctrls[i][0] = (pts[i] + pts[mod(i - 1, n)]) / 2;
			ctrls[i][1] = pts[i];
			ctrls[i][2] = (pts[i] + pts[mod(i + 1, n)]) / 2;
		}

		std::vector<double> t(n), ld(n);
		Eigen::MatrixXd matA = Eigen::MatrixXd::Zero(n, n);
		Eigen::MatrixX2d matB(n, 2);
		for (int i = 0; i < n; ++i)
			matB.row(i) = pts[i];

		for (int iter = 0; iter < 50; ++iter)
		{
			// 50 iterations by default. Chanege to different iterations on demand
			for (int i = 0; i < n; ++i)
			{
				t[i] = maxParam(ctrls[i][0], pts[i], ctrls[i][2]);
				ld[i] = lambda(ctrls[i][0], ctrls[i][1], ctrls[mod(i + 1, n)][1], ctrls[mod(i + 1, n)][2]);
			}
			for (int i = 0; i < n; ++i)
			{
				matA(i, mod(i - 1, n)) = (1 - ld[mod(i - 1, n)]) * (1 - t[i]) * (1 - t[i]);
				matA(i, i) = ld[mod(i - 1, n)] * (1 - t[i]) * (1 - t[i]) + 2 * (1 - t[i]) * t[i] + (1 - ld[i]) * t[i] * t[i];
				matA(i, mod(i + 1, n)) = ld[i] * t[i] * t[i];
			}

			// linear solver of Ax=B
			Eigen::MatrixX2d corners = matA.colPivHouseholderQr().solve(matB);

			for (int i = 0; i < n; ++i)
				ctrls[i][1] = corners.row(i);
			for (int i = 0; i < n; ++i)
				ctrls[mod(i + 1, n)][0] = ctrls[i][2] = (1 - ld[i]) * ctrls[i][1] + ld[i] * ctrls[mod(i + 1, n)][1];
		}
		return ctrls;
	}

	std::vector<std::vector<Eigen::Vector2d>> kCurveOpen(std::vector<Eigen::Vector2d> pts)
	{
		int n = pts.size() - 2;

		// return bezier control points, a vector of triple of 2d points, {{q0,q1,q2}, {q0,q1,q2}, ...}
		std::vector<std::vector<Eigen::Vector2d>> ctrls(n, std::vector<Eigen::Vector2d>(3, Eigen::Vector2d::Zero()));
		for (int i = 1; i <= n; ++i)
		{
			ctrls[i - 1][0] = (pts[i] + pts[i - 1]) / 2;
			ctrls[i - 1][1] = pts[i];
			ctrls[i - 1][2] = (pts[i] + pts[i + 1]) / 2;
		}
		ctrls[0][0] = pts[0];
		ctrls[n - 1][2] = pts[n + 1];

		std::vector<double> t(n), ld(n - 1);
		Eigen::MatrixXd matA = Eigen::MatrixXd::Zero(n, n);
		Eigen::MatrixX2d matB(n, 2);
		for (int i = 0; i < n; ++i)
			matB.row(i) = pts[i + 1];

		for (int iter = 0; iter < 50; ++iter)
		{
			// 50 iterations by default. Chanege to different iterations on demand
			for (int i = 0; i < n; ++i)
				t[i] = maxParam(ctrls[i][0], pts[i + 1], ctrls[i][2]);
			for (int i = 0; i < n - 1; ++i)
				ld[i] = lambda(ctrls[i][0], ctrls[i][1], ctrls[i + 1][1], ctrls[i + 1][2]);

			matA(0, 0) = 2 * (1 - t[0]) * t[0] + (1 - ld[0]) * t[0] * t[0];
			matA(0, 1) = ld[0] * t[0] * t[0];
			for (int i = 1; i < n - 1; ++i)
			{
				matA(i, i - 1) = (1 - ld[i - 1]) * (1 - t[i]) * (1 - t[i]);
				matA(i, i) = ld[i - 1] * (1 - t[i]) * (1 - t[i]) + 2 * (1 - t[i]) * t[i] + (1 - ld[i]) * t[i] * t[i];
				matA(i, i + 1) = ld[i] * t[i] * t[i];
			}
			matA(n - 1, n - 2) = (1 - ld[n - 2]) * (1 - t[n - 1]) * (1 - t[n - 1]);
			matA(n - 1, n - 1) = ld[n - 2] * (1 - t[n - 1]) * (1 - t[n - 1]) + 2 * (1 - t[n - 1]) * t[n - 1];
			matB.row(0) = pts[1] - (1 - t[0]) * (1 - t[0]) * pts[0];
			matB.row(n - 1) = pts[n] - t[n - 1] * t[n - 1] * pts[n + 1];

			// linear solver of Ax=B
			Eigen::MatrixX2d corners = matA.colPivHouseholderQr().solve(matB);

			for (int i = 0; i < n; ++i)
				ctrls[i][1] = corners.row(i);
			for (int i = 0; i < n - 1; ++i)
				ctrls[i + 1][0] = ctrls[i][2] = (1 - ld[i]) * ctrls[i][1] + ld[i] * ctrls[i + 1][1];
		}
		return ctrls;
	}
}

// Usage:
/* int main()
{
	std::vector<Eigen::Vector2d> input(4);
	input[0] << 0, 0;
	input[1] << 0, 1;
	input[2] << 2, 2;
	input[3] << 1, 0;
	std::vector<std::vector<Eigen::Vector2d>> bezierClosed = KCURVE::kCurveClosed(input);
	std::vector<std::vector<Eigen::Vector2d>> bezierOpen = KCURVE::kCurveOpen(input);
	return 0;
} */