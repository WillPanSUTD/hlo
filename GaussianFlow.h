// This class implements the Gaussian Curvature Flow from the following paper:
//    Huanxi Zhao, Guoliang Xu: "Triangular surface mesh fairing via Gaussian curvature flow"
//    Journal of Computational and Applied Mathematics, Volume 195, Issues 1?, 15 October 2006, Pages 300-311
///https://ac.els-cdn.com/S0377042705004942/1-s2.0-S0377042705004942-main.pdf?_tid=16cc3b3c-428c-4954-84b7-5569b21160af&acdnat=1536217275_8650bf2e1a2e67511746f9f995a84ddc

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/gaussian_curvature.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

#include <vector>
#include <set>

using namespace std;

struct GaussianFlow
{
	const char* name = "GaussianFlow";
	double timestep = 0.01;
	const double EPSILON = 0.001;
	const double ALPHA = 0.0005;
	const double BETA = 2.0;

	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd & V, Eigen::MatrixXi& F, vector<int> & domain,
		Eigen::MatrixXd & VNormal, vector<vector<int> > & NeighborList, Eigen::MatrixXd & derivative)
	{		
		igl::gaussian_curvature(V, F, curvatureGaussian);

		Eigen::SparseMatrix<double> cotangentMatrix, mass, invMass;
		igl::cotmatrix(V, F, cotangentMatrix);

		igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, mass);
		igl::invert_diag(mass, invMass);
		meanCurvatureNormal = invMass * (cotangentMatrix*V);

		for (int i = 0; i < domain.size(); ++i)
		{
			int vIdx = domain[i];

			double gaussian = curvatureGaussian.coeff(vIdx, 0) * invMass.coeff(vIdx, vIdx);
			Eigen::Vector3d normal = VNormal.row(vIdx);
			Eigen::Vector3d mcNormal = meanCurvatureNormal.row(vIdx); // = -H n
			double meanCurvature = mcNormal.norm();
			if (normal.dot(mcNormal) < 0)
				meanCurvature *= -1.0f;
			
			double chi = 1.0 / (1 + pow(gaussian, BETA));
			chi *= timestep;
			if (abs(meanCurvature) < EPSILON)
				derivative.row(vIdx) = chi * meanCurvature * normal;
			else if (gaussian > 0)
				derivative.row(vIdx) = chi * (meanCurvature > 0 ? 1.0 : -1.0) * gaussian * normal;
			else
				derivative.row(vIdx) = ALPHA * chi * gaussian * normal;
			//derivative *= -1;
		}
	}

private:
	Eigen::MatrixXd curvatureGaussian, meanCurvatureNormal;
};