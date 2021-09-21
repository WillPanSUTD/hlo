#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/gaussian_curvature.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>

// Provides various methods to analyze the quality of a triangle mesh
class MeshAnalysis
{
public:

	// Initialize the analyzer. The references to V and F will be kept during the life time of the analyzer
	MeshAnalysis(Eigen::MatrixXd & V, Eigen::MatrixXi & F)
		: V(V), F(F), curvatureGaussian(V.rows(), 1), curvatureAbsoluteMean(V.rows(), 1), meanCurvatureNormal(V.rows(), 3)
	{ }

	// Perform the analyzation based on the current state of V and F
	void analyze()
	{
		//Calculate Gaussian curvature
		// !! This is not the actual Gaussian curvature. It is the integral of Gaussian curvature over the mesh surface
		igl::gaussian_curvature(V, F, curvatureGaussian);
		gaussianEnergy = energy(curvatureGaussian);

		//Calculate Mean curvature		
		Eigen::SparseMatrix<double> cotangentMatrix;
		igl::cotmatrix(V, F, cotangentMatrix);

		////这两行为什么被注释掉了？absolute?
		//igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, mass);
		//igl::invert_diag(mass, invMass);
		// Laplace-Beltrami of position
		meanCurvatureNormal = -/*invMass*/(cotangentMatrix*V);
		//This is the absolute value!
		curvatureAbsoluteMean = (meanCurvatureNormal.rowwise().norm());
		meanEnergy = energy(curvatureAbsoluteMean);
	}

	//Call analyze() to update this value
	double getGaussianCurvatureIntegralOverArea() const { return gaussianEnergy; }



	//Call analyze() to update this value
	double getMeanCurvatureIntegralOverArea() const { return meanEnergy; }

private:

	double energy(Eigen::MatrixXd & A)
	{
		double e = 0;
		for (int i = 0; i < A.rows(); ++i)
		{
			if (!std::isnan(A(i, 0)))
			{
				e += fabs(A(i, 0));
			}
		}
		return e;
	}

	Eigen::MatrixXd & V;
	Eigen::MatrixXi & F;

	Eigen::MatrixXd curvatureGaussian;
	Eigen::MatrixXd curvatureAbsoluteMean, meanCurvatureNormal;

	double gaussianEnergy, meanEnergy;
};