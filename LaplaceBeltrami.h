#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/doublearea.h>
#include <igl/repdiag.h>
#include <igl/grad.h>

#include <vector>
#include <set>

using namespace std;
using namespace Eigen;


struct PointWiseLaplace
{
	double timestep;
	const char* name = "Pointwise Uniform Laplace";
	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd & V, Eigen::MatrixXi& F, vector<int> & domain,
		Eigen::MatrixXd & VNormal, vector<vector<int> > & NeighborList,
		Eigen::MatrixXd & derivative)
	{
		//compute the distance to tangent plane of its neighbors
		int N = domain.size();  ///number of the vertex in the same color(domain)
		int V_index;
		Eigen::Vector3d Pos0, pos, pos_normal, neighbor, neighbor_normal, L, L_normal, d, d_tmp;
		double Energy = 0, Etemp;
		double t;


		vector<int>::iterator it;

		for (int i = 0; i < N; ++i)
		{
			V_index = domain[i];
			Pos0 = V0.row(V_index);
			pos = V.row(V_index);
			pos_normal = VNormal.row(V_index);

			d.x() = 0;
			d.y() = 0;
			d.z() = 0;

			bool found = false;
			if (NeighborList[V_index].size() >= 3)
			{
				std::vector<int> neighbors(NeighborList[V_index].begin(), NeighborList[V_index].end());
				Eigen::Vector3d center(0, 0, 0);
				//for (size_t i = 0; i < N_half; i++)
				//{
				//	center += V.row(*(first + i));
				//}
				std::vector<int>::iterator first = neighbors.begin();
				for (size_t i = 0; i < neighbors.size(); i++)
				{
					center += V.row(*(first + i));
				}
				center = center / neighbors.size();
				d_tmp = center - pos;
				d_tmp *= timestep;
#ifdef USE_BARYCENTRIC_CHECK_HW
				if (N_half >= 3)
				{
					Eigen::Vector3d a = V.row(*first);
					Eigen::Vector3d b = V.row(*(first + 1));
					Eigen::Vector3d c = V.row(*(first + 2));
					auto tri_normal = (b - a).cross(c - a);
					tri_normal.normalize();

					Eigen::Vector3d proj = pos + d_tmp;
					//Calculate barycentric coordinates
					double bary1, bary2, bary3;
					if (abs(tri_normal.x()) > abs(tri_normal.y()) && abs(tri_normal.x()) > abs(tri_normal.z()))
					{
						double denom = (b.y() - c.y()) * (a.z() - c.z()) + (c.z() - b.z()) * (a.y() - c.y());
						bary1 = ((b.y() - c.y()) * (proj.z() - c.z()) + (c.z() - b.z()) * (proj.y() - c.y())) / denom;
						bary2 = ((c.y() - a.y())) * (proj.z() - c.z()) + (a.z() - c.z()) * (proj.y() - c.y()) / denom;
						bary3 = 1 - bary1 - bary2;
					}
					else if (abs(tri_normal.y()) > abs(tri_normal.z()))
					{
						double denom = (b.z() - c.z()) * (a.x() - c.x()) + (c.x() - b.x()) * (a.z() - c.z());
						bary1 = ((b.z() - c.z()) * (proj.x() - c.x()) + (c.x() - b.x()) * (proj.z() - c.z())) / denom;
						bary2 = ((c.z() - a.z())) * (proj.x() - c.x()) + (a.x() - c.x()) * (proj.z() - c.z()) / denom;
						bary3 = 1 - bary1 - bary2;
					}
					else
					{
						double denom = (b.y() - c.y()) * (a.x() - c.x()) + (c.x() - b.x()) * (a.y() - c.y());
						bary1 = ((b.y() - c.y()) * (proj.x() - c.x()) + (c.x() - b.x()) * (proj.y() - c.y())) / denom;
						bary2 = ((c.y() - a.y())) * (proj.x() - c.x()) + (a.x() - c.x()) * (proj.y() - c.y()) / denom;
						bary3 = 1 - bary1 - bary2;
					}

					//If the projection does not fall within the triangle, then it is invalid
					if (bary1 < -0.05 || bary2 < -0.05 || bary3 < -0.05)
						flag = false;
				}
#endif

#ifdef USE_ENERGY_VERTEX
				Etemp = (Pos0 - (pos + d_tmp)).norm() + d_tmp.norm();
#endif

#ifdef USE_ENERGY_DISPOSITION
				Etemp = d_tmp.norm();
#endif


				if (!found || Etemp < Energy) ///seek the smalleast
				{
					found = true;
					Energy = Etemp;
					d = d_tmp;
				}

				derivative.row(V_index) = d;

			}
			else
				derivative.row(V_index) = d;
		}
	}
};




//////all functions dependent on libigl
///Previous version... not make sense
struct LaplaceBeltrami
{
	const char* name = "LB";
	double timestep=1;

	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd & V, Eigen::MatrixXi& F, vector<int> & domain, Eigen::MatrixXd & VNormal, vector<vector<int> > & NeighborList, Eigen::MatrixXd & derivative)
	{
		Eigen::SparseMatrix<double> cotangentMatrix, mass, invMass;
		igl::cotmatrix(V0, F, cotangentMatrix);

		igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, mass);
		igl::invert_diag(mass, invMass);
		derivative = timestep * invMass * (cotangentMatrix*V);
	}
};

///Updated by Wei Pan, cotangent Laplace-Beltrami
///derivative is initial V0. 
///Laplasian Operator: cotangentMatrix
struct LaplaceBeltrami_v1
{
	const char* name = "cotLB";
	double timestep = 0.1;

	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd & V, Eigen::MatrixXi& F, vector<int> & domain,
		Eigen::MatrixXd & VNormal, vector<vector<int> > & NeighborList, Eigen::MatrixXd & derivative)
	{
		// Compute Laplace-Beltrami operator: #V by #V
		Eigen::SparseMatrix<double> cotangentMatrix, mass, invMass;
		igl::cotmatrix(V0, F, cotangentMatrix);
		Eigen::MatrixXd U = V;

		// Recompute just mass matrix on each step
		igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

		// Solve (M-delta*L) U = M*U
		const auto & S = (mass - timestep * cotangentMatrix);
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
		assert(solver.info() == Eigen::Success);
		U = solver.solve(mass*U).eval();

		// Compute centroid and subtract (also important for numerics)
		VectorXd dblA;
		igl::doublearea(U, F, dblA);
		// Compute centroid and subtract (also important for numerics)
		double area = 0.5*dblA.sum();

		MatrixXd BC;
		igl::barycenter(U, F, BC);
		RowVector3d centroid(0, 0, 0);
		for (int i = 0; i < BC.rows(); i++)
		{
			centroid += 0.5*dblA(i) / area * BC.row(i);
		}
		U.rowwise() -= centroid*timestep;
		// Normalize to unit surface area (important for numerics)
		U.array() /= sqrt(area);
		derivative = U - V;
		V = U;

		//igl::invert_diag(mass, invMass);
		//derivative = timestep * invMass * (cotangentMatrix*V);
	}
};


///Laplasian Operator: K
struct LaplaceBeltrami_v2
{
	const char* name = "LBv2";
	double timestep = 1;

	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd & V, Eigen::MatrixXi& F, vector<int> & domain,
		Eigen::MatrixXd & VNormal, vector<vector<int> > & NeighborList, Eigen::MatrixXd & derivative)
	{
		// Compute Laplace-Beltrami operator: #V by #V
		Eigen::SparseMatrix<double>  mass, invMass;

		// Alternative construction of same Laplacian
		SparseMatrix<double> G, K;
		// Gradient/Divergence
		igl::grad(V0, F, G);
		// Diagonal per-triangle "mass matrix"
		VectorXd dblA;
		igl::doublearea(V0, F, dblA);
		// Place areas along diagonal #dim times
		const auto & T = 1.*(dblA.replicate(3, 1)*0.5).asDiagonal();
		// Laplacian K built as discrete divergence of gradient or equivalently
		// discrete Dirichelet energy Hessian
		K = -G.transpose() * T * G;

		Eigen::MatrixXd U = V;

		// Recompute just mass matrix on each step
		igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);

		// Solve (M-delta*L) U = M*U
		const auto & S = (mass - timestep * K);
		Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
		assert(solver.info() == Eigen::Success);
		U = solver.solve(mass*U).eval();

		// Compute centroid and subtract (also important for numerics)

		igl::doublearea(U, F, dblA);
		// Compute centroid and subtract (also important for numerics)
		double area = 0.5*dblA.sum();

		MatrixXd BC;
		igl::barycenter(U, F, BC);
		RowVector3d centroid(0, 0, 0);
		for (int i = 0; i < BC.rows(); i++)
		{
			centroid += 0.5*dblA(i) / area * BC.row(i);
		}

		U.rowwise() -= centroid;
		// Normalize to unit surface area (important for numerics)
		U.array() /= sqrt(area);

		derivative = U- V;

		V = U;

		//igl::invert_diag(mass, invMass);
		//derivative = timestep * invMass * (cotangentMatrix*V);
	}
};

