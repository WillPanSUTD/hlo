#pragma once

namespace std
{
#include <stdint.h>
}



#include <Eigen/Dense>

#include <vector>
#include <set>
#include "Mesh.h"
#include "Geometry.h"




using namespace std;
using namespace Eigen;


double angleCOS(Vector3d a, Vector3d b)
{
	return a.dot(b) / (a.norm()*b.norm());
}


struct HWLaplace
{
	double timestep;
	char* name = "HWLaplace";
	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd & V, Eigen::MatrixXi& F, vector<int> & domain,
		Eigen::MatrixXd & VNormal, vector<vector<int> > & NeighborList,
		Eigen::MatrixXd & derivative)
	{
		//compute the distance to tangent plane of its neighbors
		int N = domain.size();  ///number of the vertex in the same color(domain)
		int V_index;
		Eigen::Vector3d Pos0, pos, pos_normal, neighbor, neighbor_normal, L, L_normal, d, d_tmp1,d_tmp2;
		double Energy = 0, Etemp;

		vector<int>::iterator it;
		////per vertex
		for (int i = 0; i < N; ++i)
		{
			V_index = domain[i];
			Pos0 = V0.row(V_index);
			pos = V.row(V_index);
			pos_normal = VNormal.row(V_index);

			d.x() = 0;
			d.y() = 0;
			d.z() = 0;
			bool flag = true;
			bool found = false;
			if (NeighborList[V_index].size() >= 3)
			{
				std::vector<int> neighbors(NeighborList[V_index].begin(), NeighborList[V_index].end());
				int NN = neighbors.size();
				std::vector<Eigen::Vector3d> n_vec(NN);
				vector<double> angles(NN);

				////calculate Neighbourhood----Neighbour vector
				Eigen::Vector3d LO(0, 0, 0),X;
				for (size_t i = 0; i < NN; i++)
				{
					LO += V.row(neighbors.at(i));
					n_vec.at(i) = V.row(neighbors.at(i));
				}
				for (size_t i = 0; i < NN; i++)
				{
					//X = V.row(neighbors.at(i)) - LO;
					n_vec.at(i) -= LO;

				}

				//cos matrix M
				vector<int> ind(NN);
				vector<vector<double>> m(NN, vector<double>(NN, 0));
				for (int i = 0; i < NN; i++)
				{
					for (int j = 0; j < NN; j++)
						m[i][j] = angleCOS(n_vec.at(i), n_vec.at(j));

					ind[i] = min_element(m[i].begin(), m[i].end()) - m[i].begin();
				}


				for (int i = 0; i < NN; i++)
				{
					int xxx, yyy;
					if (i < ind[i])
					{
						xxx = i;
						yyy = ind[i];
					}
					else
					{
						xxx = ind[i];
						yyy = i;
					}

					Eigen::Vector3d center1(0, 0, 0);
					Eigen::Vector3d center2(0, 0, 0);
					for (size_t t=xxx; t <= yyy; t++)
						center1 += V.row(neighbors.at(t));
					for(size_t t=0;t<=xxx;t++)
						center2 += V.row(neighbors.at(t));
					for(size_t t=yyy;t<NN;t++)
						center2 += V.row(neighbors.at(t));
					Eigen::Vector3d move1 = center1 - pos;
					Eigen::Vector3d move2 = center2 - pos;
					move1.normalize();
					move2.normalize();

					d_tmp1 = timestep * move1.dot(pos_normal)*pos_normal;
					d_tmp1 = timestep * move2.dot(pos_normal)*pos_normal;

					//if (d_tmp.norm() < 1e-5)
					//{
					//	found = true;
					//	d = d_tmp;
					//	break;
					//}
#ifdef USE_BARYCENTRIC_CHECK
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
#endif


#ifdef USE_ENERGY_VERTEX
					Etemp = (Pos0 - (pos + d_tmp1)).norm() + d_tmp1.norm();
#endif

#ifdef USE_ENERGY_DISPOSITION
					Etemp = d_tmp.norm();
#endif
					if (flag && (!found || Etemp < Energy)) ///seek the smalleast
					{
						found = true;
						Energy = Etemp;
						d = d_tmp1;
					}



#ifdef USE_ENERGY_VERTEX
					Etemp = (Pos0 - (pos + d_tmp2)).norm() + d_tmp2.norm();
#endif

#ifdef USE_ENERGY_DISPOSITION
					Etemp = d_tmp.norm();
#endif
					if (flag && (!found || Etemp < Energy)) ///seek the smalleast
					{
						found = true;
						Energy = Etemp;
						d = d_tmp2;
					}

					derivative.row(V_index) = d;
				}

			}
			else
				cout << "on no.... not manifold..." << endl;
		}
	}
};


struct PointWiseLaplace
{
	double timestep;
	char* name = "PWLaplace";
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
				cout << "on no..... not manifold" << endl;
	}
}
};
