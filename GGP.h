#pragma once

namespace std
{
#include <stdint.h>
}



#include <Eigen/Dense>

#include <vector>
#include <set>


//#define USE_TRI_PROJ_GGP
//#define USE_LINE_PROJ_GGP
#define USE_HWLAPLACE_GGP
#define USE_HWLAPLACE_GGP_COS
#define USE_PWLAPLACE_GGP

using namespace std;


struct GGPFilter
{
	double timestep = 1.0;
	const char* name = "GGP";
	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd& V, Eigen::MatrixXi& F, vector<int>& domain,
		Eigen::MatrixXd& VNormal, vector<vector<int> >& NeighborList,
		Eigen::MatrixXd& derivative)
	{
		//compute the distance to tangent plane of its neighbors
		int N = domain.size();  ///number of the vertex in the same color(domain)
		int V_index;
		Eigen::Vector3d Pos0, pos, pos_normal, neighbor, neighbor_normal, L, L_normal, d, d_tmp;
		double Energy = 0, Etemp;



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
				////三角投影
#ifdef USE_TRI_PROJ_GGP		
				for_each_combination(neighbors.begin(), neighbors.begin() + 3, neighbors.end(), [&](vector<int>::iterator first, vector<int>::iterator last) {
					Eigen::Vector3d a = V.row(*first);
					Eigen::Vector3d b = V.row(*(first + 1));
					Eigen::Vector3d c = V.row(*(first + 2));

					auto tri_normal = (b - a).cross(c - a);
					tri_normal.normalize();

					//orthogonal projection onto plane
					//d_tmp = -tri_normal * tri_normal.dot(pos - a);

					//projection along vertex normal
					d_tmp = -tri_normal.dot(pos - a) / tri_normal.dot(pos_normal) * pos_normal;
					d_tmp *= timestep;

#ifdef USE_ENERGY_VERTEX
					Etemp = (Pos0 - (pos + d_tmp)).norm() + d_tmp.norm();
#endif

#ifdef USE_ENERGY_DISPOSITION
					Etemp = d_tmp.norm();
#endif


					if (!found || Etemp < Energy) ///seek the smalleast
					{

#ifdef USE_BARYCENTRIC_CHECK
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
							return false;
#endif

						found = true;
						Energy = Etemp;
						d = d_tmp;
					}

					return false;
					});
#endif
				////线投影
#ifdef USE_LINE_PROJ_GGP
				for_each_combination(neighbors.begin(), neighbors.begin() + 2, neighbors.end(), [&](vector<int>::iterator first, vector<int>::iterator last) {
					Eigen::Vector3d a = V.row(*first);
					Eigen::Vector3d b = V.row(*(first + 1));

					Eigen::Vector3d o = pos - a;
					Eigen::Vector3d ab = b - a;

					double proj = o.dot(ab) / ab.squaredNorm();
					if (proj < 0 || proj > 1)
						return false;

					d_tmp = proj * ab - o;
					d_tmp *= timestep;
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
					return false;
					});
#endif

#ifdef USE_UNKNOWN_GGP
				////不知道这个是什么鬼。。
				for (std::vector<int>::iterator it = NeighborList[V_index].begin(); it != NeighborList[V_index].end(); it++)
				{
					neighbor_normal = VNormal.row(*it);
					neighbor = V.row(*it);

					L = neighbor - pos;
					L_normal = L.normalized();

					double t = pos_normal.dot(neighbor_normal);
					if (fabs(t) < 1E-4)
					{
						d_tmp = Eigen::Vector3d::Zero();
						d_tmp *= timestep;
					}
					else
					{
						d_tmp = -neighbor_normal.dot(L_normal) / t * pos_normal;
						d_tmp *= timestep;
					}

#ifdef USE_ENERGY_VERTEX
					Etemp = (Pos0 - (pos + d_tmp)).norm() + d.norm();
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
				}
#endif
#ifdef USE_HWLAPLACE_GGP

				Vector3d pos0_normal = VNormal.row(V_index);//顶点法线

				pos0_normal.normalize();

				bool found = false;

				if (NeighborList[V_index].size() >= 3)
				{
					vector<int> neighbors(NeighborList[V_index].begin(), NeighborList[V_index].end());

					vector<Vector3d> neighbors_3dpoints;
					Vector3d laplace_neighbor = Vector3d::Zero();
					//vector<int>neighborse_index;
					//if (neighbors.size() == 0)
					//	cout << "jjdfafd";
					for (int i = 0; i < neighbors.size(); i++)
					{
						//Vector3d temp;
						//temp[0] = V.row(neighbors[i]).x();
						//temp[1] = V.row(neighbors[i]).y();
						//temp[2] = V.row(neighbors[i]).z();

						//laplace_neighbor[0] += temp[0];
						//laplace_neighbor[1] += temp[1];
						//laplace_neighbor[2] += temp[2];
						laplace_neighbor += V.row(neighbors[i]);
						neighbors_3dpoints.push_back(V.row(neighbors[i]));

						//neighborse_index.push_back(neighbors[i]);

					}
					laplace_neighbor = laplace_neighbor / neighbors.size() - pos;
					if (!laplace_neighbor.isZero())
						laplace_neighbor.normalize();


					double mindistance = getMinDistance(neighbors_3dpoints, pos0_normal, pos);
					//if (fabs(mindistance) > 1e3))
					//	cout << mindistance << endl;
					//cout << "mindistance:" << mindistance<< endl;
					Vector3d d_tmp = timestep * mindistance * laplace_neighbor;
					//Vector3d d_tmp = -1.0*timestep*mindistance*pos0_normal;
					//cout << "mindistance:" << mindistance << endl;
					//cout << "d_tmp:" << d_tmp << endl;
					//Etemp = d_tmp.norm();


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
				}
				else
					d_tmp.setZero();
#endif
				///半拉普拉斯COS
#ifdef USE_HWLAPLACE_GGP_COS
				int NN = neighbors.size();
				std::vector<Eigen::Vector3d> n_vec(NN);
				vector<double> angles(NN);

				////calculate Neighbourhood----Neighbour vector
				Eigen::Vector3d LO(0, 0, 0), X;
				for (size_t i = 0; i < NN; i++)
				{
					LO += V.row(neighbors.at(i));
					n_vec.at(i) = V.row(neighbors.at(i));
				}
				LO /= NN;
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

					int N1 = yyy - xxx + 1;
					int N2 = NN - N1 + 2;
					Eigen::Vector3d center1(0, 0, 0);
					Eigen::Vector3d center2(0, 0, 0);
					for (size_t t = xxx; t <= yyy; t++)
						center1 += V.row(neighbors.at(t));
					center1 /= N1;
					for (size_t t = 0; t <= xxx; t++)
						center2 += V.row(neighbors.at(t));
					for (size_t t = yyy; t < NN; t++)
						center2 += V.row(neighbors.at(t));
					center2 /= N2;
					Eigen::Vector3d move1 = center1 - pos;
					Eigen::Vector3d move2 = center2 - pos;

					d_tmp = timestep * move1.dot(pos_normal) * pos_normal;


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


					d_tmp = timestep * move2.dot(pos_normal) * pos_normal;
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
				}

#endif




#ifdef USE_HWLAPLACE_GGP_HN
				int N = neighbors.size();
				int N_half = N / 2 + 1;
				//if (N_half < 3)
				//	N_half = 3;
				for (std::vector<int>::iterator first = neighbors.begin(); first < neighbors.end(); first++)
				{
					Eigen::Vector3d center(0, 0, 0);

					for (size_t i = 0; i < N_half; i++)
					{
						int j = i;
						if ((first - neighbors.begin()) + i >= N)
							j = i - N;

						center += V.row(*(first + j));
					}
					center = center / N_half;
					Eigen::Vector3d move = center - pos;
					move.normalize();

					d_tmp = timestep * move.dot(pos_normal) * pos_normal;
					//int N_half = neighbors.size() / 2;
					//for (std::vector<int>::iterator first = neighbors.begin(); (first + N_half) < neighbors.end(); first++)
					//{
					//	Eigen::Vector3d center(0, 0, 0);
					//	//for (size_t i = 0; i < N_half; i++)
					//	//{
					//	//	center += V.row(*(first + i));
					//	//}

					//	for (size_t i = 0; i < N_half; i++)
					//	{
					//		center += V.row(*(first + i));
					//	}
					//	center = center / N_half;
					//	d_tmp = center - pos;
					//	d_tmp *= timestep;
					//

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
				}
#endif


#ifdef USE_PWLAPLACE_GGP
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


#endif


			}
			derivative.row(V_index) = d;
		} //for every vertex
	} //calculateDerivative()
};