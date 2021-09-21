#pragma once

namespace std
{
#include <stdint.h>
}

#include "combinations.h"

#include <Eigen/Dense>

#include <vector>
#include <set>



using namespace std;



struct GCF
{
	double timestep = 0.5;
	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd& V, Eigen::MatrixXi& F, vector<int>& domain,
		Eigen::MatrixXd& VNormal, vector<vector<int> >& NeighborList,
		Eigen::MatrixXd& derivative)
	{
		//compute the distance to tangent plane of its neighbors
		int N = domain.size();  ///number of the vertex in the same color(domain)
		int V_index;
		Eigen::Vector3d Pos0, pos, pos_normal, neighbor, neighbor_normal, L, L_normal, d, d_tmp;
		double Energy = 100, Etemp;

		//set<int>::iterator it;

		for (int i = 0; i < N; ++i)
		{
			V_index = domain[i];
			Pos0 = V0.row(V_index);
			pos = V.row(V_index);
			pos_normal = VNormal.row(V_index);
			pos_normal.normalize();

			d.x() = 0;
			d.y() = 0;
			d.z() = 0;

			bool found = false;
			if (NeighborList[V_index].size() >= 3)
			{
				vector<int> neighbors(NeighborList[V_index].begin(), NeighborList[V_index].end());
				vector<Vector3d> neighbors_Vnormal, neighbors_P0P;
				vector<double>neighbors_P0DistanceTAN;

				Vector3d laplace_neighborse = Vector3d::Zero();

				for (int i = 0; i < neighbors.size(); i++)
				{

					L_normal += V.row(neighbors[i]);


				}
				L_normal = L_normal / neighbors.size() - pos;
				if (!L_normal.isZero())
					L_normal.normalize();

				for (int x = 0; x < neighbors.size(); x++)
					for (int y = x + 1; y < neighbors.size(); y++)
						for (int z = y + 1; z < neighbors.size(); z++)
						{
							Eigen::Vector3d a = V.row(neighbors[x]);
							Eigen::Vector3d b = V.row(neighbors[y]);
							Eigen::Vector3d c = V.row(neighbors[z]);

							auto tri_normal = (a - b).cross(c - b);
							if (!tri_normal.isZero())
								tri_normal.normalize();
							double proj = (pos - a).dot(tri_normal);
							//if (proj < 0 || proj > 1)
							//	return false;
						/*	if (fabs(proj) <= 1)*/
							d_tmp = timestep * proj * L_normal;


//#ifdef USE_ENERGY_VERTEX
							Etemp = (Pos0 - (pos + d_tmp)).norm() + d.norm();
//#endif

//#ifdef USE_ENERGY_DISPOSITION
							//Etemp = d_tmp.norm();

//#endif

							if (!found || Etemp < Energy) ///seek the smalleast
							{
								found = true;
								Energy = Etemp;
								d = d_tmp;
							}


						}



#ifdef USE_LINE_PROJ1
				for_each_combination(neighbors.begin(), neighbors.begin() + 2, neighbors.end(), [&](vector<int>::iterator first, vector<int>::iterator last)
					{
						Eigen::Vector3d a = V.row(*first);
						Eigen::Vector3d b = V.row(*(first + 1));

						Eigen::Vector3d o = pos - a;
						Eigen::Vector3d ab = b - a;
						double proj = 0;
						if (!ab.isZero())
							double proj = o.dot(ab) / ab.squaredNorm();
						//proj = o.dot(ab) / ab.norm();

						if (proj < 0 || proj > 1)
							return false;

						d_tmp = proj * ab - o;
						//		d_tmp = d_tmp.dot(pos_normal) * pos_normal;

								//d_tmp *= timestep;
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

						return false;
					});
#endif

			}
			else
			{
				derivative.row(V_index).setZero();
				continue;
			}
		
		} //for every vertex



		derivative.row(V_index) = d;
	} //calculateDerivative()
};




struct GCF_WM
{
	double timestep = 1;
	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd& V, Eigen::MatrixXi& F, vector<int>& domain,
		Eigen::MatrixXd& VNormal, vector<vector<int> >& NeighborList,
		Eigen::MatrixXd& derivative)
	{
		//compute the distance to tangent plane of its neighbors
		int N = domain.size();  ///number of the vertex in the same color(domain)
		int V_index;
		Eigen::Vector3d Pos0, pos, pos_normal, neighbor, neighbor_normal, L, L_normal, d, d_tmp;
		double Energy = 1e10, Etemp;

		//set<int>::iterator it;

		for (int i = 0; i < N; ++i)
		{
			V_index = domain[i];
			Pos0 = V0.row(V_index);
			pos = V.row(V_index);
			pos_normal = VNormal.row(V_index);
			pos_normal.normalize();

			d.x() = 0;
			d.y() = 0;
			d.z() = 0;

			bool found = false;
			if (NeighborList[V_index].size() >= 3)
			{
				vector<int> neighbors(NeighborList[V_index].begin(), NeighborList[V_index].end());
				vector<Vector3d> neighbors_Vnormal, neighbors_P0P;
				vector<double>neighbors_P0DistanceTAN;

				Vector3d laplace_neighborse = Vector3d::Zero();
				for (int j = 0; j < neighbors.size(); j++)
				{
					Vector3d neighborstemp, normaltemp;
					neighborstemp = V.row(neighbors[j]);
					laplace_neighborse += neighborstemp;
					Vector3d p0ptemp = neighborstemp - Pos0;

					normaltemp = VNormal.row(neighbors[j]);

					normaltemp.normalize();
					double tempdistance = p0ptemp.dot(normaltemp);
					neighbors_Vnormal.push_back(normaltemp);
					neighbors_P0P.push_back(p0ptemp);
					neighbors_P0DistanceTAN.push_back(abs(tempdistance));


				}
				laplace_neighborse = laplace_neighborse / neighbors.size() - pos;
				laplace_neighborse.normalize();
				double minD = 100.0;
				int minindex = 0;
				for (int k = 0; k < neighbors.size(); k++)
				{
					if (neighbors_P0DistanceTAN[k] < minD)
					{
						minD = neighbors_P0DistanceTAN[k];
						minindex = k;
					}
					else
					{
						continue;
					}
				}

				d_tmp = minD * laplace_neighborse;
				d_tmp *= timestep;

				d.norm();

				Etemp = (Pos0 - (pos + d_tmp)).norm() + d.norm();

				if (!found || Etemp < Energy) ///seek the smalleast
				{
					found = true;
					Energy = Etemp;
					d = d_tmp;
				}

			}



			//				std::vector<int> neighbors(NeighborList[V_index].begin(), NeighborList[V_index].end());
			//				for (int i = 0; i < neighbors.size(); i++)
			//				{
			//
			//					L_normal += V.row(neighbors[i]);
			//
			//
			//				}
			//				L_normal = L_normal / neighbors.size() - pos;
			//				if (!L_normal.isZero())
			//					L_normal.normalize();
			//
			//				for (int x = 0; x < neighbors.size(); x++)
			//					for (int y = x + 1; y < neighbors.size(); y++)
			//						for (int z = y + 1; z < neighbors.size(); z++)
			//						{
			//							Eigen::Vector3d a = V.row(neighbors[x]);
			//							Eigen::Vector3d b = V.row(neighbors[y]);
			//							Eigen::Vector3d c = V.row(neighbors[z]);
			//
			//							auto tri_normal = (a - b).cross(c - b);
			//							if (!tri_normal.isZero())
			//								tri_normal.normalize();
			//							double proj = (pos - a).dot(tri_normal);
			//							//if (proj < 0 || proj > 1)
			//							//	return false;
			//							if (fabs(proj) <= 1)
			//								d_tmp = timestep * proj * pos_normal;
			//
			//
			//#ifdef USE_ENERGY_VERTEX
			//							Etemp = (Pos0 - (pos + d_tmp)).norm() + d.norm();
			//#endif
			//
			//#ifdef USE_ENERGY_DISPOSITION
			//							Etemp = d_tmp.norm();
			//
			//
			//							if (!found || Etemp < Energy) ///seek the smalleast
			//							{
			//								found = true;
			//								Energy = Etemp;
			//								d = d_tmp;
			//							}
			//#endif
			//
			//
			//						}

						//}
			else
			{
				derivative.row(V_index).setZero();
				continue;
			}
			derivative.row(V_index) = d;
		} //for every vertex
	} //calculateDerivative()
};




//计算顶点邻域分法向量 类似于半Laplace思想；
vector<Vector3d> calculateVertexPartNormals(vector<Vector3d> P0Neighbors)
{
	const int NUM = P0Neighbors.size();
	vector<Vector3d> resultvector;
	for (int i = 0; i < NUM; i++)
	{
		Vector3d p1p0 = Vector3d::Zero();
		Vector3d p1p2 = Vector3d::Zero();
		int first = i;
		int second = i + 1;
		int third = i + 2;
		first = first % NUM;
		second = second % NUM;
		third = third % NUM;

		p1p0 = P0Neighbors[first] - P0Neighbors[second];
		p1p2 = P0Neighbors[third] - P0Neighbors[second];
		p1p0.normalize();
		p1p2.normalize();
		Vector3d tempnormal = Vector3d::Zero();
		tempnormal = p1p0.cross(p1p2);
		if (tempnormal.norm() < 0.01)//排除共线
			continue;
		else
		{
			tempnormal.normalize();
			resultvector.push_back(tempnormal);
		}


	}
	return resultvector;
};

struct GCF_WM_v2
{
	double timestep = 1.0;
	//const char* name = "ProjectiveGCFilter";
	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd& V, Eigen::MatrixXi& F, vector<int>& domain,
		Eigen::MatrixXd& VNormal, vector<vector<int> >& NeighborList,
		Eigen::MatrixXd& derivative)
	{
		//compute the distance to tangent plane of its neighbors
		int N = domain.size();  ///number of the vertex in the same color(domain)
		int V_index;
		Eigen::Vector3d Pos0, pos, pos_normal, neighbor, neighbor_normal, L, L_normal, d, d_tmp;
		double Energy = 1e10, Etemp;

		//set<int>::iterator it;

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


#ifdef USE_LINE_PROJ
				for_each_combination(neighbors.begin(), neighbors.begin() + 2, neighbors.end(), [&](vector<int>::iterator first, vector<int>::iterator last)
					{
						Eigen::Vector3d a = V.row(*first);
						Eigen::Vector3d b = V.row(*(first + 1));

						Eigen::Vector3d o = pos - a;
						Eigen::Vector3d ab = b - a;
						double proj = 0;
						if (!ab.isZero())
							double proj = o.dot(ab) / ab.squaredNorm();
						//proj = o.dot(ab) / ab.norm();

						if (proj < 0 || proj > 1)
							return false;

						d_tmp = proj * ab - o;
						//		d_tmp = d_tmp.dot(pos_normal) * pos_normal;

								//d_tmp *= timestep;
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

						return false;
					});
#endif



			}

			
			derivative.row(V_index) = d;
		} //for every vertex
	} //calculateDerivative()
};
 

struct ProjectiveGCFilter
{
	double timestep = 1.0;
	//const char* name = "ProjectiveGCFilter";
	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd& V, Eigen::MatrixXi& F, vector<int>& domain,
		Eigen::MatrixXd& VNormal, vector<vector<int> >& NeighborList,
		Eigen::MatrixXd& derivative)
	{
		//compute the distance to tangent plane of its neighbors
		int N = domain.size();  ///number of the vertex in the same color(domain)
		int V_index;
		Eigen::Vector3d Pos0, pos, pos_normal, neighbor, neighbor_normal, L, L_normal, d, d_tmp;
		double Energy = 1e10, Etemp;

		//set<int>::iterator it;

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


				//

				//				for (int i = 0; i < neighbors.size(); i++)
				//					for (int j = i + 1; j < neighbors.size(); j++)
				//						for (int k = j + 1; k < neighbors.size(); k++)
				//						{
				//Eigen::Vector3d a = V.row(neighbors[i]);
				//Eigen::Vector3d b = V.row(neighbors[j]);
				//Eigen::Vector3d c = V.row(neighbors[k]);
//
//#ifdef USE_TRI_PROJ			
//				for_each_combination(neighbors.begin(), neighbors.begin() + 3, neighbors.end(), [&](vector<int>::iterator first, vector<int>::iterator last)
//					{
//						Eigen::Vector3d a = V.row(*first);
//						Eigen::Vector3d b = V.row(*(first+1));
//						Eigen::Vector3d c = V.row(*(first+2));
//		
//
//						auto tri_normal = (b - a).cross(c - a);
//						tri_normal.normalize();
//						d_tmp = tri_normal.dot(pos_normal) * pos_normal;
//						//orthogonal projection onto plane
//						//d_tmp = -tri_normal * tri_normal.dot(pos - a);
//
//						//projection along vertex normal
//						//d_tmp = -tri_normal.dot(pos - a) / tri_normal.dot(pos_normal) * pos_normal;
//						//d_tmp *= timestep;
//#ifdef USE_ENERGY_VERTEX
//						Etemp = (Pos0 - (pos + d_tmp)).norm() + d_tmp.norm();
//#endif
//
//#ifdef USE_ENERGY_DISPOSITION
//						Etemp = d_tmp.norm();
//#endif
//
//
//						if (!found || Etemp < Energy) ///seek the smalleast					
//						{
//#ifdef USE_BARYCENTRIC_CHECK
//							Eigen::Vector3d proj = pos + d_tmp;
//
//							//Calculate barycentric coordinates
//							double bary1, bary2, bary3;
//							if (abs(tri_normal.x()) > abs(tri_normal.y()) && abs(tri_normal.x()) > abs(tri_normal.z()))
//							{
//								double denom = (b.y() - c.y()) * (a.z() - c.z()) + (c.z() - b.z()) * (a.y() - c.y());
//								bary1 = ((b.y() - c.y()) * (proj.z() - c.z()) + (c.z() - b.z()) * (proj.y() - c.y())) / denom;
//								bary2 = ((c.y() - a.y())) * (proj.z() - c.z()) + (a.z() - c.z()) * (proj.y() - c.y()) / denom;
//								bary3 = 1 - bary1 - bary2;
//							}
//							else if (abs(tri_normal.y()) > abs(tri_normal.z()))
//							{
//								double denom = (b.z() - c.z()) * (a.x() - c.x()) + (c.x() - b.x()) * (a.z() - c.z());
//								bary1 = ((b.z() - c.z()) * (proj.x() - c.x()) + (c.x() - b.x()) * (proj.z() - c.z())) / denom;
//								bary2 = ((c.z() - a.z())) * (proj.x() - c.x()) + (a.x() - c.x()) * (proj.z() - c.z()) / denom;
//								bary3 = 1 - bary1 - bary2;
//							}
//							else
//							{
//								double denom = (b.y() - c.y()) * (a.x() - c.x()) + (c.x() - b.x()) * (a.y() - c.y());
//								bary1 = ((b.y() - c.y()) * (proj.x() - c.x()) + (c.x() - b.x()) * (proj.y() - c.y())) / denom;
//								bary2 = ((c.y() - a.y())) * (proj.x() - c.x()) + (a.x() - c.x()) * (proj.y() - c.y()) / denom;
//								bary3 = 1 - bary1 - bary2;
//							}
//
//							//If the projection does not fall within the triangle, then it is invalid
//							if (bary1 < -0.05 || bary2 < -0.05 || bary3 < -0.05)
//								return false;
//#endif
//
//							found = true;
//							Energy = Etemp;
//							d = d_tmp;
//						}
//
//						return false;
//					});
//#endif
//


#ifdef USE_LINE_PROJ
				for_each_combination(neighbors.begin(), neighbors.begin() + 2, neighbors.end(), [&](vector<int>::iterator first, vector<int>::iterator last)
					{
						Eigen::Vector3d a = V.row(*first);
						Eigen::Vector3d b = V.row(*(first + 1));

						Eigen::Vector3d o = pos - a;
						Eigen::Vector3d ab = b - a;
						double proj = 0;
						if (!ab.isZero())
							double proj = o.dot(ab) / ab.squaredNorm();
						//proj = o.dot(ab) / ab.norm();

						if (proj < 0 || proj > 1)
							return false;

						d_tmp = proj * ab - o;
						//		d_tmp = d_tmp.dot(pos_normal) * pos_normal;

								//d_tmp *= timestep;
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

						return false;
					});
#endif



			}

			/*bool found = false;
			for (it = NeighborList[V_index].begin(); it != NeighborList[V_index].end(); it++)
			{
			neighbor_normal = VNormal.row(*it);
			neighbor = V.row(*it);

			L = neighbor - pos;
			L_normal = L.normalized();

			t = pos_normal.dot(neighbor_normal);
			if (fabs(t)<1E-4)
			{
			d_tmp = Eigen::Vector3d::Zero();
			}else
			{
			d_tmp = -neighbor_normal.dot(L_normal)/t*pos_normal;
			}

			if (!found || d_tmp.norm() < d.norm())
			{
			d = d_tmp;
			found = true;
			}
			}*/

			derivative.row(V_index) = d;
		} //for every vertex
	} //calculateDerivative()
};

