#pragma once

namespace std
{
#include <stdint.h>
}



#include <Eigen/Dense>

#include <vector>
#include <set>

#include <time.h>


using namespace std;
using namespace Eigen;


double angleCOS(Vector3d a, Vector3d b)
{
	return a.dot(b) / (a.norm() * b.norm());
}

double area(Vector3d a, Vector3d b)
{
	return a.cross(b).norm();
}



double getMinDistance(std::vector<Vector3d> &points, Vector3d &normal, Vector3d& pos0)
{

	Vector4d par;//点法式平面方程(顶点与邻域平面法线组成的第一个平面点法式方程:A(x-x0)+B(y-y0)+C(z-z0)=0)
	par[0] = normal.x();
	par[1] = normal.y();
	par[2] = normal.z();
	par[3] = -1.0 * (pos0.x() * normal.x() + pos0.y() * normal.y() + pos0.z() * normal.z());



	Vector3d point1, p0p1, normal1, p0p2, normal2, p0p3, normal3, p0p4, normal4;
	if (abs(par[2]) < 0.000001)
	{
		point1[2] = 0.0;
		point1[0] = 1.0;
		point1[1] = (-1.0 * (point1[0] * par[0] + par[3])) / par[1];

	}
	else if (abs(par[1]) < 0.000001)
	{
		point1[1] = 0.0;
		point1[0] = 1.0;
		point1[2] = (-1.0 * (point1[0] * par[0] + par[3])) / par[2];


	}
	else if (abs(par[0]) < 0.000001)
	{
		point1[0] = 0.0;
		point1[1] = 1.0;
		point1[2] = (-1.0 * (point1[1] * par[1] + par[3])) / par[2];
	}
	else
	{
		point1[0] = 1.0;
		point1[1] = 1.0;
		point1[2] = (-1.0 * (par[0] * point1[0] + par[1] * point1[1] + par[3])) / par[2];
	}

	p0p1 = point1 - pos0;
	p0p1.normalize();
	normal1 = p0p1.cross(normal);
	p0p2 = normal1 + p0p1;
	normal2 = p0p2.cross(normal);
	normal2.normalize();
	normal3 = p0p1;
	p0p4 = -1.0 * p0p1 + normal1;
	normal4 = p0p4.cross(normal);
	normal4.normalize();
	//////////切平面构造
	vector<Vector3d>normallist;
	normallist.push_back(normal1);
	normallist.push_back(normal2);
	normallist.push_back(normal3);
	normallist.push_back(normal4);
	/*cout << "0-1:" << acos(normallist[0].dot( normallist[1]))*180.0/3.14 << endl;
	cout << "0-2:" << acos(normallist[0].dot( normallist[2]))*180.0 / 3.14 << endl;
	cout << "0-3:" << acos(normallist[0].dot( normallist[3]))*180.0/3.14 << endl;*/

	vector<double>result_distance;
	for (int i = 0; i < 4; i++)
	{
		Vector4d ABCDtemp;
		ABCDtemp[0] = normallist[i][0];
		ABCDtemp[1] = normallist[i][1];
		ABCDtemp[2] = normallist[i][2];
		ABCDtemp[3] = -1.0 * (pos0.x() * normallist[i].x() + pos0.y() * normallist[i].y() + pos0.z() * normallist[i].z());

		std::vector<Vector3d>points1, points2;//通过平面切出来的二分类邻域点
		for (int j = 0; j < points.size(); j++)
		{
			double flag = ABCDtemp[0] * points[j][0] + ABCDtemp[1] * points[j][1] + ABCDtemp[2] * points[j][2] + ABCDtemp[3];//判决条件，通过点在平面方程的值
			if (flag >= 0.0)
			{
				points1.push_back(points[j]);
			}
			else
			{
				points2.push_back(points[j]);

			}
		}




		double projection1 = 0.0;
		if (points1.size() > 0)
		{

			Vector3d plane1_laplace1 = Vector3d::Zero();//HLO
			for (int i = 0; i < points1.size(); i++)
			{
				plane1_laplace1[0] += points1[i][0];
				plane1_laplace1[1] += points1[i][1];
				plane1_laplace1[2] += points1[i][2];

			}
			plane1_laplace1 = plane1_laplace1 / points1.size();

			plane1_laplace1 = plane1_laplace1 - pos0;
			projection1 = plane1_laplace1.dot(normal);
			result_distance.push_back(abs(projection1));

		}
		double projection2 = 0.0;
		if (points2.size() > 0)
		{
			Vector3d plane1_laplace2 = Vector3d::Zero();//HLO
			for (int i = 0; i < points2.size(); i++)
			{
				plane1_laplace2[0] += points2[i][0];
				plane1_laplace2[1] += points2[i][1];
				plane1_laplace2[2] += points2[i][2];

			}
			plane1_laplace2 = plane1_laplace2 / points2.size();
			plane1_laplace2 = plane1_laplace2 - pos0;
			projection2 = plane1_laplace2.dot(normal);
			result_distance.push_back(abs(projection2));
		}






		//平面法线投影
		/*double projection1 = 0.0;
		if (points1.size() >= 3)
		{
			vector<cv::Point3d> plane1_points;
			Vector3d plane1_laplace1 = Vector3d::Zero();//半拉布拉斯算子
			for (int i = 0; i < points1.size(); i++)
			{
				cv::Point3d temp;
				temp.x = points1[i][0];
				temp.y = points1[i][1];
				temp.z = points1[i][2];

				plane1_laplace1[0] += points1[i][0];
				plane1_laplace1[1] += points1[i][1];
				plane1_laplace1[2] += points1[i][2];

				plane1_points.push_back(temp);

			}
			plane1_laplace1 = plane1_laplace1 / points1.size();

			plane1_laplace1 = plane1_laplace1 - pos0;
			cv::Mat result1 = getPlane(plane1_points);
			Vector3d plane1_points_normal;//半平面法线

			plane1_points_normal[0] = result1.at<double>(0, 0);
			plane1_points_normal[1] = result1.at<double>(1, 0);
			plane1_points_normal[2] = -1.0;
			plane1_points_normal.normalize();
			projection1 = plane1_laplace1.dot(plane1_points_normal);
			//result_distance.push_back(abs(projection1));

		}
		double projection2 = 0.0;
		if (points2.size() >= 3)
		{
			vector<cv::Point3d> plane2_points;
			Vector3d plane1_laplace2 = Vector3d::Zero();//半拉布拉斯算子
			for (int i = 0; i < points2.size(); i++)
			{
				cv::Point3d temp;
				temp.x = points2[i][0];
				temp.y = points2[i][1];
				temp.z = points2[i][2];

				plane1_laplace2[0] += points2[i][0];
				plane1_laplace2[1] += points2[i][1];
				plane1_laplace2[2] += points2[i][2];

				plane2_points.push_back(temp);

			}
			plane1_laplace2 = plane1_laplace2 / points2.size();
			plane1_laplace2 = plane1_laplace2 - pos0;

			cv::Mat result2 = getPlane(plane2_points);
			Vector3d plane1_points_norma2;
			plane1_points_norma2[0] = result2.at<double>(0, 0);
			plane1_points_norma2[1] = result2.at<double>(1, 0);
			plane1_points_norma2[2] = -1.0;
			plane1_points_norma2.normalize();
			projection2 = plane1_laplace2.dot(plane1_points_norma2);
			//result_distance.push_back(abs(projection2));
		}

		if (points1.size() >= 3 && points2.size() >= 3)
		{
			result_distance.push_back(abs(projection1));
			result_distance.push_back(abs(projection2));
		}
		else
		{
			continue;
		}*/


	}
	double ReMinDistance = 100.0;
	if (result_distance.size())
	{
		sort(result_distance.begin(), result_distance.end());
		ReMinDistance = result_distance[0];
		//for (int i = 0; i < result_distance.size(); i++)
		//{
		//	if (result_distance[i] < ReMinDistance)
		//	{
		//		ReMinDistance = result_distance[i];
		//	}
		//	else
		//	{
		//		continue;
		//	}
		//}

	}
	else
	{
		ReMinDistance = 0.0;
	}
	return ReMinDistance;


}

///
struct HWLaplace_4Cuts
{
	double timestep;
	const char* name = "HWLaplace";
	void calculateDerivative(const Eigen::MatrixXd & V0, Eigen::MatrixXd& V, Eigen::MatrixXi& F, vector<int>& domain,
		Eigen::MatrixXd& VNormal, vector<vector<int> >& NeighborList,
		Eigen::MatrixXd& derivative)
	{
		//compute the distance to tangent plane of its neighbors
		int N = domain.size();  ///number of the vertex in the same color(domain)
		int V_index;
		Eigen::Vector3d Pos0, pos, pos_normal;
		bool flag = true;
		double Energy = 1e10, Etemp;
	
		for (int i = 0; i < N; ++i)
		{
			V_index = domain[i];
			//cout << "V_index:" << V_index << endl;
			Pos0 = V0.row(V_index);
			pos = V.row(V_index);
			//pos_normal = VNormal.row(V_index);

			//vector<int> neighbors_face;//顶点与邻域点组成的三角面片
			//vector<double> neighbors_face_area;//顶点与邻域点组成的三角面片
			//vector<Vector3d> neighbors_face_normals;//顶点与邻域点组成的面的法线

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
				Etemp = (Pos0 - (pos + d_tmp)).norm() + d_tmp.norm();

				if (!found || Etemp < Energy) ///seek the smalleast
				{
					found = true;
					Energy = Etemp;
				}

				//if (fabs(d_tmp.norm()) > 1e5)
				//{
				//	cout << Etemp << endl;
				//	cout << laplace_neighbor.x() << laplace_neighbor.y() << laplace_neighbor.z() << endl;
				//	cout << mindistance << endl;

				//}
				//cout << d_tmp.norm() << endl;

				derivative.row(V_index) = d_tmp;
				//if (isnan(derivative.row(V_index).x()))
				//	cout << "dfafd" << endl;

				/*
								vector<vector<int>>index = getMinDistancePoints(neighborse_index,neighbors_3dpoints, pos0_normal, Pos0);
								vector<Vector3d>laplace_points_distance;

								if (index.size() == 2) {
									//cout << "index.size() == 2" << endl;
									vector<double> area;
									for (int i = 0; i < index.size(); i++)
									{
										Vector3d la_point = Vector3d::Zero();
										Vector3d AB = Vector3d::Zero();
										Vector3d AC = Vector3d::Zero();
										for (int j = 0; j < index[i].size(); j++)
										{
											la_point[0] += V.row(index[i][j]).x();
											la_point[1] += V.row(index[i][j]).y();
											la_point[2] += V.row(index[i][j]).z();

											if (j == 0)
											{
												AB = la_point - pos;

											}
											if (j == index[i].size() - 1)
											{
												AC = la_point - pos;
											}

										}
										area.push_back(AB.cross(AC).norm());
										la_point = la_point / index[i].size();
										laplace_points_distance.push_back(la_point);
									}
									Vector3d move1 = laplace_points_distance[0] - pos;
									Vector3d move2 = laplace_points_distance[1] - pos;
									Vector3d d_tmp1 = timestep * move1.dot(pos_normal)*pos_normal;
									Vector3d d_tmp2 = timestep * move2.dot(pos_normal)*pos_normal;

									if (area[0 > area[1]]) {
				#ifdef USE_ENERGY_VERTEX
										Etemp = (Pos0 - (pos + d_tmp1)).norm() + d_tmp1.norm();
				#endif

				#ifdef USE_ENERGY_DISPOSITION
										Etemp = d_tmp1.norm();
				#endif
										if (flag && (!found || Etemp < Energy)) ///seek the smalleast
										{
											found = true;
											Energy = Etemp;
											d = d_tmp1;
										}


									}
									else
									{


				#ifdef USE_ENERGY_VERTEX
										Etemp = (Pos0 - (pos + d_tmp2)).norm() + d_tmp2.norm();
				#endif

				#ifdef USE_ENERGY_DISPOSITION
										Etemp = d_tmp2.norm();
				#endif
										if (flag && (!found || Etemp < Energy)) ///seek the smalleast
										{
											found = true;
											Energy = Etemp;
											d = d_tmp2;
										}
									}

								derivative.row(V_index) = d;

								}
								else
								{
									Vector3d center = Vector3d::Zero();

									Vector3d d_tmp = Vector3d::Zero();
									for (int i = 0; i < neighbors.size(); i++)
									{
										Vector3d temp;
										temp[0] = V.row(neighbors[i]).x();
										temp[1] = V.row(neighbors[i]).y();
										temp[2] = V.row(neighbors[i]).z();

										center[0] += temp[0];
										center[1] += temp[1];
										center[2] += temp[2];

									}
									center = center / neighbors.size();
									d_tmp = center - pos;
									d_tmp *= timestep;
				#ifdef USE_ENERGY_VERTEX
									Etemp = (Pos0 - (pos + d_tmp)).norm() + d_tmp.norm();
				#endif
									if (!found || Etemp < Energy) ///seek the smalleast
									{
										found = true;
										Energy = Etemp;
										d = d_tmp;
									}

									derivative.row(V_index) = d;

								}
								*/
			}
			else
			{
				//cout << "dfjadf" << endl;
				derivative.row(V_index).setZero();

			}

			//cout << derivative.row(V_index).norm() << endl;

		}

	}
};


//
//
/////cos\alpha = -1, Weighted by Area (not ready yet)
//struct HWLaplaceCOS_AW
//{
//	double timestep;
//	char* name = "HWLaplace";
//	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd & V, Eigen::MatrixXi& F, vector<int> & domain,
//		Eigen::MatrixXd & VNormal, vector<vector<int> > & NeighborList,
//		Eigen::MatrixXd & derivative)
//	{
//		//compute the distance to tangent plane of its neighbors
//		int N = domain.size();  ///number of the vertex in the same color(domain)
//		int V_index;
//		Eigen::Vector3d Pos0, pos, pos_normal, neighbor, neighbor_normal, L, L_normal, d, d_tmp;
//		double Energy = 0, Etemp;
//
//		vector<int>::iterator it;
//		////per vertex
//		for (int i = 0; i < N; ++i)
//		{
//			V_index = domain[i];
//			Pos0 = V0.row(V_index);
//			pos = V.row(V_index);
//			pos_normal = VNormal.row(V_index);
//
//			d.x() = 0;
//			d.y() = 0;
//			d.z() = 0;
//			bool flag = true;
//			bool found = false;
//			if (NeighborList[V_index].size() < 3)
//			{
//				cout << "on no.... not manifold..." << endl;
//				continue;
//
//			}
//			else
//			{
//				std::vector<int> neighbors(NeighborList[V_index].begin(), NeighborList[V_index].end());
//				int NN = neighbors.size();
//				std::vector<Eigen::Vector3d> n_vec(NN);
//				vector<double> angles(NN);
//
//				////calculate Neighbourhood----Neighbour vector
//				Eigen::Vector3d LO(0, 0, 0), X;
//				for (size_t i = 0; i < NN; i++)
//				{
//					LO += V.row(neighbors.at(i));
//					n_vec.at(i) = V.row(neighbors.at(i));
//				}
//				LO /= NN;
//				for (size_t i = 0; i < NN; i++)
//				{
//					//X = V.row(neighbors.at(i)) - LO;
//					n_vec.at(i) -= LO;
//
//				}
//
//				//cos matrix M
//				vector<int> ind(NN);
//				vector<vector<double>> m(NN, vector<double>(NN, 0));
//				for (int i = 0; i < NN; i++)
//				{
//					for (int j = 0; j < NN; j++)
//						m[i][j] = angleCOS(n_vec.at(i), n_vec.at(j));
//
//					ind[i] = min_element(m[i].begin(), m[i].end()) - m[i].begin();
//				}
//
//
//				for (int i = 0; i < NN; i++)
//				{
//					int xxx, yyy;
//					if (i < ind[i])
//					{
//						xxx = i;
//						yyy = ind[i];
//					}
//					else
//					{
//						xxx = ind[i];
//						yyy = i;
//					}
//
//					int N1 = yyy - xxx + 1;
//					int N2 = NN - N1 + 2;
//					Eigen::Vector3d center1(0, 0, 0);
//					Eigen::Vector3d center2(0, 0, 0);
//
//
//					for (size_t t = xxx; t <= yyy; t++)
//						center1 += V.row(neighbors.at(t));
//					center1 /= N1;
//					for (size_t t = 0; t <= xxx; t++)
//						center2 += V.row(neighbors.at(t));
//					for (size_t t = yyy; t < NN; t++)
//						center2 += V.row(neighbors.at(t));
//					center2 /= N2;
//					Eigen::Vector3d move1 = center1 - pos;
//					Eigen::Vector3d move2 = center2 - pos;
//
//					d_tmp = timestep * move1.dot(pos_normal)*pos_normal;
//#ifdef USE_ENERGY_VERTEX
//					Etemp = (Pos0 - (pos + d_tmp)).norm() + d_tmp.norm();
//#endif
//
//#ifdef USE_ENERGY_DISPOSITION
//					Etemp = d_tmp.norm();
//#endif
//					if (!found || Etemp < Energy) ///seek the smalleast
//					{
//						found = true;
//						Energy = Etemp;
//						d = d_tmp;
//					}
//
//
//					d_tmp = timestep * move2.dot(pos_normal)*pos_normal;
//#ifdef USE_ENERGY_VERTEX
//					Etemp = (Pos0 - (pos + d_tmp)).norm() + d_tmp.norm();
//#endif
//
//#ifdef USE_ENERGY_DISPOSITION
//					Etemp = d_tmp.norm();
//#endif
//					if (!found || Etemp < Energy) ///seek the smalleast
//					{
//						found = true;
//						Energy = Etemp;
//						d = d_tmp;
//					}
//
//					derivative.row(V_index) = d;
//				}
//
//			}
//
//		}
//	}
//};




struct HWLaplaceCOS
{
	double timestep=1;
	const char* name = "HWLaplace";
	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd& V, Eigen::MatrixXi& F, vector<int>& domain,
		Eigen::MatrixXd& VNormal, vector<vector<int> >& NeighborList,
		Eigen::MatrixXd& derivative)
	{
		//compute the distance to tangent plane of its neighbors
		int N = domain.size();  ///number of the vertex in the same color(domain)
		int V_index;
		Eigen::Vector3d Pos0, pos, pos_normal, neighbor, neighbor_normal, L, L_normal,d_tmp;
		double Energy = 0, Etemp;

		vector<int>::iterator it;
		////per vertex
		for (int i = 0; i < N; ++i)
		{
			V_index = domain[i];
			Pos0 = V0.row(V_index);
			pos = V.row(V_index);
			//pos_normal = VNormal.row(V_index);

			bool flag = true;
			bool found = false;
			if (NeighborList[V_index].size() < 3)
			{
				derivative.row(V_index).setZero();
				//cout << "on no.... not manifold..." << endl;
				continue;
			}
			else
			{
				std::vector<int> neighbors(NeighborList[V_index].begin(), NeighborList[V_index].end());
				int NN = neighbors.size();
				std::vector<Eigen::Vector3d> n_vec(NN);
				vector<double> angles(NN);

				////calculate Neighbourhood----Neighbour vector
				Eigen::Vector3d LO(0, 0, 0);
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

					LO.normalize();

					d_tmp = timestep * move1.dot(LO) * LO;

					//d_tmp = timestep * move1.dot(pos_normal)*pos_normal;
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
					}

					derivative.row(V_index).setZero();
				}

			}

		}
	}
};




struct HWLaplaceCOS_NP
{
	double timestep;
	const char* name = "HWLaplace";
	void calculateDerivative(const Eigen::MatrixXd & V0, Eigen::MatrixXd& V, Eigen::MatrixXi& F, vector<int>& domain,
		Eigen::MatrixXd& VNormal, vector<vector<int> >& NeighborList,
		Eigen::MatrixXd& derivative)
	{
		//compute the distance to tangent plane of its neighbors
		int N = domain.size();  ///number of the vertex in the same color(domain)
		int V_index;
		Eigen::Vector3d Pos0, pos, pos_normal, neighbor, neighbor_normal, L, L_normal, d, d_tmp;
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
			if (NeighborList[V_index].size() < 3)
			{
				derivative.row(V_index) = d;
				//cout << "on no.... not manifold..." << endl;
				continue;
			}
			else
			{
				std::vector<int> neighbors(NeighborList[V_index].begin(), NeighborList[V_index].end());
				int NN = neighbors.size();
				std::vector<Eigen::Vector3d> n_vec(NN);
				vector<double> angles(NN);

				////calculate Neighbourhood----Neighbour vector
				Eigen::Vector3d LO(0, 0, 0);
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

					LO.normalize();

					d_tmp = timestep * move1;

					//d_tmp = timestep * move1.dot(pos_normal)*pos_normal;
#ifdef USE_ENERGY_VERTEX
					Etemp = /*(Pos0 - (pos + d_tmp)).norm() + */d_tmp.norm();
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

					derivative.row(V_index) = d;
				}

			}

		}
	}
};




struct HWLaplaceHN
{
	double timestep;
	const char* name = "HWLaplace";
	void calculateDerivative(const Eigen::MatrixXd  V0, Eigen::MatrixXd& V, Eigen::MatrixXi& F, vector<int>& domain,
		Eigen::MatrixXd& VNormal, vector<vector<int> >& NeighborList,
		Eigen::MatrixXd& derivative)
	{
		//compute the distance to tangent plane of its neighbors
		int N = domain.size();  ///number of the vertex in the same color(domain)
		int V_index;
		Eigen::Vector3d Pos0, pos, pos_normal, neighbor, neighbor_normal, L, L_normal, d, d_tmp, d_tmp2;
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
				int N_half = NN / 2 + 1;
				//if (N_half < 3)
				//	N_half = 3;
				for (std::vector<int>::iterator first = neighbors.begin(); first < neighbors.end(); first++)
				{
					Eigen::Vector3d center(0, 0, 0);

					for (size_t t = 0; t < N_half; t++)
					{
						int j = t;
						if ((first - neighbors.begin()) + t >= NN)
							j = t - NN;

						center += V.row(*(first + j));
					}
					center = center / N_half;
					Eigen::Vector3d move = center - pos;

					d_tmp = timestep * move.dot(pos_normal) * pos_normal;


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

			}
			else
				derivative.row(V_index) = d;
			//cout << "on no.... not manifold..." << endl;
			//continue;
		}
	}
};

