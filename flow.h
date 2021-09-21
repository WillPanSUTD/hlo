#pragma once

// This file contains methods for generalized flows on triangle meshes.

#define USE_DOMAIN_DECOMPOSITION

#include "MeshAnalysis.h"

#include <vector>
#include <set>
#include <iostream>
#include <fstream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/per_vertex_normals.h>
#include "LaplaceBeltrami.h"

////for projective GCfilter
#define USE_TRI_PROJ
#define USE_LINE_PROJ


//#define USE_BARYCENTRIC_CHECK

////Energy definition
#define USE_ENERGY_VERTEX
//#define USE_ENERGY_DISPOSITION



using namespace Eigen;
using namespace std;
//#define OUTPUT_PATH "/output"




vector<Eigen::MatrixXd> Vlist;
//vector<Eigen::MatrixXi> Flist;
Eigen::MatrixXd Viewer_C;
//Greedily solves the vertex coloring problem
//returns number of used colors
int ColorIt(Eigen::MatrixXd &V, Eigen::VectorXi &color, vector<vector<int> > & NeighborList)
{
	int maxColor = -1;

	for (int i = 0; i < V.rows(); ++i)
	{
		//color vertex i

		//retrieve all adjacent colors
		std::set<int> adjacentColors;
		for (auto n : NeighborList[i])
		{
			if (color(n) >= 0)
				adjacentColors.insert(color(n));
		}

		//choose the first unused color
		int c = 0;
		while (adjacentColors.find(c) != adjacentColors.end())
			++c;

		color(i) = c;

		if (c > maxColor)
			maxColor = c;
	}

	return maxColor + 1;
}


int ColorIt(Eigen::MatrixXd &V, Eigen::VectorXi &color, vector<set<int> > & NeighborList)
{
	int maxColor = -1;

	for (int i = 0; i < V.rows(); ++i)
	{
		//color vertex i

		//retrieve all adjacent colors
		std::set<int> adjacentColors;
		for (auto n : NeighborList[i])
		{
			if (color(n) >= 0)
				adjacentColors.insert(color(n));
		}

		//choose the first unused color
		int c = 0;
		while (adjacentColors.find(c) != adjacentColors.end())
			++c;

		color(i) = c;

		if (c > maxColor)
			maxColor = c;
	}

	return maxColor + 1;
}


//Checks for a valid coloring
//Returns true iff all neighbors of a vertex v have a different color than v
bool checkColoring(const vector<set<int> >& neighborList, const Eigen::VectorXi &  colors)
{
	for (int i = 0; i < neighborList.size(); ++i)
	{
		///if a neighbour of vertex i has the same color with v_i, return false
		for (auto n : neighborList.at(i))
		{
			if (colors(i) == colors(n))
				return false;
		}
	}
	return true;
}

bool checkColoring(const vector<vector<int> >& neighborList, const Eigen::VectorXi &  colors)
{
	for (int i = 0; i < neighborList.size(); ++i)
	{
		///if a neighbour of vertex i has the same color with v_i, return false
		for (auto n : neighborList.at(i))
		{
			if (colors(i) == colors(n))
				return false;
		}
	}
	return true;
}




template <typename DerivedV, typename DerivedF>
bool writeOBJWithColor(
	const std::string str,
	const Eigen::PlainObjectBase<DerivedV>& V,
	const Eigen::PlainObjectBase<DerivedF>& F,
	const Eigen::VectorXi &  C)
{
	FILE * obj_file = fopen(str.c_str(), "w");
	if (NULL == obj_file)
	{
		printf("IOError: %s could not be opened for writing...", str.c_str());
		return false;
	}
	// Loop over V
	for (int i = 0; i < (int)V.rows(); i++)
	{
		int r = (C(i) % 2) * 255;
		int g = ((C(i) / 2) % 2) * 255;
		int b = ((C(i) / 4) % 2) * 255;
		fprintf(obj_file, "v %0.15g %0.15g %0.15g %i %i %i\n",
			V(i, 0),
			V(i, 1),
			V(i, 2),
			r, g, b
		);
	}

	// loop over F
	for (int i = 0; i < (int)F.rows(); ++i)
	{
		fprintf(obj_file, "f");
		for (int j = 0; j < (int)F.cols(); ++j)
		{
			// OBJ is 1-indexed
			fprintf(obj_file, " %u", F(i, j) + 1);
		}
		fprintf(obj_file, "\n");
	}

	return true;
}


void VertexNormal(Eigen::MatrixXd & V, Eigen::MatrixXi & F, Eigen::MatrixXd & VNormal)
{
	int N = V.rows();  ///nunber of vertices
	int M = F.rows();	////number of faces

	VNormal = Eigen::MatrixXd::Constant(N, 3, 0);
	Eigen::Vector3d norm_tmp1, norm_tmp2;

	for (int i = 0; i < M; ++i)
	{
		auto a = V.row(F(i, 0));
		auto b = V.row(F(i, 1));
		auto c = V.row(F(i, 2));

		norm_tmp1 = V.row(F(i, 0)) - V.row(F(i, 1));
		norm_tmp2 = V.row(F(i, 1)) - V.row(F(i, 2));
		norm_tmp1 = norm_tmp1.cross(norm_tmp2);
		//norm_tmp1.normalize();

		for (int j = 0; j < 3; ++j)
		{
			VNormal.row(F(i, j)) += norm_tmp1;
		}
	}

	for (int i = 0; i < N; ++i)
	{
		VNormal.row(i).normalize();
	}
}

template <typename TFlow>
void flow(Eigen::MatrixXd & V, Eigen::MatrixXi & F, int ItNum, double Tstep)
{
	Eigen::MatrixXd V0 = V;
	const int nVertices = V.rows();
	const int nFaces = F.rows();

	Eigen::VectorXi degree(nVertices);  ///degree for each vertex
	Eigen::VectorXi color(nVertices);
	Eigen::MatrixXd VNormal(nVertices, 3);

	//vertex normal for tangent plane
	//Eigen::MatrixXd VNormal = Eigen::MatrixXd::Constant(N,3,0);
	//// Constant(rows,cols,value) sets all coefficients to value
	VNormal = Eigen::MatrixXd::Constant(nVertices, 3, 0);

	//Eigen::MatrixXd K_before(N, 1), K_after(N, 1);

	//vertex vertex neighborlist
	///set.insert can avoid duplicate elements
	vector<vector<int> > NeighborList(nVertices);

	//init
	for (int i = 0; i < nVertices; ++i) degree(i) = 0;


	//vertex degree, vertex-vertex neighbor list, vertex normal
	for (int i = 0; i < nFaces; ++i)
	{
		std::vector<int>::iterator it1 = find(NeighborList[F(i, 0)].begin(),
			NeighborList[F(i, 0)].end(), F(i, 1));
		std::vector<int>::iterator it2 = find(NeighborList[F(i, 0)].begin(),
			NeighborList[F(i, 0)].end(), F(i, 2));
		if (it1 == NeighborList[F(i, 0)].end() && it2 == NeighborList[F(i, 0)].end())
		{
			NeighborList[F(i, 0)].push_back(F(i, 1));
			NeighborList[F(i, 0)].push_back(F(i, 2));
		}
		else if (it1 != NeighborList[F(i, 0)].end() && it2 == NeighborList[F(i, 0)].end())
		{
			NeighborList[F(i, 0)].insert(it1, F(i, 2));
		}
		else if (it1 == NeighborList[F(i, 0)].end() && it2 != NeighborList[F(i, 0)].end())
		{
			if (it2 - NeighborList[F(i, 0)].begin() == 0)
				it2 = NeighborList[F(i, 0)].begin();
			else
				it2--;
			NeighborList[F(i, 0)].insert(it2, F(i, 1));
		}


		it1 = find(NeighborList[F(i, 1)].begin(),
			NeighborList[F(i, 1)].end(), F(i, 2));
		it2 = find(NeighborList[F(i, 1)].begin(),
			NeighborList[F(i, 1)].end(), F(i, 0));
		if (it1 == NeighborList[F(i, 1)].end() && it2 == NeighborList[F(i, 1)].end())
		{
			NeighborList[F(i, 1)].push_back(F(i, 2));
			NeighborList[F(i, 1)].push_back(F(i, 0));
		}
		else if (it1 != NeighborList[F(i, 1)].end() && it2 == NeighborList[F(i, 1)].end())
		{
			NeighborList[F(i, 1)].insert(it1, F(i, 0));
		}
		else if (it1 == NeighborList[F(i, 1)].end() && it2 != NeighborList[F(i, 1)].end())
		{
			if (it2 - NeighborList[F(i, 1)].begin() == 0)
				it2 = NeighborList[F(i, 1)].begin();
			else
				it2--;
			NeighborList[F(i, 1)].insert(it2, F(i, 2));
		}


		it1 = find(NeighborList[F(i, 2)].begin(),
			NeighborList[F(i, 2)].end(), F(i, 0));
		it2 = find(NeighborList[F(i, 2)].begin(),
			NeighborList[F(i, 2)].end(), F(i, 1));
		if (it1 == NeighborList[F(i, 2)].end() && it2 == NeighborList[F(i, 2)].end())
		{
			NeighborList[F(i, 2)].push_back(F(i, 0));
			NeighborList[F(i, 2)].push_back(F(i, 1));
		}
		else if (it1 != NeighborList[F(i, 2)].end() && it2 == NeighborList[F(i, 2)].end())
		{
			NeighborList[F(i, 2)].insert(it1, F(i, 1));
		}
		else if (it1 == NeighborList[F(i, 2)].end() && it2 != NeighborList[F(i, 2)].end())
		{
			if (it2 - NeighborList[F(i, 2)].begin() == 0)
				it2 = NeighborList[F(i, 2)].begin();
			else
				it2--;
			NeighborList[F(i, 2)].insert(it2, F(i, 0));
		}






		//NeighborList[F(i, 0)].insert(F(i, 1));
		//NeighborList[F(i, 0)].insert(F(i, 2));

		//NeighborList[F(i, 1)].insert(F(i, 0));
		//NeighborList[F(i, 1)].insert(F(i, 2));

		//NeighborList[F(i, 2)].insert(F(i, 0));
		//NeighborList[F(i, 2)].insert(F(i, 1));

		for (int j = 0; j < 3; ++j)
		{
			degree(F(i, j)) += 1;
		}
	}

	/************************ color for decomposition *******************************/
	//init color 
	for (int i = 0; i < nVertices; ++i)
	{
		color(i) = -1; //not colored
	}


	int nLabels = ColorIt(V, color, NeighborList);    ////return the number of colors
	//cout << "Colored vertices using " << nLabels << " colors." << endl;
	//std::string filename = OUTPUT_PATH "coloring.obj";
	//writeOBJWithColor(filename, V, F, color);
	bool validColoring = checkColoring(NeighborList, color);
	if (!validColoring)
	{
		std::cout << "Coloring is invalid." << std::endl;
		return;
	}

	/************************ domain decomposition *******************************/
#ifdef USE_DOMAIN_DECOMPOSITION
	int nDomains = nLabels;
#else
	int nDomains = 1;
#endif
	vector<vector<int> > Domain(nDomains);
	vector<int> Domain_count(nDomains, 0);

	for (int i = 0; i < nVertices; ++i)
	{
#ifdef USE_DOMAIN_DECOMPOSITION
		int domain = color(i);
#else
		int domain = 0;
#endif
		Domain[domain].push_back(i);
		Domain_count[domain]++;
	}
	//filename = OUTPUT_PATH "/stats.txt";
	//std::ofstream stats(filename);
	//stats.imbue(locale("de-DE"));	

	/************************ perform flow *******************************/
	TFlow concreteFlow;
	Eigen::MatrixXd derivative(nVertices, 3);
	MeshAnalysis analyzer(V, F);
	concreteFlow.timestep = Tstep;
	//stats << concreteFlow.name << endl;
	//stats << "Energy profile: Iteration, Gaussian Curvature Energy, Mean Curvature Energy; " << endl;

	//analyzer.analyze();
	//cout << 0 << ", " << analyzer.getGaussianCurvatureIntegralOverArea() << ", " << analyzer.getMeanCurvatureIntegralOverArea() << ";" << endl;
	//stats << 0 << "; " << analyzer.getGaussianCurvatureIntegralOverArea() << "; " << analyzer.getMeanCurvatureIntegralOverArea() << ";" << endl;
	//igl::opengl::glfw::Viewer viewer;
	Vlist.push_back(V);
	//Flist.push_back(F);;
	//Viewer_F = F;
	//Viewer_V = V;
	//viewer.data().set_mesh(Viewer_V, Viewer_F);
	//viewer.launch_shut();

	//viewer.callback_key_down = &key_down;
	// Allocate temporary buffers for 1280x800 image
	//Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R(1280, 800);
	//Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> G(1280, 800);
	//Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> B(1280, 800);
	//Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> A(1280, 800);

	////// Draw the scene in the buffers
	//viewer.core.draw_buffer(viewer.data(), false, R, G, B, A);
	//
	//////// Save it to a PNG
	//save_png(R, G, B, A, "screenshot" + std::to_string(0) + ".png");

	//cout << "Energy profile: Iteration, Gaussian Curvature Energy, Mean Curvature Energy; " << endl;



	//derivative = V;

	for (int i = 0; i < ItNum; ++i)
	{
		igl::per_vertex_normals(V, F, VNormal);
		///Global matrix is not applicable for Domain computation: LBv1,v2
		if (std::is_same_v<TFlow, LaplaceBeltrami_v1> || std::is_same_v<TFlow, LaplaceBeltrami_v2>)
		{
			concreteFlow.calculateDerivative(V0, V, F, Domain[0], VNormal, NeighborList, derivative);

		}

		else {
			VertexNormal(V, F, VNormal);
			//is it same to igl normal function?
			for (int j = 0; j < nDomains; ++j)
			{
				concreteFlow.calculateDerivative(V0, V, F, Domain[j], VNormal, NeighborList, derivative);
				for (int k = 0; k < Domain[j].size(); ++k)
					V.row(Domain[j][k]) += derivative.row(Domain[j][k]);
			}
		}


		Vlist.push_back(V);
		//Flist.push_back(F);

		//Viewer_F = F;
		//Viewer_V = V;
		//viewer.data().set_mesh(Viewer_V, Viewer_F);

		////// Draw the scene in the buffers
		//viewer.core.draw_buffer(viewer.data(), false, R, G, B, A);

		//////// Save it to a PNG
		//save_png(R, G, B, A, "screenshot" + std::to_string(i) + ".png");

		////igl::png::writePNG(R, G, B, A, "screenshot" + std::to_string(i) + ".png");
		//std::string filename = OUTPUT_PATH;
		//igl::writeOBJ(filename + std::to_string(i) + ".obj", V, F);

		//analyzer.analyze();
		//cout << i + 1 << ", " << analyzer.getGaussianCurvatureIntegralOverArea() << ", " << analyzer.getMeanCurvatureIntegralOverArea() << ";" << endl;
		//stats << i + 1 << "; " << analyzer.getGaussianCurvatureIntegralOverArea() << "; " << analyzer.getMeanCurvatureIntegralOverArea() << ";" << endl;
	}



	//stats.close();

	//ofstream result_file;
	//filename = OUTPUT_PATH "/mesh_color.txt";
	//result_file.open(filename);
	//for (int i = 0; i < nVertices; ++i)
	//{
	//	result_file << color(i) << endl;
	//}
	//result_file.close();
	//viewer.launch();

}