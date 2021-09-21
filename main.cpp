//#include <igl/gaussian_curvature.h>
//#include <igl/principal_curvature.h>
//#include <igl/massmatrix.h>
//#include <igl/invert_diag.h>
//#include <igl/readPLY.h>
//#include <igl/readOFF.h>
//#include <iostream>
////#include <igl/opengl/glfw/Viewer.h>
//#include <igl/jet.h>
//#include <igl/avg_edge_length.h>
//#include <igl/cotmatrix.h>
//#include <igl/invert_diag.h>
//#include <igl/massmatrix.h>
//#include <igl/parula.h>
//#include <igl/per_corner_normals.h>
//#include <igl/per_face_normals.h>
//#include <igl/principal_curvature.h>

#include <igl/read_triangle_mesh.h>
#include "personal_path.h"
#ifdef __APPLE__
#include <igl/OpenGL_convenience.h>
#endif
#include <Eigen/Dense>
#include <igl/readOBJ.h>
#include <set>
//#include "MeshAnalysis.h"
#include "flow.h"
#include "GCfilter.h"
#include "GaussianFlow.h"
#include "HWLaplace.h"
#include "LaplaceBeltrami.h"
#include "GGP.h"

#include "lodepng.h"
#include <thread>
#include <chrono>
//#include <igl/opengl/glfw/Viewer.h>
//#include <imgui/imgui.h>
//#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
//#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
//#include <imgui/imgui.h>
#include <direct.h>
#include <iostream>
#include <time.h>

using namespace Eigen;
using namespace std;

int Viewer_i = 0, ItNum;  ////////////number of iteration
Eigen::MatrixXd Viewer_V;
Eigen::MatrixXi Viewer_F;
double zoom;


inline unsigned save_png(const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>& R,
	const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>& G,
	const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>& B,
	const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic>& A,
	const std::string png_file)
{
	assert((R.rows() == G.rows()) && (G.rows() == B.rows()) && (B.rows() == A.rows()));
	assert((R.cols() == G.cols()) && (G.cols() == B.cols()) && (B.cols() == A.cols()));
	int Height = R.rows(), Width = R.cols();
	const int comp = 4;                                  // 4 Channels Red, Green, Blue, Alpha
	const int stride_in_bytes = R.rows()*comp;           // Length of one row in bytes
	std::vector<unsigned char> data(R.size()*comp, 0);     // The image itself;

	for (std::int32_t i = 0; i < Height; ++i)
	{
		for (std::int32_t j = 0; j < Width; ++j)
		{
			data[(j * R.rows() * comp) + (i * comp) + 0] = R(i, R.cols() - 1 - j);
			data[(j * R.rows() * comp) + (i * comp) + 1] = G(i, R.cols() - 1 - j);
			data[(j * R.rows() * comp) + (i * comp) + 2] = B(i, R.cols() - 1 - j);
			data[(j * R.rows() * comp) + (i * comp) + 3] = A(i, R.cols() - 1 - j);
		}
	}
	//std::vector<std::uint8_t> PngBuffer(data.size());

	//for (std::int32_t I = 0; I < Height; ++I)
	//{
	//	for (std::int32_t J = 0; J < Width; ++J)
	//	{
	//		std::size_t OldPos = (Height - I - 1) * (Width * 4) + 4 * J;
	//		std::size_t NewPos = I * (Width * 4) + 4 * J;
	//		PngBuffer[NewPos + 0] = data[OldPos + 2]; //B is offset 2
	//		PngBuffer[NewPos + 1] = data[OldPos + 1]; //G is offset 1
	//		PngBuffer[NewPos + 2] = data[OldPos + 0]; //R is offset 0
	//		PngBuffer[NewPos + 3] = data[OldPos + 3]; //A is offset 3
	//	}
	//}

	std::vector<std::uint8_t> ImageBuffer;
	unsigned error = lodepng::encode(ImageBuffer, data, Height, Width);
	lodepng::save_file(ImageBuffer, png_file);

	return error;
}



//// This function is called every time a keyboard button is pressed
//bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
//{
//
//	//std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
//	if (key == '1') ////display the Gaussian Curvature
//	{
//		Viewer_i++;
//		if (Viewer_i > ItNum)
//			Viewer_i--;
//		std::cout << "current frame: " << Viewer_i << std::endl;
//		//cout << ItNum << endl;
//		// Clear should be called before drawing the mesh
//		viewer.data().clear();
//		// Draw_mesh creates or updates the vertices and faces of the displayed mesh.
//		// If a mesh is already displayed, draw_mesh returns an error if the given V and
//		// F have size different than the current ones
//		//viewer.core().background_color.setConstant(1);
//		Viewer_V = Vlist[Viewer_i];
//		//Viewer_F = Flist[Viewer_i];
//		viewer.data().set_mesh(Viewer_V, Viewer_F);
//		viewer.data().compute_normals();
//		//viewer.data().set_colors(Viewer_C);
//
//		viewer.core().align_camera_center(Viewer_V, Viewer_F);
//
//	}
//	else if (key == '2') ////display the Gaussian Curvature
//	{
//		Viewer_i--;
//		if (Viewer_i < 0)
//			Viewer_i++;
//		std::cout << "current frame: " << Viewer_i << std::endl;
//		// Clear should be called before drawing the mesh
//		viewer.data().clear();
//		// Draw_mesh creates or updates the vertices and faces of the displayed mesh.
//		// If a mesh is already displayed, draw_mesh returns an error if the given V and
//		// F have size different than the current ones
//		viewer.core().background_color.setConstant(1);
//		Viewer_V = Vlist[Viewer_i];
//		//Viewer_F = Flist[Viewer_i];
//		viewer.data().set_mesh(Viewer_V, Viewer_F);
//		viewer.data().compute_normals();
//		//viewer.data().set_colors(Viewer_C);
//		viewer.core().align_camera_center(Viewer_V, Viewer_F);
//
//	}
//	else if (key == '3') ////display the Gaussian Curvature
//	{
//		// Allocate temporary buffers
//		//Viewer_V = Vlist[Viewer_i];
//		//Viewer_F = Flist[Viewer_i];
//		Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R(1280, 800);
//		Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> G(1280, 800);
//		Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> B(1280, 800);
//		Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> A(1280, 800);
//
//		// Draw the scene in the buffers
//		viewer.core().draw_buffer(
//			viewer.data(), false, R, G, B, A);
//
//		//////// Save it to a PNG
//		std::string filename = OUTPUT_PATH "/screenshot" + std::to_string(Viewer_i) + ".png";
//		unsigned err = save_png(R, G, B, A, filename);
//
//
//		//if there's an error, display it
//		if (err) std::cout << "encoder error: " << lodepng_error_text(err) << std::endl;
//		else
//			std::cout << "screenshot" + std::to_string(Viewer_i) + ".png saved successfully!" << std::endl;
//	}
//	else if (key == ' ')
//	{
//		Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R(1280, 800);
//		Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> G(1280, 800);
//		Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> B(1280, 800);
//		Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> A(1280, 800);
//		for (size_t i = 0; i < Vlist.size(); i++)
//		{
//			std::cout << "current frame: " << i << std::endl;
//			viewer.data().clear();
//			//viewer.core().clear_framebuffers();
//			Viewer_V = Vlist[i];
//			viewer.core().background_color.setConstant(1);
//			viewer.data().set_mesh(Viewer_V, Viewer_F);
//			viewer.data().compute_normals();
//			//viewer.data().set_colors(Viewer_C);
//			//viewer.core().model_zoom = zoom;
//			viewer.core().align_camera_center(Viewer_V, Viewer_F);
//
//			//std::this_thread::sleep_for(std::chrono::milliseconds(1000));
//			viewer.core().draw_buffer(
//				viewer.data(), false, R, G, B, A);
//			//////// Save it to a PNG
//			std::string filename = OUTPUT_PATH "/auto_screenshot" + std::to_string(i) + ".png";
//			unsigned err = save_png(R, G, B, A, filename);
//
//
//			//if there's an error, display it
//			if (err) std::cout << "encoder error: " << lodepng_error_text(err) << std::endl;
//			else
//				std::cout << "auto_screenshot" + std::to_string(i) + ".png saved successfully!" << std::endl;
//		}
//	}
//	return false;
//}

#define HLO
int main(int argc, char * argv[])
//////////////////////////////////////////////////////////////
////////////////////////main part
{
	std::string filename = USER_DEFINED_PATH OBJNAME;


	if (argc != 3)
	{
		cout << "Usage: "<<endl<<"Surface_Smoothing.exe filename ItNum" << endl;
		//cout << "Filter Type: " << endl;
		//cout << "1: GGP" << endl;
		//cout << "11: Half-kernel Laplace by COS approximtes to -1" << endl;
		//cout << "111: Half-kernel Laplace by COS no projection" <<endl;
		//cout << "12: Half-kernel Laplace by HALF Neighboorhood" << endl;
		//cout << "13: Half-kernel Laplace by cut and mindistance" << endl;
		//cout << "3: Projective Gaussian Curvature Filter" << endl;
		//cout << "4: PointWise Laplace" << endl;
		//cout << "5: Cotangent Laplace" << endl;
		//cout << "6: Gaussian Flow" << endl;
		//cout << "7: Uniform Laplace" << endl;
		return -1;
		
	}

	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	Eigen::MatrixXd Vn2, VNormal;
	string objname(argv[1]);
	//string objname("Noisy(1).obj");
	filename = /*USER_DEFINED_PATH +*/ objname;

	igl::read_triangle_mesh(filename, V, F);
	if (F.cols() == 0 || V.rows() == 0)
	{
		cout << "Empty Data!" << endl;
		return -1;
	}


	//igl::readOBJ(argv[1], V, F);
	//igl::readOFF(argv[1],V,F);
	ItNum = atoi(argv[2]);
	//const double timestep = atof(argv[4]);
	const double timestep = 1;
	cout << V.rows() << " vertices, " << F.rows() << " faces" << endl;

	// add noise
	//noise = Eigen::MatrixXd::Random(V.rows(),3)*0.05;
	//cout<<noise.rows()<<endl;
	//V = V + noise;
	//igl::per_vertex_normals(V,F, Vn2);    

	//Eigen::VectorXd curvatureGaussian;
	//igl::gaussian_curvature(V, F, curvatureGaussian);
	_mkdir(OUTPUT_PATH);
	filename = OUTPUT_PATH "/stats_compare.txt";
	std::ofstream log(filename);

	/************************/
	//MeshAnalysis analyzer(V, F);
	//analyzer.analyze();
	//cout << "Gaussian curvature integral before flow: " << analyzer.getGaussianCurvatureIntegralOverArea() << endl;
	//cout << "Mean curvature energy before flow: " << analyzer.getMeanCurvatureIntegralOverArea() << endl;
	//log << "Gaussian curvature integral before flow: " << analyzer.getGaussianCurvatureIntegralOverArea() << endl;
	//log << "Mean curvature energy before flow: " << analyzer.getMeanCurvatureIntegralOverArea() << endl;

	string filtertype;
	filtertype = "Half Window Laplace";
	clock_t time1 = clock();
	flow<HWLaplace_4Cuts>(V, F, ItNum, timestep);

	//switch (atoi(argv[1]))
	//{
	//case 1:
	//	flow<GGPFilter>(V, F, ItNum, timestep);
	//	filtertype = "GGP Filter";
	//	break;
	//case 11:
	//	flow<HWLaplaceCOS>(V, F, ItNum, timestep);
	//	filtertype = "Half Window Laplace by COS";
	//	break;
	//case 111:
	//	flow<HWLaplaceCOS_NP>(V, F, ItNum, timestep);
	//	filtertype = "Half Window Laplace by COS no projection";
	//	break;
	//case 12:
	//	flow<HWLaplaceHN>(V, F, ItNum, timestep);
	//	filtertype = "Half Window Laplace by HALF Neighboorhood";
	//	break;
	//case 13:
		//flow<HWLaplace_4Cuts>(V, F, ItNum, timestep);
	//	filtertype = "Half Window Laplace";
	//	break;
	//case 3:
	//	flow<ProjectiveGCFilter>(V, F, ItNum, timestep);
	//	filtertype = "Projective GCFilter";
	//	break;
	//case 4:
	//	flow<PointWiseLaplace>(V, F, ItNum, timestep);
	//	filtertype = "PointWise Laplace";
	//	break;
	//case 5:
	//	flow<LaplaceBeltrami_v1>(V, F, ItNum, timestep);
	//	filtertype = "cotangent Laplace";
	//	break;
	//case 6:
	//	flow<GaussianFlow>(V, F, ItNum, timestep);
	//	filtertype = "Gaussian Flow";
	//	break;
	//case 7:
	//	flow<LaplaceBeltrami_v2>(V, F, ItNum, timestep);
	//	filtertype = "uniform Laplace";
	//	break;
	//default:
	//	flow<GGPFilter>(V, F, ItNum, timestep);
	//	filtertype = "GGP Filter";
	//	break;
	//}
	clock_t time2 = clock();

	cout << "denoising time in total(ms): " << time2 - time1 << endl;
	filename = OUTPUT_PATH + objname + "_denoised.obj";
	igl::writeOBJ(filename, V, F);

	//analyzer.analyze();
	//cout << "Gaussian curvature integral after flow: " << analyzer.getGaussianCurvatureIntegralOverArea() << endl;
	//cout << "Mean curvature energy after flow: " << analyzer.getMeanCurvatureIntegralOverArea() << endl;
	//log << "Gaussian curvature integral before flow: " << analyzer.getGaussianCurvatureIntegralOverArea() << endl;
	//log << "Mean curvature energy before flow: " << analyzer.getMeanCurvatureIntegralOverArea() << endl;

	//log << "Complete!" << endl;
	log << "denoising time in total(ms): " << time2 - time1 << endl;

	log.close();

	//// Use original normals as pseudo-co22lors
	//MatrixXd N;
	//igl::per_vertex_normals(V, F, N);
	//MatrixXd C = N.rowwise().normalized().array()*0.5 + 0.5;
	//cout << "#####################################" << endl;
	//cout << "Usage: " << endl;
	//cout << "1: " << "Next Frame;" << endl;
	//cout << "2: " << "Previous Frame" << endl;
	//cout << "3: " << "Screenshot on current frame" << endl;
	//cout << "Space: " << "Automatically display and screenshot all frames" << endl;
	//cout << "#####################################" << endl;

	//igl::opengl::glfw::Viewer viewer;


	//Viewer_V = Vlist[0];
	//Viewer_F = F;
	//viewer.callback_key_down = &key_down;


	//// Attach a menu plugin

	//igl::opengl::glfw::imgui::ImGuiMenu menu;
	//viewer.plugins.push_back(&menu);


	//// Customize the menu
	//float floatVariable = 0.1f; // Shared between two menus

	//							// Add content to the default menu window
	//menu.callback_draw_viewer_menu = [&]()
	//{
	//	// Draw parent menu content
	//	menu.draw_viewer_menu();

	//	// Add new group
	//	if (ImGui::CollapsingHeader("Set Cursor Position", ImGuiTreeNodeFlags_DefaultOpen))
	//	{
	//		// Expose variable directly ...
	//		ImGui::InputFloat("float", &floatVariable, 0, 0, 3);

	//		// ... or using a custom callback
	//		static bool boolVariable = true;
	//		if (ImGui::Checkbox("bool", &boolVariable))
	//		{
	//			// do something
	//			std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
	//		}

	//		// Expose an enumeration type
	//		enum Orientation { Up = 0, Down, Left, Right };
	//		static Orientation dir = Up;
	//		ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");

	//		// We can also use a std::vector<std::string> defined dynamically
	//		static int num_choices = 3;
	//		static std::vector<std::string> choices;
	//		static int idx_choice = 0;
	//		if (ImGui::InputInt("Num letters", &num_choices))
	//		{
	//			num_choices = std::max(1, std::min(26, num_choices));
	//		}
	//		if (num_choices != (int)choices.size())
	//		{
	//			choices.resize(num_choices);
	//			for (int i = 0; i < num_choices; ++i)
	//				choices[i] = std::string(1, 'A' + i);
	//			if (idx_choice >= num_choices)
	//				idx_choice = num_choices - 1;
	//		}
	//		ImGui::Combo("Letter", &idx_choice, choices);

	//		// Add a button
	//		if (ImGui::Button("Print Hello", ImVec2(-1, 0)))
	//		{
	//			std::cout << "Hello\n";
	//		}
	//	}
	//};

	//// Draw additional windows
	//menu.callback_draw_custom_window = [&]()
	//{
	//	// Define next window position + size
	//	ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
	//	ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
	//	ImGui::Begin(
	//		"Mesh Info", nullptr,
	//		ImGuiWindowFlags_NoSavedSettings
	//	);

	//	// Expose the same variable directly ...
	//	ImGui::PushItemWidth(-80);
	//	ImGui::DragFloat("float", &floatVariable, 0.0, 0.0, 3.0);
	//	ImGui::PopItemWidth();

	//	//static std::string str = "bunny";
	//	ImGui::InputText("Name", objname);
	//	ImGui::InputText("Vertices", to_string(V.rows()));
	//	ImGui::InputText("Faces", to_string(F.rows()));
	//	ImGui::InputText("Filter", filtertype);


	//	ImGui::End();
	//};


	//viewer.core().background_color.setConstant(1);
	////viewer.data().set_colors(C);
	//viewer.data().set_mesh(Viewer_V, Viewer_F);

	////viewer.core().align_camera_center(Vlist[0]);
	////zoom= viewer.core().model_zoom;
	//viewer.launch();

	return 0;

}


/////////////////testing part: for paper figure
///////////////////////////////////////////////////////////////////
////
//{

//	Eigen::MatrixXd V;
//	Eigen::MatrixXi F;
//
//	Eigen::MatrixXd Vn2, VNormal;
//	//string objname(argv[2]);	
//	//string objname = "F:\\meshfilter\\meshfilter\\x64\\Debug\\Original.obj";
//	//string objname = "C:\\Users\\Qiu\\Desktop\\test3dfile\\PWL\\0.obj";
//	//filename = USER_DEFINED_PATH+objname;
//	//string filename = objname;
//	//igl::read_triangle_mesh(filename, V, F);
//
//
//	Eigen::MatrixXd V0, V1, V2, V3, V4, V5, V6;
//	Eigen::MatrixXi F0, F1, F2, F3, F4, F5, F6;
//
//
//	string objname0 = "I:\\libigl\\output\\20181130\\Noisy.obj";
//	igl::read_triangle_mesh(objname0, V0, F0);
//
//	VectorXd K;
//	MatrixXd C0;
//	const int n_Vertices = V0.rows();
//	Eigen::MatrixXd VerticesNormal(n_Vertices, 3);
//	K.resize(V0.rows(), 1);
//	igl::per_vertex_normals(V0, F0, VerticesNormal);
//
//	int i = 9;
//
//	//for (int i = 0; i < 20; i++)
//	//{
//#ifdef LBO
//	string objname6 = "I:\\libigl\\output\\20181130\\4\\" + to_string(i) + ".obj";
//
//#endif // LBO
//
//#ifdef HLO
//	string objname6 = "I:\\libigl\\output\\20181130\\11\\" + to_string(i) + ".obj";
//#endif // HLO
//
//
//	igl::read_triangle_mesh(objname6, V6, F6);
//	cout << "K .rows():" << K.rows() << endl;
//	cout << "V0 .rows():" << V0.rows() << endl;
//	cout << "F0 .rows():" << F0.rows() << endl;
//	cout << "V4.rows():" << V6.rows() << endl;
//	cout << "F4.rows():" << F6.rows() << endl;
//	double minnumber = 0.0, maxnumber = 0.0;
//	double sum = 0.0;
//	for (int i = 0; i < n_Vertices; i++)
//	{
//		//cout << "i: " << i << endl;
//		Vector3d tempDistance = Vector3d::Zero();
//		Vector3d tempVerticesNormal = Vector3d::Zero();
//		Vector3d tempDistanceNormalize = Vector3d::Zero();
//		tempVerticesNormal = VerticesNormal.row(i);
//		tempVerticesNormal.normalize();
//		tempDistance = V0.row(i) - V6.row(i);
//		tempDistanceNormalize = tempDistance;
//		tempDistanceNormalize.normalize();
//		double D = 0.0;
//
//		if ((tempDistanceNormalize.dot(tempVerticesNormal)) > 0.0)
//		{
//			D = tempDistance.norm();
//		}
//		else
//		{
//			D = -1.0*tempDistance.norm();
//		}
//		sum = sum + abs(D);
//		K(i, 0) = D;
//		if (maxnumber < D)
//		{
//			maxnumber = D;
//		}
//		if (minnumber > D)
//		{
//			minnumber = D;
//		}
//
//	}
//	cout << "sum:" << sum << endl;
//	cout << "maxnumber:" << maxnumber << endl;
//	cout << "minnumber:" << minnumber << endl;
//
//
//#ifdef LBO
//	igl::jet(K, true, C0);
//#endif // LBO
//
//
//#ifdef HLO
//	fstream dataFile;
//	string output = "I:\\libigl\\output\\20181130\\11\\" + to_string(i) + ".txt";
//	dataFile.open(output);
//	//for (int i = 0; i < n_Vertices; i++)
//	//{
//		//dataFile << V0.row(i).x()<<" "<< V0.row(i).y()<<" "<< V0.row(i).z()<<" "<<C0.row(i).x()<<" "<< C0.row(i).y()<<" "<< C0.row(i).z()<< "\n";
//		//dataFile << V0.row(i).x() << " " << V0.row(i).y() << " " << V0.row(i).z() << " " << K.row(i) << "\n";
//		//dataFile << K.row(i) << "\n";
//	dataFile >> maxnumber;
//	dataFile >> minnumber;
//	//}
//
//	dataFile.close();
//	//}
//	igl::jet(K, minnumber, maxnumber, C0);
//#endif // 
//
//
//
//	// ¹Ø±ÕÎÄµµ
//
//	// Attach a menu plugin
//
//
//	igl::opengl::glfw::Viewer viewer0;
//	igl::opengl::glfw::imgui::ImGuiMenu menu0;
//	viewer0.plugins.push_back(&menu0);
//
//	// Customize the menu
//	float floatVariable0 = 0.1f; // Shared between two menus
//
//								// Add content to the default menu window
//	menu0.callback_draw_viewer_menu = [&]()
//	{
//		// Draw parent menu content
//		menu0.draw_viewer_menu();
//
//		// Add new group
//		if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
//		{
//			// Expose variable directly ...
//			ImGui::InputFloat("float", &floatVariable0, 0, 0, 3);
//
//			// ... or using a custom callback
//			static bool boolVariable = true;
//			if (ImGui::Checkbox("bool", &boolVariable))
//			{
//				// do something
//				std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
//			}
//
//			// Expose an enumeration type
//			enum Orientation { Up = 0, Down, Left, Right };
//			static Orientation dir = Up;
//			ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");
//
//			// We can also use a std::vector<std::string> defined dynamically
//			static int num_choices = 3;
//			static std::vector<std::string> choices;
//			static int idx_choice = 0;
//			if (ImGui::InputInt("Num letters", &num_choices))
//			{
//				num_choices = std::max(1, std::min(26, num_choices));
//			}
//			if (num_choices != (int)choices.size())
//			{
//				choices.resize(num_choices);
//				for (int i = 0; i < num_choices; ++i)
//					choices[i] = std::string(1, 'A' + i);
//				if (idx_choice >= num_choices)
//					idx_choice = num_choices - 1;
//			}
//			ImGui::Combo("Letter", &idx_choice, choices);
//
//			// Add a button
//			if (ImGui::Button("Print Hello", ImVec2(-1, 0)))
//			{
//				std::cout << "Hello\n";
//			}
//		}
//	};
//
//	// Draw additional windows
//	menu0.callback_draw_custom_window = [&]()
//	{
//		// Define next window position + size
//		ImGui::SetNextWindowPos(ImVec2(180.f * menu0.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
//		ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
//		ImGui::Begin(
//			"Mesh Info", nullptr,
//			ImGuiWindowFlags_NoSavedSettings
//		);
//
//		// Expose the same variable directly ...
//		ImGui::PushItemWidth(-80);
//		ImGui::DragFloat("float", &floatVariable0, 0.0, 0.0, 3.0);
//		ImGui::PopItemWidth();
//
//		//static std::string str = "bunny";
//		ImGui::InputText("Name", objname0);
//		ImGui::InputText("Vertices", to_string(V.rows()));
//		ImGui::InputText("Faces", to_string(F.rows()));
//
//		ImGui::End();
//	};
//
//	//Eigen::MatrixXd Vharlf;
//	//Vharlf.resize(V0.rows(), 3);
//	//Eigen::MatrixXi Fharlf;
//	//Fharlf.resize(F0.rows(), 3);
//	//const int NumbersVHarlf = int(Vharlf.rows() / 2);
//	//const int NumbersFHarlf = int(Fharlf.rows() / 2);
//	//Vector3d zerotv = Vector3d::Zero();
//	//Vector3i zerotf = Vector3i::Zero();
//	//for (int i = 0; i < V0.rows(); i++)
//	//{
//	//	if (i < NumbersVHarlf)
//	//	{
//	//		Vharlf.row(i)[0] = V0.row(i).x();
//	//		Vharlf.row(i)[1] = V0.row(i).y();
//	//		Vharlf.row(i)[2] = V0.row(i).z();
//	//	}
//	//	else
//	//	{
//	//		Vharlf.row(i)[0] = 0.0;
//	//		Vharlf.row(i)[1] = 0.0;
//	//		Vharlf.row(i)[2] = 0.0;
//	//	}
//
//	//}
//
//	//for (int i = 0; i < Fharlf.rows(); i++)
//	//{
//	//	if (i < NumbersFHarlf)
//	//	{
//	//		Fharlf.row(i) = F0.row(i);
//	//	}
//	//	else
//	//	{
//	//		Fharlf.row(i) = zerotf;
//	//	}
//
//	//}
//
//	viewer0.core().background_color.setConstant(1);
//	viewer0.data().set_mesh(V0, F0);
//	igl::jet(K, -100, 100, C0);
//	viewer0.data().set_colors(C0);
//
//	//// Allocate temporary buffers for 1280x800 image
//	//Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R(1280, 800);
//	//Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> G(1280, 800);
//	//Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> B(1280, 800);
//	//Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> A(1280, 800);
//
//	//// Draw the scene in the buffers
//	//viewer0.core().draw_buffer(viewer0.data(), false, R, G, B, A);
//
//	//// Save it to a PNG
//	//igl::png::writePNG(R, G, B, A, "out1.png");
//
//	//viewer.core().align_camera_center(Vlist[0]);
//	//zoom= viewer.core().model_zoom;
//	viewer0.launch();
//
//
//
//	return 0;
//
//}




