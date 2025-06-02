#include <gtest/gtest.h>
#include "ReactionDiffusionEquation/GrayScottRescaledModel/GrayScottRescaledModel.hpp"
#include <iomanip>  // for std::setw, std::setfill
#include <sstream>  // for std::ostringsrtream

class GrayScottRescaledModelTest: public ::testing::Test
{
};

TEST_F(GrayScottRescaledModelTest, TurningPatternTest)
{
	const uint32_t row = 201;
	const uint32_t col = 201;

	// === 初始条件 ===
	Matrix<double> initialA(row, col, 1.0);  // 整个区域A=1.0
	Matrix<double> initialB(row, col, 0.0);  // 整个区域B=0.0

	const double phi = 0.01;
	const int    qx  = 1;
	const int    qy  = 1;

	for(uint32_t i = 0; i < row; ++i) {
		for(uint32_t j = 0; j < col; ++j) {
			double x       = static_cast<double>(i);
			double y       = static_cast<double>(j);
			double perturb = phi * std::cos(2 * 3.1415926 * qx * x / row) * std::cos(2 * 3.1415926 * qy * y / col);
			initialA.setData(i, j, perturb + initialA.getData(i, j));
			initialB.setData(i, j, perturb + initialB.getData(i, j));
		}
	}

	initialA.toGrayScaleImage().savePNG("./build/GrayScottRescaledModelTest_TurningPatternTest_initialA.png");
	initialB.toGrayScaleImage().savePNG("./build/GrayScottRescaledModelTest_TurningPatternTest_initialB.png");

	// === 格子类型 ===
	Matrix<int> gridType(row, col, static_cast<int>(GrayScottRescaledModel::GRID::NONE));
	gridType.setRow(0, static_cast<int>(GrayScottRescaledModel::GRID::NEUMANN));
	gridType.setRow(row - 1, static_cast<int>(GrayScottRescaledModel::GRID::NEUMANN));
	gridType.setCol(0, static_cast<int>(GrayScottRescaledModel::GRID::NEUMANN));
	gridType.setCol(col - 1, static_cast<int>(GrayScottRescaledModel::GRID::NEUMANN));
	gridType.toGrayScaleImage().savePNG("./build/GrayScottRescaledModelTest_TurningPatternTest_gridType.png");

	// === 创建并运行 ===
	GrayScottRescaledModel gsr(2.7339, 2.0139, 10.06, initialA, initialB, gridType);
	gsr.getMatrixA().toGrayScaleImage().savePNG(
		"./build/ReactionDiffusionEquation/GrayScottRescaledModelShader/GrayScottRescaledModelTest_TurningPatternTest_A_00000.png");
	gsr.getMatrixB().toGrayScaleImage().savePNG(
		"./build/ReactionDiffusionEquation/GrayScottRescaledModelShader/GrayScottRescaledModelTest_TurningPatternTest_B_00000.png");
	const int stride = 200;
	for(size_t i = 0; i < 10000; i += stride) {
		gsr.step(stride);

		std::ostringstream filenameA;
		filenameA
			<< "./build/ReactionDiffusionEquation/GrayScottRescaledModelShader/GrayScottRescaledModelTest_TurningPatternTest_A_"
			<< std::setw(5) << std::setfill('0') << i << ".png";
		std::ostringstream filenameB;
		filenameB
			<< "./build/ReactionDiffusionEquation/GrayScottRescaledModelShader/GrayScottRescaledModelTest_TurningPatternTest_B_"
			<< std::setw(5) << std::setfill('0') << i << ".png";
		gsr.getMatrixA().toGrayScaleImage().savePNG(filenameA.str());
		gsr.getMatrixB().toGrayScaleImage().savePNG(filenameB.str());

		if(i % (stride * 10) == 0 || i < 10) {
			std::cerr << gsr.getMatrixA().getMin() << "," << gsr.getMatrixA().getMax() << ","
					  << gsr.getMatrixB().getMin() << "," << gsr.getMatrixB().getMax() << std::endl;
		}
	}
}
