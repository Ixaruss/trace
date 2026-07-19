#define ODLIB_EXPORTS
#include "export_adapt.h" // Compiling inside the same folder
#include <iostream>
#include <vector>

// Global 2D Matrix
std::vector<std::vector<int>> g_matrix;

extern "C" {

ODLIB_API int LoadData(const int* flat_data, int rows, int cols) {
    if (!flat_data || rows <= 0 || cols <= 0) {
        std::cerr << "[DLL Error] Invalid data or dimensions." << std::endl;
        return -1;
    }

    g_matrix.clear();
    g_matrix.resize(rows, std::vector<int>(cols, 0));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            g_matrix[i][j] = flat_data[i * cols + j];
        }
    }

    std::cout << "[DLL] Successfully loaded " << rows << "x" << cols << " matrix." << std::endl;
    return 0;
}

ODLIB_API int DisplayData() {
    if (g_matrix.empty()) {
        std::cout << "[DLL] Matrix is empty. Please load data first." << std::endl;
        return -1;
    }

    std::cout << "\n--- Global 2D Matrix (from DLL memory) ---\n";
    for (const auto& row : g_matrix) {
        for (int val : row) {
            std::cout << val << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "------------------------------------------\n\n";

    return 0;
}

} // extern "C"
