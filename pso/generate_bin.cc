/**
 * @file GridMatrixEngine.cpp
 * @brief Integrated pipeline to generate, stream, and display an abstract grid matrix.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <string>

#pragma pack(push, 1)
struct ContinuousPoint {
    double y; // Row index tracker (0 to N)
    double x; // Column index tracker (0 to M)
};
#pragma pack(pop)

// ============================================================================
// GENERATION ENGINE FUNCTION
// ============================================================================
bool generateGridMatrix(const std::string& filename, size_t rows, size_t cols) {
    std::ofstream binFile(filename, std::ios::binary | std::ios::trunc);
    if (!binFile.is_open()) {
        std::cerr << ">> [ERROR] Failed to open write stream target: " << filename << "\n";
        return false;
    }

    std::cout << ">> [GENERATOR] Building abstract " << rows << " x " << cols
              << " matrix (" << (rows * cols) << " elements)...\n";

    for (size_t r = 0; r < rows; ++r) {
        for (size_t c = 0; c < cols; ++c) {
            ContinuousPoint point = { static_cast<double>(r), static_cast<double>(c) };
            binFile.write(reinterpret_cast<const char*>(&point), sizeof(ContinuousPoint));
        }
    }

    binFile.close();
    std::cout << ">> [GENERATOR] Successfully flushed grid stream to disks.\n\n";
    return true;
}

// ============================================================================
// DISPLAY SCANNER FUNCTION
// ============================================================================
void displayGridMatrix(const std::string& filename, size_t limit = 15) {
    std::ifstream binFile(filename, std::ios::binary | std::ios::ate);
    if (!binFile.is_open()) {
        std::cerr << ">> [ERROR] Target binary asset unreadable.\n";
        return;
    }

    std::streamsize fileSize = binFile.tellg();
    binFile.seekg(0, std::ios::beg);

    size_t totalElements = fileSize / sizeof(ContinuousPoint);
    std::cout << "=== BINARY STORAGE MAP DISCOVERY ===\n";
    std::cout << "Located File Target  : " << filename << "\n";
    std::cout << "Total Bytes Stored   : " << fileSize << " bytes\n";
    std::cout << "Total Array Nodes    : " << totalElements << " cells\n";
    std::cout << "====================================\n";

    // Restrict output loop to prevent terminal scrolling overload
    size_t scanWindow = std::min(totalElements, limit);
    std::vector<ContinuousPoint> displayBuffer(scanWindow);

    binFile.read(reinterpret_cast<char*>(displayBuffer.data()), scanWindow * sizeof(ContinuousPoint));
    binFile.close();

    for (size_t i = 0; i < scanWindow; ++i) {
        std::cout << "Node Offset [" << std::setw(3) << i << "] -> "
                  << "Row Index (Y): " << std::setw(4) << displayBuffer[i].y << " | "
                  << "Col Index (X): " << std::setw(4) << displayBuffer[i].x << "\n";
    }

    if (totalElements > limit) {
        std::cout << "... [Truncated payload output window view for display protection] ...\n";
    }
}

// ============================================================================
// EXECUTION PIPELINE ENTRYPOINT
// ============================================================================
int main() {
    const std::string targetFile = "pso/feasible.bin";

    // Set abstract matrix resolution constraints (Independent of maps)
    const size_t MATRIX_ROWS = 1500; // N range
    const size_t MATRIX_COLS = 1500; // M range

    // 1. Fire generation loop to populate file
    if (generateGridMatrix(targetFile, MATRIX_ROWS, MATRIX_COLS)) {

        // 2. Fire reader parsing to display the top records
        displayGridMatrix(targetFile, 12);
    }

    return 0;
}
