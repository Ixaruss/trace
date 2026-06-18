#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <png.h>

void readFirstResultAndSavePNG(const std::string& binFilename, const std::string& outPngFilename, int maxRadiusCells) {
    // 1. Read metadata for the first radar layout
    std::ifstream binFile(binFilename, std::ios::in | std::ios::binary);
    if (!binFile.is_open()) {
        std::cerr << "[ERROR] Could not open cache file: " << binFilename << "\n";
        return;
    }

    int srcR = 0, srcC = 0, numCovered = 0;
    if (!binFile.read(reinterpret_cast<char*>(&srcR), sizeof(int)) ||
        !binFile.read(reinterpret_cast<char*>(&srcC), sizeof(int)) ||
        !binFile.read(reinterpret_cast<char*>(&numCovered), sizeof(int))) {
        std::cerr << "[ERROR] File empty or corrupted.\n";
        return;
    }

    std::vector<int> coveredRows(numCovered);
    std::vector<int> coveredCols(numCovered);

    if (numCovered > 0) {
        binFile.read(reinterpret_cast<char*>(coveredRows.data()), numCovered * sizeof(int));
        binFile.read(reinterpret_cast<char*>(coveredCols.data()), numCovered * sizeof(int));
    }
    binFile.close();

    std::cout << "[READER] Visualizing Radar Base at Row: " << srcR << ", Col: " << srcC
              << " (" << numCovered << " covered cells)\n";

    // 2. Define the Local Image Dimensions based on the Radar Range bounds
    int rMin = srcR - maxRadiusCells;
    int rMax = srcR + maxRadiusCells;
    int cMin = srcC - maxRadiusCells;
    int cMax = srcC + maxRadiusCells;

    int imgWidth  = (cMax - cMin) + 1;
    int imgHeight = (rMax - rMin) + 1;

    // Allocate continuous block for 24-bit RGB canvas (initialized to 0/Black background)
    std::vector<uint8_t> imageBuffer(imgWidth * imgHeight * 3, 0);

    // 3. Populate visibility indices onto the canvas image
    for (int i = 0; i < numCovered; ++i) {
        int rLocal = coveredRows[i] - rMin;
        int cLocal = coveredCols[i] - cMin;

        if (rLocal >= 0 && rLocal < imgHeight && cLocal >= 0 && cLocal < imgWidth) {
            // Inverting the row height flip matches standard GIS coordinate mapping
            size_t pixelOffset = (size_t)((imgHeight - 1 - rLocal) * imgWidth + cLocal) * 3;
            imageBuffer[pixelOffset + 0] = 50;   // Red
            imageBuffer[pixelOffset + 1] = 205;  // Green (Lime Green)
            imageBuffer[pixelOffset + 2] = 50;   // Blue
        }
    }

    // Draw a quick Red Crosshair (+) over the source station center
    int localSrcR = srcR - rMin;
    int localSrcC = srcC - cMin;
    for (int i = -4; i <= 4; ++i) {
        if (localSrcC + i >= 0 && localSrcC + i < imgWidth) {
            size_t off = (size_t)((imgHeight - 1 - localSrcR) * imgWidth + (localSrcC + i)) * 3;
            imageBuffer[off + 0] = 255; imageBuffer[off + 1] = 0; imageBuffer[off + 2] = 0;
        }
        if (localSrcR + i >= 0 && localSrcR + i < imgHeight) {
            size_t off = (size_t)((imgHeight - 1 - (localSrcR + i)) * imgWidth + localSrcC) * 3;
            imageBuffer[off + 0] = 255; imageBuffer[off + 1] = 0; imageBuffer[off + 2] = 0;
        }
    }

    // 4. Initialize libpng components
    FILE* fp = fopen(outPngFilename.c_str(), "wb");
    if (!fp) {
        std::cerr << "[ERROR] Could not open file for writing PNG: " << outPngFilename << "\n";
        return;
    }

    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png_ptr) {
        fclose(fp);
        return;
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
        png_destroy_write_struct(&png_ptr, nullptr);
        fclose(fp);
        return;
    }

    // Error handling setjmp block required by libpng specifications
    if (setjmp(png_jmpbuf(png_ptr))) {
        png_destroy_write_struct(&png_ptr, &info_ptr);
        fclose(fp);
        std::cerr << "[ERROR] Code triggered an internal libpng failure trap.\n";
        return;
    }

    png_init_io(png_ptr, fp);

    // Write PNG file header metrics
    png_set_IHDR(
        png_ptr,
        info_ptr,
        imgWidth, imgHeight,
        8,                          // Bit depth per channel
        PNG_COLOR_TYPE_RGB,         // Color layout
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
    );
    png_write_info(png_ptr, info_ptr);

    // 5. Construct row pointer array mapping directly down your contiguous vector sequence
    std::vector<png_bytep> rowPointers(imgHeight);
    for (int y = 0; y < imgHeight; ++y) {
        rowPointers[y] = (png_bytep)&imageBuffer[y * imgWidth * 3];
    }

    // Finalize stream to disk
    png_write_image(png_ptr, rowPointers.data());
    png_write_end(png_ptr, nullptr);

    // Clean up file and structural memory footprints
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);

    std::cout << "[READER] Native png.h export successful: " << outPngFilename << "\n";
}
