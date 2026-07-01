/**
 * @file RadarCoverage.cpp
 *
 * RADAR COVERAGE ENGINE — with sector coverage support
 *
 * SECTOR COVERAGE:
 *   A sector is an arbitrary arc, given as two azimuths (degrees,
 *   North = 0°, clockwise positive):
 *
 *     leftDeg   = start of the arc (left edge, looking outward from source)
 *     rightDeg  = end of the arc   (right edge)
 *
 *   The scanned region is everything swept going clockwise from leftDeg
 *   to rightDeg. If leftDeg > rightDeg the arc wraps through North
 *   (e.g. leftDeg=350°, rightDeg=20° scans 350°→360°→0°→20°, a 30°-wide
 *   arc straddling North).
 *
 *   Pass to compute()/computeFast():
 *     leftDeg = -1, rightDeg = -1 (default)  → full 360°
 *     leftDeg = 300, rightDeg = 60           → 120°-wide arc centred on North
 *     leftDeg = 80,  rightDeg = 100          → narrow 20°-wide arc around East
 *
 *   Any azimuth value is accepted (negative or >360 are normalised).
 */

#include "../common.h"
#include "../operations.h"

#if defined(_OPENMP)
#include <omp.h>
#endif




/******************************************************************************
 * CIRCLE BOUNDARY
 ******************************************************************************/

std::vector<Point2D>
RadarEngine::circleBoundary(int cr, int cc, int radius)
{
    std::vector<Point2D> pts;
    pts.reserve(8 * radius);
    int x=0,y=radius,d=3-2*radius;
    auto push8=[&](int px,int py){
        pts.push_back({cr+px,cc+py}); pts.push_back({cr-px,cc+py});
        pts.push_back({cr+px,cc-py}); pts.push_back({cr-px,cc-py});
        pts.push_back({cr+py,cc+px}); pts.push_back({cr-py,cc+px});
        pts.push_back({cr+py,cc-px}); pts.push_back({cr-py,cc-px});
    };
    while(y>=x){
        push8(x,y); x++;
        if(d>0){y--;d+=4*(x-y)+10;}else{d+=4*x+6;}
    }
    std::sort(pts.begin(),pts.end(),[](const Point2D&a,const Point2D&b){
        return a.x<b.x||(a.x==b.x&&a.y<b.y);});
    pts.erase(std::unique(pts.begin(),pts.end(),
        [](const Point2D&a,const Point2D&b){return a.x==b.x&&a.y==b.y;}),pts.end());
    return pts;
}

/******************************************************************************
 * CORE RAY — DIAGNOSTIC
 ******************************************************************************/

static void castRadarRayDiag(
    int srcRow, int srcCol,
    int endRow, int endCol,
    double srcElev, double targetHeight,
    double minSlope, double maxSlope,
    double maxRangeM, double horizonM,
    int maxRadiusCells,
    const short*  __restrict__ dem,
    uint8_t*      __restrict__ stat,
    int COLS, int ROWS)
{
    int adr=std::abs(endRow-srcRow), adc=std::abs(endCol-srcCol);
    int sr=(srcRow<endRow)?1:-1, sc=(srcCol<endCol)?1:-1;
    int err=adr-adc;

    const std::ptrdiff_t strideR=static_cast<std::ptrdiff_t>(sr)*COLS;
    const std::ptrdiff_t strideC=sc;
    std::ptrdiff_t off=static_cast<std::ptrdiff_t>(srcRow)*COLS+srcCol;

    const short* demCell  = dem  + off;
    uint8_t*     statCell = stat + off;

    int row=srcRow,col=srcCol;
    double maxTerrainSlope=-1e18;
    int cellCount=0;

    while(true)
    {
        if(row==endRow && col==endCol) break;

        int e2=2*err;
        bool mr=(e2>-adc), mc=(e2<adr);
        if(mr){err-=adc;row+=sr;demCell+=strideR;statCell+=strideR;}
        if(mc){err+=adr;col+=sc;demCell+=strideC;statCell+=strideC;}

        ++cellCount;
        if(cellCount>maxRadiusCells) break;

        // True straight-line ground distance from source (not the
        // Bresenham step-accumulated "octile" distance, which traces
        // an octagon instead of a circle at any fixed threshold).
        double drow = static_cast<double>(row-srcRow);
        double dcol = static_cast<double>(col-srcCol);
        double dist = std::sqrt(drow*drow+dcol*dcol) * RESOLUTION;
        double curvDrop = (dist*dist) * INV_2_KEFFR;

        if(dist>maxRangeM){ *statCell=(uint8_t)CoverageStatus::OUTSIDE_RANGE; continue; }
        if(dist>horizonM) break;  // pre-filled as BEYOND_HORIZON

        short raw=*demCell;
        if(raw<-9000) break;

        double terrainElev  = static_cast<double>(raw) - curvDrop;
        double terrainSlope = (terrainElev - srcElev) / dist;
        double targetSlope  = terrainSlope + targetHeight / dist;

        uint8_t s;
        if     (maxTerrainSlope > targetSlope) s=(uint8_t)CoverageStatus::TERRAIN_MASKED;
        else if(targetSlope < minSlope)        s=(uint8_t)CoverageStatus::BELOW_MIN_ELEV;
        else if(targetSlope > maxSlope)        s=(uint8_t)CoverageStatus::ABOVE_MAX_ELEV;
        else                                   s=(uint8_t)CoverageStatus::COVERED;

        *statCell = s;
        if(terrainSlope > maxTerrainSlope) maxTerrainSlope = terrainSlope;
    }
}

/******************************************************************************
 * CORE RAY — FAST (binary)
 ******************************************************************************/

static void castRadarRayFast(
    int srcRow, int srcCol,
    int endRow, int endCol,
    double srcElev, double targetHeight,
    double minSlope, double maxSlope,
    double horizonM, int maxRadiusCells,
    const short* __restrict__ dem,
    uint8_t*     __restrict__ stat,
    int COLS, int ROWS)
{
    int adr=std::abs(endRow-srcRow), adc=std::abs(endCol-srcCol);
    int sr=(srcRow<endRow)?1:-1, sc=(srcCol<endCol)?1:-1;
    int err=adr-adc;

    const std::ptrdiff_t strideR=static_cast<std::ptrdiff_t>(sr)*COLS;
    const std::ptrdiff_t strideC=sc;
    std::ptrdiff_t off=static_cast<std::ptrdiff_t>(srcRow)*COLS+srcCol;

    const short* demCell  = dem  + off;
    uint8_t*     statCell = stat + off;

    int row=srcRow,col=srcCol;
    double maxTerrainSlope=-1e18;
    int cellCount=0;

    while(true)
    {
        if(row==endRow && col==endCol) break;

        int e2=2*err;
        bool mr=(e2>-adc), mc=(e2<adr);
        if(mr){err-=adc;row+=sr;demCell+=strideR;statCell+=strideR;}
        if(mc){err+=adr;col+=sc;demCell+=strideC;statCell+=strideC;}

        ++cellCount;
        if(cellCount>maxRadiusCells) break;

        // True straight-line ground distance from source (see castRadarRayDiag).
        double drow = static_cast<double>(row-srcRow);
        double dcol = static_cast<double>(col-srcCol);
        double dist = std::sqrt(drow*drow+dcol*dcol) * RESOLUTION;
        double curvDrop = (dist*dist) * INV_2_KEFFR;

        if(dist>horizonM) break;

        short raw=*demCell;
        if(raw<-9000) break;

        double terrainElev  = static_cast<double>(raw) - curvDrop;
        double terrainSlope = (terrainElev - srcElev) / dist;
        double targetSlope  = terrainSlope + targetHeight / dist;

        if(maxTerrainSlope<=targetSlope && targetSlope>=minSlope && targetSlope<=maxSlope)
            *statCell = 1;

        if(terrainSlope > maxTerrainSlope) maxTerrainSlope = terrainSlope;
    }
}

/******************************************************************************
 * RADAR COVERAGE ENGINE
 ******************************************************************************/



    std::vector<Point2D> RadarEngine::sectorBoundary(double leftDeg, double rightDeg) const
    {
        if(leftDeg < 0.0 || rightDeg < 0.0) return boundary;   // full 360

        std::vector<Point2D> out;
        // rough reservation based on arc width vs full circle
        double width = rightDeg - leftDeg;
        if(width < 0.0) width += 360.0;
        out.reserve(boundary.size() * (width / 360.0) + 8);

        for(const auto& pt : boundary){
            double az = cellAzimuth(srcRow, srcCol, pt.x, pt.y);
            if(inArc(az, leftDeg, rightDeg))
                out.push_back(pt);
        }
        return out;
    }

    // ------------------------------------------------------------------
    // Pre-fill the active area with BEYOND_HORIZON.
    // Cells outside the [leftDeg,rightDeg] arc get OUTSIDE_SECTOR.
    // leftDeg < 0 || rightDeg < 0 → fill full circle, no OUTSIDE_SECTOR
    // ------------------------------------------------------------------
    void RadarEngine::prefill(double leftDeg, double rightDeg, uint8_t* buf, uint8_t fillInside, uint8_t fillOutside)
    {
        const long long rSq = (long long)maxRadiusCells * maxRadiusCells;
        const bool fullCircle = (leftDeg < 0.0 || rightDeg < 0.0);

        // Normalise once outside the hot loop.
        const double L = fullCircle ? 0.0 : normAz(leftDeg);
        const double R = fullCircle ? 0.0 : normAz(rightDeg);
        const bool   wraps = (!fullCircle) && (L > R);

        #pragma omp parallel for schedule(static,32)
        for(int r=rowMin; r<rowMax; r++){
            int dr = r - srcRow;
            uint8_t* rowPtr = buf + (size_t)r * gridCols + colMin;

            for(int c=colMin; c<colMax; c++){
                int dc = c - srcCol;
                if((long long)dr*dr + (long long)dc*dc > rSq){
                    rowPtr[c - colMin] = (uint8_t)CoverageStatus::OUTSIDE_RANGE;
                    continue;
                }
                // inside circle
                if(fullCircle){
                    rowPtr[c - colMin] = fillInside;
                } else {
                    double az = cellAzimuth(srcRow, srcCol, r, c);
                    bool inside = wraps ? (az >= L || az <= R)
                                        : (az >= L && az <= R);
                    rowPtr[c - colMin] = inside ? fillInside : fillOutside;
                }
            }
        }
    }


    // ------------------------------------------------------------------
    // DIAGNOSTIC COMPUTE
    //
    // leftDeg/rightDeg: arc swept clockwise from leftDeg to rightDeg,
    //                    North = 0°. leftDeg < 0 || rightDeg < 0 → full 360.
    //                    leftDeg > rightDeg wraps through North.
    // ------------------------------------------------------------------
    RadarCoverageResult RadarEngine::compute(double leftDeg, double rightDeg)
    {
        bool fullCircle = (leftDeg < 0.0 || rightDeg < 0.0);

        auto activeBoundary = sectorBoundary(leftDeg, rightDeg);

        std::cout<<"[RADAR] Sector        : ";
        if(fullCircle) std::cout<<"360° (full)\n";
        else std::cout<<normAz(leftDeg)<<"° -> "<<normAz(rightDeg)
                      <<"° (clockwise)  rays="<<activeBoundary.size()<<"\n";

        auto t0=std::chrono::high_resolution_clock::now();

        // pre-fill: inside arc → BEYOND_HORIZON, outside → OUTSIDE_SECTOR
        prefill(leftDeg, rightDeg, statusBuf.data(),
                (uint8_t)CoverageStatus::BEYOND_HORIZON,
                (uint8_t)CoverageStatus::OUTSIDE_SECTOR);

        statusBuf[(size_t)srcRow*gridCols+srcCol]=(uint8_t)CoverageStatus::COVERED;

        const size_t N      = activeBoundary.size();
        const short* demPtr = g.bin.elev_matrix.data();
        uint8_t*     statPtr= statusBuf.data();
        const int COLS=gridCols,ROWS=gridRows;
        const double SE=srcElev,TH=targetHeight;
        const double MNS=minSlopeRatio,MXS=maxSlopeRatio;
        const double MXR=maxRangeM,HOR=horizonM;
        const int    MRC=maxRadiusCells;

        #pragma omp parallel for schedule(dynamic,64) default(none) \
            shared(demPtr,statPtr,activeBoundary) \
            firstprivate(N,COLS,ROWS,SE,TH,MNS,MXS,MXR,HOR,MRC)
        for(size_t i=0;i<N;++i){
            int er=activeBoundary[i].x, ec=activeBoundary[i].y;
            if(er>=0&&ec>=0&&er<ROWS&&ec<COLS)
                castRadarRayDiag(
                    srcRow,srcCol,er,ec,
                    SE,TH,MNS,MXS,MXR,HOR,MRC,
                    demPtr,statPtr,COLS,ROWS);
        }

        auto t1=std::chrono::high_resolution_clock::now();

        // stats — only count cells in active area
        RadarCoverageResult res;
        res.srcRow=srcRow; res.srcCol=srcCol;
        res.radiusCells=maxRadiusCells;
        res.gridCols=gridCols; res.gridRows=gridRows;
        res.sectorLeftDeg  = fullCircle ? -1.0 : normAz(leftDeg);
        res.sectorRightDeg = fullCircle ? -1.0 : normAz(rightDeg);
        res.countCovered=res.countTerrainMasked=0;
        res.countBelowMinElev=res.countAboveMaxElev=0;
        res.countBeyondHorizon=res.countOutsideRange=res.countOutsideSector=0;

        long long rSq=(long long)maxRadiusCells*maxRadiusCells;
        for(int r=rowMin;r<rowMax;r++){
        for(int c=colMin;c<colMax;c++){
            int dr=r-srcRow,dc=c-srcCol;
            if((long long)dr*dr+(long long)dc*dc>rSq) continue;
            switch((CoverageStatus)statusBuf[(size_t)r*gridCols+c]){
                case CoverageStatus::COVERED:        res.countCovered++;        break;
                case CoverageStatus::TERRAIN_MASKED: res.countTerrainMasked++;  break;
                case CoverageStatus::BELOW_MIN_ELEV: res.countBelowMinElev++;   break;
                case CoverageStatus::ABOVE_MAX_ELEV: res.countAboveMaxElev++;   break;
                case CoverageStatus::BEYOND_HORIZON: res.countBeyondHorizon++;  break;
                case CoverageStatus::OUTSIDE_SECTOR: res.countOutsideSector++;  break;
                default:                             res.countOutsideRange++;   break;
            }
        }
        }

        res.status.resize((size_t)gridRows*gridCols);
        std::memcpy(res.status.data(),statusBuf.data(),(size_t)gridRows*gridCols);

        int total=res.countCovered+res.countTerrainMasked+
                  res.countBelowMinElev+res.countAboveMaxElev+
                  res.countBeyondHorizon;
        double ms=std::chrono::duration<double,std::milli>(t1-t0).count();

        std::cout<<"[RADAR] Compute       : "<<ms<<" ms\n";
        std::cout<<"  COVERED        : "<<res.countCovered<<"\n";
        std::cout<<"  TERRAIN_MASKED : "<<res.countTerrainMasked<<"\n";
        std::cout<<"  BELOW_MIN_ELEV : "<<res.countBelowMinElev<<"\n";
        std::cout<<"  ABOVE_MAX_ELEV : "<<res.countAboveMaxElev<<"\n";
        std::cout<<"  BEYOND_HORIZON : "<<res.countBeyondHorizon<<"\n";
        std::cout<<"  OUTSIDE_SECTOR : "<<res.countOutsideSector<<"\n";
        if(total>0) std::cout<<"  Coverage %     : "
                             <<(100.0*res.countCovered/total)<<"%\n";
        return res;
    }

    // ------------------------------------------------------------------
    // FAST COMPUTE — binary, no diagnostic
    // ------------------------------------------------------------------
    std::vector<uint8_t> RadarEngine::computeFast(double leftDeg, double rightDeg)
    {
        bool fullCircle = (leftDeg < 0.0 || rightDeg < 0.0);
        auto activeBoundary = sectorBoundary(leftDeg, rightDeg);

        auto t0=std::chrono::high_resolution_clock::now();

        // reset bounding box to 0
        for(int r=rowMin;r<rowMax;r++){
            uint8_t* p=fastBuf.data()+(size_t)r*gridCols+colMin;
            std::memset(p,0,colMax-colMin);
        }
        fastBuf[(size_t)srcRow*gridCols+srcCol]=1;

        const size_t N      = activeBoundary.size();
        const short* demPtr = g.bin.elev_matrix.data();
        uint8_t*     statPtr= fastBuf.data();
        const int COLS=gridCols,ROWS=gridRows;
        const double SE=srcElev,TH=targetHeight;
        const double MNS=minSlopeRatio,MXS=maxSlopeRatio;
        const double HOR=horizonM;
        const int    MRC=maxRadiusCells;

        #pragma omp parallel for schedule(dynamic,64) default(none) \
            shared(demPtr,statPtr,activeBoundary) \
            firstprivate(N,COLS,ROWS,SE,TH,MNS,MXS,HOR,MRC)
        for(size_t i=0;i<N;++i){
            int er=activeBoundary[i].x,ec=activeBoundary[i].y;
            if(er>=0&&ec>=0&&er<ROWS&&ec<COLS)
                castRadarRayFast(
                    srcRow,srcCol,er,ec,
                    SE,TH,MNS,MXS,HOR,MRC,
                    demPtr,statPtr,COLS,ROWS);
        }

        auto t1=std::chrono::high_resolution_clock::now();
        int vis=0;
        for(int r=rowMin;r<rowMax;r++){
            for(int c=colMin;c<colMax;c++)
                vis+=fastBuf[(size_t)r*gridCols+c];
        }

        std::cout<<"[RADAR] Fast ("
                 <<(fullCircle ? std::string("360°")
                               : std::to_string(normAz(leftDeg))+"->"+std::to_string(normAz(rightDeg))+"°")
                 <<")  : "
                 <<std::chrono::duration<double,std::milli>(t1-t0).count()
                 <<" ms  covered="<<vis<<"\n";

        return fastBuf;
    }

/******************************************************************************
 * MAIN
 ******************************************************************************/

// int main()
// {
//     auto t0=std::chrono::high_resolution_clock::now();
//     Global g;

//     RadarEngine engine(
//         g,
//         33.987, 74.7742, // Srinagar airport
//         20.0,    // antenna height AGL (m)
//         0.1,  // target height AGL (m) — e.g. aircraft at 1km
//         0.5,     // min elevation angle (degrees)
//         75.0,    // max elevation angle (degrees)
//         300.0    // max range (km)
//     );

//     // // full 360
//     auto res360 = engine.compute();
//     engine.writeToSVG(res360, "radar_360.svg");

//     // // 60°-wide arc centred on NE (20° -> 80°)
//     // auto resNE = engine.compute(20, 80);
//     // engine.writeToSVG(resNE, "radar_sector_ne.svg");

//     // 40°-wide arc straddling North (340° -> 20°)
//     auto resN = engine.compute(300, 60);
//     engine.writeToSVG(resN, "330.svg");

//     // // fast mode
//     // engine.computeFast();          // 360
//     // engine.computeFast(110, 200);  // SE-ish arc

//     auto t1=std::chrono::high_resolution_clock::now();
//     std::cout<<"\n==============================================\n"
//              <<">> TOTAL: "
//              <<std::chrono::duration<double,std::milli>(t1-t0).count()/1000.0
//              <<" s\n==============================================\n";
// }
