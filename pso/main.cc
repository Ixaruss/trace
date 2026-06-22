
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>
#include <iomanip>
#include <string>
#include <fstream>
#include <mutex>
#include <memory>
#include <queue>
#include <thread>
#include <condition_variable>
#include <functional>
#include <atomic>

// ============================================================================
// DATA STRUCTURES & DATA MODELS
// ============================================================================
struct BoundingBox {
    double minX, maxX, minY, maxY;
};

struct RadarPhysicsConfig {
    double maxRangeUnits;
};

struct PsoHyperparameters {
    double c1;
    double c2;
    double defaultW;
};

#pragma pack(push, 1)
struct GridPoint {
    double x;
    double y;

    double distanceTo(const GridPoint& other) const {
        double dx = other.x - x;
        double dy = other.y - y;
        return std::sqrt(dx * dx + dy * dy);
    }
};
#pragma pack(pop)

struct ContinuousParticle {
    GridPoint position;
    double velX;
    double velY;
    GridPoint pBestPosition;
    double pBestFitness;
    double currentFitness;
};

struct PSOBenchmarkResult {
    std::string profileName;
    std::string topology;
    std::string inertia;
    std::string boundary;
    int swarmSize;
    int iterations;
    double maxCoveragePercentage;
    GridPoint bestLocation;
    double totalTimeMs;
    std::vector<GridPoint> initialPositions;
};

// ============================================================================
// LIGHTWEIGHT WORKER THREAD POOL ENGINE
// ============================================================================
class ThreadPool {
private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;
    std::mutex queueMutex;
    std::condition_variable cv;
    bool stop = false;

public:
    ThreadPool(size_t threads) {
        for (size_t i = 0; i < threads; ++i) {
            workers.emplace_back([this] {
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock<std::mutex> lock(this->queueMutex);
                        this->cv.wait(lock, [this] { return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty()) return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    task();
                }
            });
        }
    }

    template<class F>
    void enqueue(F&& f) {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            tasks.emplace(std::forward<F>(f));
        }
        cv.notify_one();
    }

    ~ThreadPool() {
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            stop = true;
        }
        cv.notify_all();
        for (std::thread &worker : workers) {
            if (worker.joinable()) worker.join();
        }
    }
};

// ============================================================================
// NATIVE CONFIGURATION LOADER
// ============================================================================
struct ConfigLoader {
    static bool loadConfig(const std::string& filename, RadarPhysicsConfig& physics,
                            PsoHyperparameters& hyper, int& swarmSize, int& iterations) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << ">> [CONFIG ERROR] Could not open " << filename << ".\n";
            return false;
        }

        std::string line;
        auto extractDouble = [](const std::string& l) {
            size_t colon = l.find(':');
            if (colon == std::string::npos) return 0.0;
            size_t comma = l.find(',');
            size_t endIdx = (comma == std::string::npos) ? l.find_last_not_of(" \t\r\n}") + 1 : comma;
            return std::stod(l.substr(colon + 1, endIdx - colon - 1));
        };

        while (std::getline(file, line)) {
            if (line.find("max_range_units") != std::string::npos || line.find("max_range_km") != std::string::npos) {
                physics.maxRangeUnits = extractDouble(line);
            }
            else if (line.find("cognitive_weight_c1") != std::string::npos) hyper.c1 = extractDouble(line);
            else if (line.find("social_weight_c2") != std::string::npos) hyper.c2 = extractDouble(line);
            else if (line.find("inertia_weight_w") != std::string::npos) hyper.defaultW = extractDouble(line);
            else if (line.find("swarm_size") != std::string::npos) swarmSize = static_cast<int>(extractDouble(line));
            else if (line.find("iterations") != std::string::npos) iterations = static_cast<int>(extractDouble(line));
        }

        std::cout << ">> [CONFIG LOADER] Engine Parameters Configured:\n"
                  << "   - Max Range Radius: " << physics.maxRangeUnits << " grid units\n"
                  << "   - Swarm Allocation: " << swarmSize << " particles | Horizon: " << iterations << " generations\n";
        return true;
    }
};

// ============================================================================
// ABSTRACT POLYMORPHIC STRATEGY INTERFACES
// ============================================================================
class IInertiaStrategy {
public:
    virtual ~IInertiaStrategy() = default;
    virtual double computeWeight(int currentIter, int maxIter, double defaultW) = 0;
    virtual std::string getName() const = 0;
};

class IBoundaryHandler {
public:
    virtual ~IBoundaryHandler() = default;
    virtual void handle(double& x, double& y, double& vx, double& vy, const BoundingBox& bbox) = 0;
    virtual std::string getName() const = 0;
};

class ITopologySelector {
public:
    virtual ~ITopologySelector() = default;
    virtual GridPoint selectTarget(int idx, int swarmSize, const GridPoint& gBest, const std::vector<ContinuousParticle>& swarm) = 0;
    virtual std::string getName() const = 0;
};

// ============================================================================
// CONCRETE STRATEGY IMPLEMENTATIONS
// ============================================================================
class ConstantInertia : public IInertiaStrategy {
public:
    double computeWeight(int, int, double defaultW) override { return defaultW; }
    std::string getName() const override { return "Constant"; }
};

class LinearDecayInertia : public IInertiaStrategy {
public:
    double computeWeight(int currentIter, int maxIter, double) override {
        return 0.9 - (0.5 * (static_cast<double>(currentIter) / maxIter));
    }
    std::string getName() const override { return "Linear Decay"; }
};

class RandomJitterInertia : public IInertiaStrategy {
private:
    std::mt19937 rng{1337};
public:
    double computeWeight(int, int, double) override {
        std::uniform_real_distribution<double> jitter(0.4, 0.9);
        return jitter(rng);
    }
    std::string getName() const override { return "Rand Jitter"; }
};

class AbsorbClampHandler : public IBoundaryHandler {
public:
    void handle(double& x, double& y, double&, double&, const BoundingBox& bbox) override {
        x = std::clamp(x, bbox.minX, bbox.maxX);
        y = std::clamp(y, bbox.minY, bbox.maxY);
    }
    std::string getName() const override { return "Clamp"; }
};

class ReflectBounceHandler : public IBoundaryHandler {
public:
    void handle(double& x, double& y, double& vx, double& vy, const BoundingBox& bbox) override {
        if (x < bbox.minX) { x = bbox.minX + (bbox.minX - x); vx *= -1.0; }
        if (x > bbox.maxX) { x = bbox.maxX - (x - bbox.maxX); vx *= -1.0; }
        if (y < bbox.minY) { y = bbox.minY + (bbox.minY - y); vy *= -1.0; }
        if (y > bbox.maxY) { y = bbox.maxY - (y - bbox.maxY); vy *= -1.0; }
        x = std::clamp(x, bbox.minX, bbox.maxX);
        y = std::clamp(y, bbox.minY, bbox.maxY);
    }
    std::string getName() const override { return "Bounce"; }
};

class GlobalStarTopology : public ITopologySelector {
public:
    GridPoint selectTarget(int, int, const GridPoint& gBest, const std::vector<ContinuousParticle>&) override {
        return gBest;
    }
    std::string getName() const override { return "Global (Star)"; }
};

class LocalRingTopology : public ITopologySelector {
public:
    GridPoint selectTarget(int idx, int swarmSize, const GridPoint&, const std::vector<ContinuousParticle>& swarm) override {
        int leftNeighbor = (idx - 1 + swarmSize) % swarmSize;
        int rightNeighbor = (idx + 1) % swarmSize;

        double maxLocalBest = swarm[idx].pBestFitness;
        GridPoint targetPosition = swarm[idx].pBestPosition;

        if (swarm[leftNeighbor].pBestFitness > maxLocalBest) {
            maxLocalBest = swarm[leftNeighbor].pBestFitness;
            targetPosition = swarm[leftNeighbor].pBestPosition;
        }
        if (swarm[rightNeighbor].pBestFitness > maxLocalBest) {
            targetPosition = swarm[rightNeighbor].pBestPosition;
        }
        return targetPosition;
    }
    std::string getName() const override { return "Local (Ring)"; }
};

// ============================================================================
// STRATEGY-DRIVEN PARALLEL GRID PSO ENGINE
// ============================================================================
class ParallelPolygonRadarPSO {
private:
    BoundingBox bbox;
    int swarmSize, iterations;
    RadarPhysicsConfig radarPhysics;
    PsoHyperparameters psoHyper;
    std::mt19937 rng;

    std::shared_ptr<ITopologySelector> activeTopology;
    std::shared_ptr<IInertiaStrategy> activeInertia;
    std::shared_ptr<IBoundaryHandler> activeBoundary;
    std::string profileName;

    GridPoint gBestPosition;
    double gBestFitness;
    std::mutex gBestMutex;

    // Persistent runtime thread pool initialization
    ThreadPool pool{std::max(1u, std::thread::hardware_concurrency())};

    void initializeSwarm(std::vector<ContinuousParticle>& swarm) {
        double maxVelX = (bbox.maxX - bbox.minX) * 0.01;
        double maxVelY = (bbox.maxY - bbox.minY) * 0.01;

        std::uniform_real_distribution<double> vDistX(-maxVelX, maxVelX);
        std::uniform_real_distribution<double> vDistY(-maxVelY, maxVelY);

        // Dynamic random selection mapping profile
        std::uniform_int_distribution<size_t> meshIndexDist(0, filteredSamplingMesh.size() - 1);

        for (int i = 0; i < swarmSize; ++i) {
            size_t randomIndex = meshIndexDist(rng);
            swarm[i].position = filteredSamplingMesh[randomIndex];
            swarm[i].velX = vDistX(rng);
            swarm[i].velY = vDistY(rng);
            swarm[i].currentFitness = 0.0;
            swarm[i].pBestFitness = -1.0;
        }
    }

    void parallelEvaluateSwarm(std::vector<ContinuousParticle>& swarm) {
        if (filteredSamplingMesh.empty() || swarm.empty()) return;

        unsigned int numThreads = std::thread::hardware_concurrency();
        if (numThreads == 0) numThreads = 4;

        size_t totalParticles = swarm.size();
        size_t chunkSize = totalParticles / numThreads;
        if (chunkSize == 0) chunkSize = 1;

        size_t totalPoints = filteredSamplingMesh.size();

        std::atomic<int> completedChunks(0);
        int expectedChunks = 0;

        for (size_t startParticle = 0; startParticle < totalParticles; startParticle += chunkSize) {
            size_t endParticle = std::min(startParticle + chunkSize, totalParticles);
            expectedChunks++;

            pool.enqueue([this, &swarm, startParticle, endParticle, totalPoints, &completedChunks]() {
                for (size_t p = startParticle; p < endParticle; ++p) {
                    int localCount = 0;
                    for (size_t i = 0; i < totalPoints; ++i) {
                        if (swarm[p].position.distanceTo(filteredSamplingMesh[i]) <= radarPhysics.maxRangeUnits) {
                            localCount++;
                        }
                    }
                    swarm[p].currentFitness = static_cast<double>(localCount) / totalPoints;
                }
                completedChunks++;
            });
        }

        while (completedChunks < expectedChunks) {
            std::this_thread::yield();
        }
    }

public:
    std::vector<GridPoint> filteredSamplingMesh;

    ParallelPolygonRadarPSO(RadarPhysicsConfig physics, PsoHyperparameters hyper)
        : radarPhysics(physics), psoHyper(hyper) {
        rng.seed(1337);
        generateSamplingMesh();
    }

    void generateSamplingMesh() {
        filteredSamplingMesh.clear();
        std::string binPath = "feasible.bin";
        std::ifstream binFile(binPath, std::ios::binary | std::ios::ate);

        if (!binFile.is_open()) { binPath = "pso/feasible.bin"; binFile.open(binPath, std::ios::binary | std::ios::ate); }
        if (!binFile.is_open()) { binPath = "../feasible.bin"; binFile.open(binPath, std::ios::binary | std::ios::ate); }

        if (binFile.is_open()) {
            std::streamsize fileSize = binFile.tellg();
            binFile.seekg(0, std::ios::beg);

            size_t pointCount = fileSize / sizeof(GridPoint);
            filteredSamplingMesh.resize(pointCount);
            binFile.read(reinterpret_cast<char*>(filteredSamplingMesh.data()), fileSize);
            binFile.close();

            bbox.minX = 0.0;
            bbox.minY = 0.0;
            bbox.maxX = 0.0;
            bbox.maxY = 0.0;

            for (const auto& pt : filteredSamplingMesh) {
                if (pt.x > bbox.maxX) bbox.maxX = pt.x;
                if (pt.y > bbox.maxY) bbox.maxY = pt.y;
            }

            std::cout << ">> [RECTANGLE GRID] Loaded " << pointCount << " spatial nodes from feasible.bin\n"
                      << "   - Matrix Frame Bounds: [0 to " << bbox.maxX << "] x [0 to " << bbox.maxY << "]\n";
            return;
        }

        double fallbackN = 1500.0;
        double fallbackM = 1500.0;
        std::cerr << ">> [WARNING] feasible.bin asset skipped. Initializing fallback " << fallbackN << "x" << fallbackM << " matrix.\n";
        bbox = {0.0, fallbackN, 0.0, fallbackM};

        filteredSamplingMesh.reserve(static_cast<size_t>(fallbackN * fallbackM));
        for (int x = 0; x < static_cast<int>(fallbackN); ++x) {
            for (int y = 0; y < static_cast<int>(fallbackM); ++y) {
                filteredSamplingMesh.push_back({static_cast<double>(x), static_cast<double>(y)});
            }
        }
    }

    BoundingBox getDerivedBounds() const { return bbox; }

    void setRuntimeStrategy(const std::string& name, int size, int iters,
                            std::shared_ptr<ITopologySelector> topo,
                            std::shared_ptr<IInertiaStrategy> inertia,
                            std::shared_ptr<IBoundaryHandler> boundary) {
        profileName = name; swarmSize = size; iterations = iters;
        activeTopology = topo; activeInertia = inertia; activeBoundary = boundary;
    }

    PSOBenchmarkResult runPSO() {
        gBestFitness = -1.0;
        auto totalExecutionStart = std::chrono::high_resolution_clock::now();

        std::vector<ContinuousParticle> swarm(swarmSize);
        initializeSwarm(swarm);

        std::vector<GridPoint> initialPositionsSnapshot;
        initialPositionsSnapshot.reserve(swarmSize);
        for (const auto& p : swarm) {
            initialPositionsSnapshot.push_back(p.position);
        }

        parallelEvaluateSwarm(swarm);

        gBestPosition = swarm[0].position;
        gBestFitness = swarm[0].currentFitness;

        for (int i = 0; i < swarmSize; ++i) {
            swarm[i].pBestPosition = swarm[i].position;
            swarm[i].pBestFitness = swarm[i].currentFitness;
            if (swarm[i].pBestFitness > gBestFitness) {
                gBestFitness = swarm[i].pBestFitness;
                gBestPosition = swarm[i].pBestPosition;
            }
        }

        std::uniform_real_distribution<double> dist01(0.0, 1.0);
        double maxVelX = (bbox.maxX - bbox.minX) * 0.01;
        double maxVelY = (bbox.maxY - bbox.minY) * 0.01;

        for (int iter = 1; iter <= iterations; ++iter) {
            double w = activeInertia->computeWeight(iter, iterations, psoHyper.defaultW);

            for (int i = 0; i < swarmSize; ++i) {
                GridPoint targetPosition = activeTopology->selectTarget(i, swarmSize, gBestPosition, swarm);

                double r1 = dist01(rng), r2 = dist01(rng);
                swarm[i].velX = w * swarm[i].velX + psoHyper.c1 * r1 * (swarm[i].pBestPosition.x - swarm[i].position.x) + psoHyper.c2 * r2 * (targetPosition.x - swarm[i].position.x);
                swarm[i].velY = w * swarm[i].velY + psoHyper.c1 * r1 * (swarm[i].pBestPosition.y - swarm[i].position.y) + psoHyper.c2 * r2 * (targetPosition.y - swarm[i].position.y);

                swarm[i].velX = std::clamp(swarm[i].velX, -maxVelX, maxVelX);
                swarm[i].velY = std::clamp(swarm[i].velY, -maxVelY, maxVelY);

                swarm[i].position.x += swarm[i].velX;
                swarm[i].position.y += swarm[i].velY;

                activeBoundary->handle(swarm[i].position.x, swarm[i].position.y, swarm[i].velX, swarm[i].velY, bbox);
            }

            parallelEvaluateSwarm(swarm);

            for (int i = 0; i < swarmSize; ++i) {
                if (swarm[i].currentFitness > swarm[i].pBestFitness) {
                    swarm[i].pBestFitness = swarm[i].currentFitness;
                    swarm[i].pBestPosition = swarm[i].position;
                }
                if (swarm[i].currentFitness >= gBestFitness) {
                    std::lock_guard<std::mutex> lock(gBestMutex);
                    if (swarm[i].currentFitness >= gBestFitness) {
                        gBestFitness = swarm[i].currentFitness;
                        gBestPosition = swarm[i].position;
                    }
                }
            }
        }

        auto totalExecutionEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> totalExecutionDuration = totalExecutionEnd - totalExecutionStart;

        return { profileName, activeTopology->getName(), activeInertia->getName(), activeBoundary->getName(),
                 swarmSize, iterations, (gBestFitness < 0.0 ? 0.0 : gBestFitness * 100.0), gBestPosition,
                 totalExecutionDuration.count(), initialPositionsSnapshot };
    }
};

// ============================================================================
// CLEAN GEOMETRIC AOI SVG EXPORT MACHINE
// ============================================================================
void exportTopSolutionToSVG(const std::string& filename, const BoundingBox& box,
                            const std::vector<GridPoint>&, const GridPoint& radarPos,
                            const RadarPhysicsConfig& physics, const std::vector<GridPoint>& initialParticles) {
    std::ofstream svg(filename);
    if (!svg.is_open()) return;

    int canvasWidth = 1000, canvasHeight = 700, padding = 60;
    int activeWidthPx = canvasWidth - (2 * padding);
    int activeHeightPx = canvasHeight - (2 * padding);

    double gridWidth = box.maxX - box.minX;
    double gridHeight = box.maxY - box.minY;

    if (gridWidth <= 0) gridWidth = 1.0;
    if (gridHeight <= 0) gridHeight = 1.0;

    double scaleX = static_cast<double>(activeWidthPx) / gridWidth;
    double scaleY = static_cast<double>(activeHeightPx) / gridHeight;
    double pixelsPerUnit = std::min(scaleX, scaleY);

    double offsetX = padding + (activeWidthPx - (gridWidth * pixelsPerUnit)) / 2.0;
    double offsetY = padding + (activeHeightPx - (gridHeight * pixelsPerUnit)) / 2.0;

    auto mapX = [&](double x) { return offsetX + (x - box.minX) * pixelsPerUnit; };
    auto mapY = [&](double y) { return offsetY + (gridHeight - (y - box.minY)) * pixelsPerUnit; };

    double rangeInPx = physics.maxRangeUnits * pixelsPerUnit;

    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << canvasWidth << "\" height=\"" << canvasHeight << "\" style=\"background-color:#0F172A;\">\n";
    svg << "<defs><filter id=\"glow\"><feGaussianBlur stdDeviation=\"2\" result=\"coloredBlur\"/><feMerge><feMergeNode in=\"coloredBlur\"/><feMergeNode in=\"SourceGraphic\"/></feMerge></filter></defs>\n";

    svg << "<text x=\"30\" y=\"35\" fill=\"#F8FAFC\" font-family=\"monospace\" font-size=\"16\" font-weight=\"bold\">PSO PLACEMENT MONITOR</text>\n";
    svg << "<text x=\"30\" y=\"55\" fill=\"#94A3B8\" font-family=\"monospace\" font-size=\"12\">AOI Bounds: [0 to "
        << static_cast<int>(box.maxX) << " X] x [0 to " << static_cast<int>(box.maxY) << " Y]</text>\n";

    // 1. DRAW THE AOI (Solid, clean rectangular background area)
    double bx = mapX(box.minX);
    double by = mapY(box.maxY);
    double bw = gridWidth * pixelsPerUnit;
    double bh = gridHeight * pixelsPerUnit;

    svg << "<rect x=\"" << bx << "\" y=\"" << by << "\" width=\"" << bw << "\" height=\"" << bh
        << "\" fill=\"#1E293B\" stroke=\"#475569\" stroke-width=\"2\" />\n";

    // 2. DRAW RADAR COVERAGE CIRCLE (Semi-translucent blue bubble)
    double rx = mapX(radarPos.x), ry = mapY(radarPos.y);
    svg << "<circle cx=\"" << rx << "\" cy=\"" << ry << "\" r=\"" << rangeInPx
        << "\" stroke=\"#38BDF8\" stroke-width=\"2.5\" stroke-opacity=\"0.9\" fill=\"#0284C7\" fill-opacity=\"0.2\" />\n";

    // 3. DRAW INITIAL POPULATION PARTICLES (Scattered points inside the AOI)
    for (const auto& pt : initialParticles) {
        double px = mapX(pt.x), py = mapY(pt.y);
        svg << "<circle cx=\"" << px << "\" cy=\"" << py << "\" r=\"5\" fill=\"#38BDF8\" opacity=\"0.4\" filter=\"url(#glow)\"/>\n";
        svg << "<circle cx=\"" << px << "\" cy=\"" << py << "\" r=\"2.5\" fill=\"#FFFFFF\" stroke=\"#0F172A\" stroke-width=\"1.0\"/>\n";
    }

    // 4. DRAW FINAL OPTIMIZED RADAR CENTER
    svg << "<circle cx=\"" << rx << "\" cy=\"" << ry << "\" r=\"7\" fill=\"#F59E0B\" stroke=\"#FFFFFF\" stroke-width=\"2.0\"/>\n";
    svg << "<circle cx=\"" << rx << "\" cy=\"" << ry << "\" r=\"2\" fill=\"#FFFFFF\"/>\n";

    svg << "</svg>\n";
    svg.close();

    std::cout << ">> [SVG EXPORT] Clean geometric map saved to: " << filename << "\n";
}

// ============================================================================
// SYSTEM RUNTIME EXECUTION PIPELINE
// ============================================================================
int main() {
    RadarPhysicsConfig physics;
    PsoHyperparameters hyper;
    int swarmSize = 40, iterations = 50;

    if (!ConfigLoader::loadConfig("pso/pso_config.json", physics, hyper, swarmSize, iterations)) {
        if (!ConfigLoader::loadConfig("pso_config.json", physics, hyper, swarmSize, iterations)) {
            return 1;
        }
    }

    ParallelPolygonRadarPSO psoEngine(physics, hyper);
    BoundingBox derivedBounds = psoEngine.getDerivedBounds();

    std::vector<std::shared_ptr<ITopologySelector>> topologies = {
        std::make_shared<GlobalStarTopology>(),
        std::make_shared<LocalRingTopology>()
    };

    std::vector<std::shared_ptr<IInertiaStrategy>> inertias = {
        std::make_shared<ConstantInertia>(),
        std::make_shared<LinearDecayInertia>(),
        std::make_shared<RandomJitterInertia>()
    };

    std::vector<std::shared_ptr<IBoundaryHandler>> boundaries = {
        std::make_shared<AbsorbClampHandler>(),
        std::make_shared<ReflectBounceHandler>()
    };

    std::vector<PSOBenchmarkResult> globalMetricsTable;

    int profileId = 1;
    for (auto& topo : topologies) {
        for (auto& inertia : inertias) {
            for (auto& bound : boundaries) {
                std::string label = "Profile_" + std::to_string(profileId++);
                psoEngine.setRuntimeStrategy(label, swarmSize, iterations, topo, inertia, bound);
                globalMetricsTable.push_back(psoEngine.runPSO());
            }
        }
    }

    std::sort(globalMetricsTable.begin(), globalMetricsTable.end(), [](const PSOBenchmarkResult& a, const PSOBenchmarkResult& b) {
        return a.maxCoveragePercentage > b.maxCoveragePercentage;
    });

    std::cout << "\n### Strategy Benchmarking Table (Pure Grid-Space View)\n\n";
    std::cout << "| Rank | Topology Strategy | Inertia Model | Boundary Action | Max Coverage | Compute Time |\n";
    std::cout << "|:----:|:------------------|:--------------|:----------------|:------------:|:------------:|\n";

    int rank = 1;
    for (const auto& metrics : globalMetricsTable) {
        std::cout << "| " << std::setw(4) << rank++ << " "
                  << "| " << std::setw(17) << std::left << metrics.topology
                  << "| " << std::setw(13) << std::left << metrics.inertia
                  << "| " << std::setw(14) << std::left << metrics.boundary
                  << "| " << std::setw(10) << std::right << std::fixed << std::setprecision(2) << metrics.maxCoveragePercentage << "% "
                  << "| " << std::setw(9) << std::fixed << std::setprecision(1) << metrics.totalTimeMs << " ms |\n";
    }
    std::cout << "\n";

    if (!globalMetricsTable.empty() && !psoEngine.filteredSamplingMesh.empty()) {
        exportTopSolutionToSVG("pso/top_solution.svg",
                               derivedBounds,
                               psoEngine.filteredSamplingMesh,
                               globalMetricsTable[0].bestLocation,
                               physics,
                               globalMetricsTable[0].initialPositions);
    }

    return 0;
}
