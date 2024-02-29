#ifndef GRAINPROCESSOR_H
#define GRAINPROCESSOR_H

#include <string>
#include <vector>
#include <array>

#include <spdlog/spdlog.h>
#include <Eigen/Core>
#include <gmsh.h>

#include "bvhn.h"

struct Tetra
{
    Eigen::Vector3f nds[4];
    int grain;
};


class GrainProcessor
{
    static constexpr float scale = 2;

public:
    void LoadMSH(std::string fileName);
    void Load_Points_HDF5(std::string fileName);
    void IdentifyGrains();
    void Update_HDF5(std::string fileName);
    static bool PointInsideTetrahedron(Eigen::Vector3f point, Eigen::Vector3f tetra[4]);

    std::vector<Eigen::Vector3f> vertices, vertices2;
    std::vector<std::array<int,5>> elems, elems2;   // 4 nodes + grain id
    std::vector<Tetra> tetra, tetra2;

    std::vector<std::array<float, 3>> buffer;
    std::vector<short> grainID;

    std::vector<BVHN*> leaves;
    BVHN root;

};

#endif // GRAINPROCESSOR_H
