#ifndef GRAINPROCESSOR_H
#define GRAINPROCESSOR_H

#include <string>
#include <vector>
#include <array>

#include <spdlog/spdlog.h>
#include <Eigen/Core>
#include <gmsh.h>



struct Tetra
{
    Eigen::Vector3f nds[4];
    int grain;
};


class GrainProcessor
{
public:
    GrainProcessor();
    void LoadMSH(std::string fileName);

    std::vector<Eigen::Vector3f> vertices, vertices2;
    std::vector<std::array<int,5>> elems, elems2;   // 4 nodes + grain id
    std::vector<Tetra> tetra, tetra2;

};

#endif // GRAINPROCESSOR_H
