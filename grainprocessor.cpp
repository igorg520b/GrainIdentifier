#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>

#include <H5Cpp.h>
#include <Eigen/Dense>

#include "grainprocessor.h"


void GrainProcessor::Update_HDF5(std::string fileName)
{
    H5::H5File file(fileName, H5F_ACC_RDWR);
    hsize_t dims_grains[1] = {grainID.size()};
    H5::DataSpace dataspace_points(1, dims_grains);
    H5::DataSet dataset_points = file.createDataSet("GrainIDs", H5::PredType::NATIVE_INT16, dataspace_points);
    dataset_points.write(grainID.data(), H5::PredType::NATIVE_INT16);
    file.close();
}


bool GrainProcessor::PointInsideTetrahedron(Eigen::Vector3f point, Eigen::Vector3f t[4])
{
    constexpr float eps = 1e-3;
    Eigen::Matrix3f M,B;
    M << t[1]-t[0], t[2]-t[0], t[3]-t[0];
    B = M.inverse();
    point -= t[0];
    Eigen::Vector3f b = B * point; // barycentric
    return (b[0]>-eps && b[1]>-eps && b[2]>-eps && (b.sum() < 1+eps));
}



void GrainProcessor::IdentifyGrains()
{
    spdlog::info("identify grains");

    BVHN::BVHNFactory.release(leaves);
    leaves.reserve(tetra2.size());

    for(int i=0;i<tetra2.size();i++)
    {
        Tetra &t = tetra2[i];
        BVHN* bvhn = BVHN::BVHNFactory.take();
        bvhn->isLeaf = true;
        bvhn->box.Reset();
        for(int j=0;j<4;j++) bvhn->Expand(t.nds[j]);
        bvhn->child1 = bvhn->child2 = nullptr;
        bvhn->elem = i;
        leaves.push_back(bvhn);
    }

    spdlog::info("building bvh");
    root.Build(leaves, 0);
    spdlog::info("finished building bvh");

    // identify grains
    int unidentifiedPoints = 0;
    int problematicPoints = 0;

#pragma omp parallel for reduction(+:unidentifiedPoints,problematicPoints)
    for(int i=0;i<buffer.size();i++)
    {
        std::array<float, 3> &arr = buffer[i];
        Eigen::Vector3f v(arr[0],arr[1],arr[2]);
        v*=scale;
        for(int j=0;j<3;j++) v[j] = v[j] - floor(v[j]);

        BVHN bvhn;
        bvhn.isLeaf = true;
        bvhn.elem = -1;
        bvhn.box.Reset();
        bvhn.Expand(v);

        std::vector<std::pair<BVHN*,BVHN*>> broad_list;
        broad_list.reserve(10);
        bvhn.Collide(&root, broad_list);

        grainID[i] = -1;
        for(auto &pair : broad_list)
        {
            BVHN *bvhn2 = pair.second;
            int idx = bvhn2->elem;
            if(idx == -1) spdlog::critical("elem index -1");
            Tetra &t = tetra2[idx];
            bool result = PointInsideTetrahedron(v, t.nds);
            if(result)
            {
                grainID[i] = (short)t.grain;
                break;
            }
        }

        if(grainID[i] == -1)
        {
            problematicPoints++;
            spdlog::warn("node {}; grain {}; broad_list {}",i,grainID[i],broad_list.size());
            if(broad_list.size()!=0)
            {
                auto &p = broad_list.front();
                grainID[i] = (short)tetra2[p.second->elem].grain;
            }
        }
        if(grainID[i] == -1) { grainID[i] = 0; unidentifiedPoints++; }
    }
    spdlog::info("finished processing points; problematic {}; unidentified {}",problematicPoints, unidentifiedPoints);
}


void GrainProcessor::Load_Points_HDF5(std::string fileName)
{
    spdlog::info("ReadRawPoints {}",fileName);
    if(!std::filesystem::exists(fileName)) throw std::runtime_error("error reading raw points file - no file");;

    H5::H5File file(fileName, H5F_ACC_RDONLY);

    H5::DataSet dataset = file.openDataSet("Points_Raw");
    hsize_t dims[2] = {};
    dataset.getSpace().getSimpleExtentDims(dims, NULL);
    int nPoints = dims[0];
    if(dims[1]!=3) throw std::runtime_error("error reading raw points file - dimensions mismatch");
    spdlog::info("dims[0] {}, dims[1] {}", dims[0], dims[1]);

    buffer.resize(nPoints);
    dataset.read(buffer.data(), H5::PredType::NATIVE_FLOAT);
    file.close();

    grainID.resize(nPoints);

    spdlog::info("Load_Points_HDF5 done");
}





void GrainProcessor::LoadMSH(std::string fileName)
{
    vertices.clear();
    elems.clear();

    spdlog::info("load {}", fileName);
    if(!std::filesystem::exists(fileName)) spdlog::critical("file does not exist");

    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::open(fileName);


    // get nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;

    std::unordered_map<std::size_t, int> mtags; // gmsh nodeTag -> node sequential number

    gmsh::model::mesh::getNodesByElementType(4, nodeTags, nodeCoords, parametricCoords);
    vertices.reserve(nodeTags.size()*27);

    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");
        Eigen::Vector3f coords = Eigen::Vector3f(nodeCoords[i*3+0], nodeCoords[i*3+1], nodeCoords[i*3+2]);
        mtags[tag] = vertices.size();
        vertices.push_back(coords);
    }
    spdlog::info("nodes read {}",vertices.size());
    int nVerticesOriginal = vertices.size();
    vertices2.resize(nVerticesOriginal*27);

    // create surrounding copies
    int count = 0;
    for(int i=-1;i<=1;i++)
        for(int j=-1;j<=1;j++)
            for(int k=-1;k<=1;k++)
            {
                for(int m=0;m<nVerticesOriginal;m++)
                {
                    vertices2[m + nVerticesOriginal*count] = vertices[m];
                    vertices2[m + nVerticesOriginal*count][0] += (float)i;
                    vertices2[m + nVerticesOriginal*count][1] += (float)j;
                    vertices2[m + nVerticesOriginal*count][2] += (float)k;
                }
                count++;
            }


    std::vector<std::pair<int,int>> dimTagsGrains;
    gmsh::model::getEntities(dimTagsGrains,3);
    spdlog::info("dimTagsGrains size {}",dimTagsGrains.size());

//    std::unordered_set<int> used_nodes;
    elems.reserve(dimTagsGrains.size()*50);
    tetra.clear();
    tetra.reserve(elems.size());


    std::vector<int> types;
    std::vector<std::vector<std::size_t>> elemtags, nodetags;
    gmsh::model::mesh::getElements(types, elemtags, nodetags);

    for(int i=0;i<types.size();i++)
        spdlog::info("type {}; elemtags {}", types[i], elemtags[i].size());


    for(std::size_t j=0;j<dimTagsGrains.size();j++)
    {
        std::vector<std::size_t> tetraTags, nodeTagsInTetra;
        int entityTag = dimTagsGrains[j].second;
        gmsh::model::mesh::getElementsByType(4, tetraTags, nodeTagsInTetra, entityTag);

        for(std::size_t i=0;i<tetraTags.size();i++)
        {
            std::array<int,5> elem;
            Tetra t;

            for(int k=0;k<4;k++)
            {
                elem[k] = mtags.at(nodeTagsInTetra[i*4+k]);
                t.nds[k] = vertices[elem[k]];
            }
            t.grain = elem[4] = j-1;
            elems.push_back(elem);
            tetra.push_back(t);
        }
    }
    spdlog::info("tetra read {}", tetra.size());

    int nTetraOriginal = tetra.size();
    elems2.resize(nTetraOriginal*27);
    tetra2.resize(nTetraOriginal*27);

    count = 0;
    for(int i=-1;i<=1;i++)
        for(int j=-1;j<=1;j++)
            for(int k=-1;k<=1;k++)
            {
                for(int m=0;m<tetra.size();m++)
                {
                    int idx = m + count*nTetraOriginal;
                    elems2[idx] = elems[m];
                    Tetra &t = tetra2[idx];
                    for(int n=0;n<4;n++)
                    {
                        elems2[idx][n] += count*nVerticesOriginal;
                        t.nds[n] = vertices2[elems2[idx][n]];
                    }
                    t.grain = elems2[m][4];
                }
                count++;
            }

    gmsh::clear();
}




