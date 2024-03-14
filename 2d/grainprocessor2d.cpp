#include "grainprocessor2d.h"
#include "poisson_disk_sampling.h"

#include <vector>
#include <unordered_map>
#include <filesystem>

#include <H5Cpp.h>
#include <gmsh.h>
#include <spdlog/spdlog.h>

#include <Eigen/Dense>


void GrainProcessor2D::generate_block_and_write(float scale, float bx, float by, int n, std::string msh, std::string outputFile)
{
    LoadMSH(msh);
    GenerateBlock(bx, by, n);
    IdentifyGrains(scale);
    Write_HDF5(outputFile);
}


void GrainProcessor2D::LoadMSH(std::string fileName)
{
    std::vector<Eigen::Vector2f> vertices;
    std::vector<std::array<int,4>> elems;
    std::vector<Triangle> tris1;

    spdlog::info("load {}", fileName);
    if(!std::filesystem::exists(fileName)) spdlog::critical("file does not exist");

    gmsh::clear();
    gmsh::option::setNumber("General.Terminal", 0);
    gmsh::open(fileName);

    // get nodes
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords, parametricCoords;
    std::unordered_map<std::size_t, int> mtags; // gmsh nodeTag -> node sequential number

    gmsh::model::mesh::getNodesByElementType(2, nodeTags, nodeCoords, parametricCoords);
    vertices.reserve(nodeTags.size()*9);

    // set the size of the resulting nodes array
    for(unsigned i=0;i<nodeTags.size();i++)
    {
        std::size_t tag = nodeTags[i];
        if(mtags.count(tag)>0) continue; // throw std::runtime_error("GetFromGmsh() node duplication in deformable");
        Eigen::Vector2f coords = Eigen::Vector2f(nodeCoords[i*3+0], nodeCoords[i*3+1]);
        mtags[tag] = vertices.size();
        vertices.push_back(coords);
    }
    spdlog::info("nodes read {}",vertices.size());
    int nVerticesOriginal = vertices.size();
    vertices2.resize(nVerticesOriginal*9);

    // create surrounding copies
    int count = 0;
    for(int i=-1;i<=1;i++)
        for(int j=-1;j<=1;j++)
        {
            for(int m=0;m<nVerticesOriginal;m++)
                vertices2[m + nVerticesOriginal*count] = vertices[m] + Eigen::Vector2f(i,j);
            count++;
        }

    std::vector<std::pair<int,int>> dimTagsGrains;
    gmsh::model::getEntities(dimTagsGrains,2);
    spdlog::info("dimTagsGrains size {}",dimTagsGrains.size());

    //    std::unordered_set<int> used_nodes;
    elems.reserve(dimTagsGrains.size()*50);
    tris1.clear();
    tris1.reserve(elems.size());

    std::vector<int> types;
    std::vector<std::vector<std::size_t>> elemtags, nodetags;
    gmsh::model::mesh::getElements(types, elemtags, nodetags);

    for(int i=0;i<types.size();i++)
        spdlog::info("type {}; elemtags {}", types[i], elemtags[i].size());


    for(std::size_t j=0;j<dimTagsGrains.size();j++)
    {
        std::vector<std::size_t> tetraTags, nodeTagsInTetra;
        int entityTag = dimTagsGrains[j].second;
        gmsh::model::mesh::getElementsByType(2, tetraTags, nodeTagsInTetra, entityTag);
//        spdlog::info("grain {}; tetratags {}", j, tetraTags.size());

        for(std::size_t i=0;i<tetraTags.size();i++)
        {
            std::array<int,4> elem;
            Triangle t;

            for(int k=0;k<3;k++)
            {
                elem[k] = mtags.at(nodeTagsInTetra[i*3+k]);
                t.nds[k] = vertices[elem[k]];
            }
            t.grain = elem[3] = j;
            elems.push_back(elem);
            tris1.push_back(t);
        }
    }
    spdlog::info("tris read {}", tris1.size());

    int nTetraOriginal = tris1.size();
    elems2.resize(nTetraOriginal*9);
    tris2.resize(nTetraOriginal*9);

    count = 0;
    for(int i=-1;i<=1;i++)
        for(int j=-1;j<=1;j++)
            {
                for(int m=0;m<tris1.size();m++)
                {
                    int idx = m + count*nTetraOriginal;
                    elems2[idx] = elems[m];
                    Triangle &t = tris2[idx];
                    for(int n=0;n<3;n++)
                    {
                        elems2[idx][n] += count*nVerticesOriginal;
                        t.nds[n] = vertices2[elems2[idx][n]];
                    }
                    t.grain = elems2[m][3];
                }
                count++;
            }

    gmsh::clear();
}

void GrainProcessor2D::Write_HDF5(std::string fileName)
{
    H5::H5File file(fileName, H5F_ACC_TRUNC);

    hsize_t dims_grains[1] = {grainID.size()};
    H5::DataSpace dataspace_points_grains(1, dims_grains);
    H5::DataSet dataset_grainids = file.createDataSet("GrainIDs", H5::PredType::NATIVE_INT16, dataspace_points_grains);
    dataset_grainids.write(grainID.data(), H5::PredType::NATIVE_INT16);

    hsize_t dims_points[2] = {buffer.size(), 2};
    H5::DataSpace dataspace_points(2, dims_points);
    hsize_t chunk_dims[2] = {1024*1024,2};
    if(chunk_dims[0] > buffer.size()) chunk_dims[0] = std::max(buffer.size()/10,(size_t)1);
    spdlog::info("chunk {}; dims_pts {}", chunk_dims[0], dims_points[0]);
    H5::DSetCreatPropList proplist;
    proplist.setChunk(2, chunk_dims);
    proplist.setDeflate(7);
    H5::DataSet dataset_points = file.createDataSet("Points_Raw_2D", H5::PredType::NATIVE_FLOAT, dataspace_points, proplist);
    dataset_points.write(buffer.data(), H5::PredType::NATIVE_FLOAT);

    // write volume as attribute
    hsize_t att_dim = 1;
    H5::DataSpace att_dspace(1, &att_dim);
    H5::Attribute att_volume = dataset_grainids.createAttribute("volume", H5::PredType::NATIVE_FLOAT, att_dspace);
    att_volume.write(H5::PredType::NATIVE_FLOAT, &volume);
    spdlog::info("volume written {}",volume);
    file.close();
}

bool GrainProcessor2D::PointInsideTriangle(Eigen::Vector2f point, Eigen::Vector2f t[3])
{
    constexpr float eps = 1e-3;
    Eigen::Matrix2f M,B;
    M << t[1]-t[0], t[2]-t[0];
    B = M.inverse();
    point -= t[0];
    Eigen::Vector2f b = B * point; // barycentric
    return (b[0]>-eps && b[1]>-eps && (b.sum() < 1+eps));
}

void GrainProcessor2D::GenerateBlock(float dx, float dy, int n)
{
    constexpr float magic_constant = 0.656;
    volume = dx*dy;

    const float kRadius = sqrt(magic_constant*volume/n);

    const std::array<float, 2>kXMin{0, 0};
    const std::array<float, 2>kXMax{dx, dy};
    buffer = thinks::PoissonDiskSampling(kRadius, kXMin, kXMax);
}

void GrainProcessor2D::IdentifyGrains(const float scale)
{
    spdlog::info("identify grains");

    BVHN2D::BVHNFactory.release(leaves);
    leaves.reserve(tris2.size());

    for(int i=0;i<tris2.size();i++)
    {
        Triangle &t = tris2[i];
        BVHN2D* bvhn = BVHN2D::BVHNFactory.take();
        bvhn->isLeaf = true;
        bvhn->box.Reset();
        for(int j=0;j<3;j++) bvhn->Expand(t.nds[j]);
        bvhn->child1 = bvhn->child2 = nullptr;
        bvhn->elem = i;
        leaves.push_back(bvhn);
    }

    spdlog::info("building bvh");
    root.Build(leaves, 0);
    spdlog::info("finished building bvh");

    grainID.resize(buffer.size());
    spdlog::info("tetra2 {}; grainID {}", tris2.size(), grainID.size());

    // identify grains
    int unidentifiedPoints = 0;
    int problematicPoints = 0;

#pragma omp parallel for reduction(+:unidentifiedPoints,problematicPoints)
    for(int i=0;i<buffer.size();i++)
    {
        std::array<float, 2> &arr = buffer[i];
        Eigen::Vector2f v(arr[0],arr[1]);
        v *= scale;
//        for(int j=0;j<2;j++) v[j] = v[j] - floor(v[j]);
        v -= v.array().floor().matrix();

        BVHN2D bvhn;
        bvhn.isLeaf = true;
        bvhn.elem = -1;
        bvhn.box.Reset();
        bvhn.Expand(v);

        std::vector<std::pair<BVHN2D*,BVHN2D*>> broad_list;
        broad_list.reserve(10);
        bvhn.Collide(&root, broad_list);

        grainID[i] = -1;
        for(auto &pair : broad_list)
        {
            BVHN2D *bvhn2 = pair.second;
            int idx = bvhn2->elem;
            if(idx == -1) spdlog::critical("elem index -1");
            Triangle &t = tris2[idx];
            bool result = PointInsideTriangle(v, t.nds);
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
                short grain = (short)tris2[p.second->elem].grain;
                grainID[i] = grain;
                if(grain < 0) spdlog::critical("elem {} grain {}",p.second->elem,grain);
            }
            else
            {
                spdlog::critical("broad list 0");
            }
        }
        if(grainID[i] == -1) { grainID[i] = 0; unidentifiedPoints++; }
    }
    spdlog::info("finished processing points; problematic {}; unidentified {}",problematicPoints, unidentifiedPoints);
}
