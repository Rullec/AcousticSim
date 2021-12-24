#include "TetUtil.h"
#include "utils/LogUtil.h"
#include "utils/FileUtil.h"
#include "utils/StringUtil.h"
#include "geometries/Primitives.h"
#include <iostream>

tVertexPtrVector LoadVertFile(const std::string &vert_file)
{
    SIM_INFO("load vertices from {}", vert_file);
    int num_of_vertices = -1;
    auto line_list = cFileUtil::ReadFileAllLines(vert_file);
    cStringUtil::RemoveEmptyLine(line_list);
    cStringUtil::RemoveCommentLine(line_list);
    int cur_line_idx = 0;

    tVertexPtrVector vert_vec = {};
    num_of_vertices = std::atoi(line.c_str());
    
    for (cur_line_idx; cur_line_idx < line_list.size(); cur_line_idx++)
    {
        auto line = line_list[cur_line_idx];
        if (line[0] == '#')
        {
            continue;
        }
        else
        {
            
            cur_line_idx++;
            break;
        }
    }
    for (cur_line_idx; cur_line_idx < line_list.size(); cur_line_idx++)
    {
        auto line = line_list[cur_line_idx];
        if (line[0] == '#')
        {
            continue;
        }
        else
        {
            // begin to parse
            auto splited = cStringUtil::SplitString(line, ",");
            std::vector<float> pos_list = {};
            for (auto &x : splited)
            {
                pos_list.push_back(std::atof(cStringUtil::Strip(x).c_str()));
            }
            SIM_ASSERT(pos_list.size() == 3);
            auto cur_vertex = std::make_shared<tVertex>();
            cur_vertex->mPos = tVector(pos_list[0], pos_list[1], pos_list[2], 1);
            vert_vec.push_back(cur_vertex);
        }
    }
    std::cout << "num of vertices = " << num_of_vertices << std::endl;
    return vert_vec;
}

tEigenArr<tVector4i> LoadTetFile(const std::string &vert_file)
{
    
}

/**
 * \brief               load tet mesh from ".node" files (custom format)
 * \param path          mesh file path
 * \param vertex_vec    all vertices
 * \param edge_vec      all edges
 * \param tri_vec       all triangles vector
 * \param tet_vec       all tets vector
*/
void cTetUtil::LoadTet(const std::string &path,
                       tVertexPtrVector &vertex_vec,
                       tEdgePtrVector &edge_vec,
                       tTrianglePtrVector &tri_vec,
                       tTetPtrVector &tet_vec)

{
    SIM_DEBUG("begin to load tet from {}", path);
    // SIM_ASSERT(cFileUtil::ExistsFile(path.c_str()));
    std::string vert_file = path + ".vert";
    std::string tet_file = path + ".tet";
    vertex_vec = LoadVertFile(vert_file);
    SIM_DEBUG("load {} vertices from file {}", vertex_vec.size(), vert_file);
    LoadTetFile(tet_file);

    return;
}