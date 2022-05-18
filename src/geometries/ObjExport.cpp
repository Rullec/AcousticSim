#include "ObjExport.h"
#include "geometries/Primitives.h"
#include <fstream>
#include <iostream>
template <typename... Args>
inline std::string format_string(const char *format, Args... args)
{
    constexpr size_t oldlen = BUFSIZ;
    char buffer[oldlen]; // 默认栈上的缓冲区

    size_t newlen = snprintf(&buffer[0], oldlen, format, args...);
    newlen++; // 算上终止符'\0'

    if (newlen > oldlen)
    { // 默认缓冲区不够大，从堆上分配
        std::vector<char> newbuffer(newlen);
        snprintf(newbuffer.data(), newlen, format, args...);
        return std::string(newbuffer.data());
    }

    return buffer;
}
bool cObjExporter::ExportObj(std::string export_path,
                             const std::vector<tVertexPtr> &vertices_array,
                             const std::vector<tTrianglePtr> &triangles_array,
                             bool silent /* = false*/)
{
    // 1. output the vertices info
    std::ofstream fout(export_path, std::ios::out);
    for (int i = 0; i < vertices_array.size(); i++)
    {
        auto v = vertices_array[i];
        std::string cur_str = format_string("v %.5f %.5f %.5f\n", v->mPos[0],
                                            v->mPos[1], v->mPos[2]);
        fout << cur_str;
    }
    // if (enable_texutre_output == true)
    {
        // std::cout << "cloth texture coord *= 0.3\n";
        for (int i = 0; i < vertices_array.size(); i++)
        {
            auto v = vertices_array[i];
            std::string cur_str = format_string(
                "vt %.5f %.5f\n", v->muv[0] * 0.3, v->muv[1] * 0.3);
            fout << cur_str;
        }
    }

    // 2. output the face id
    // double thre = 1e-6;
    for (int i = 0; i < triangles_array.size(); i++)
    {
        auto t = triangles_array[i];
        std::string cur_str =
            format_string("f %d/%d %d/%d %d/%d\n", t->mId0 + 1, t->mId0 + 1,
                          t->mId1 + 1, t->mId1 + 1, t->mId2 + 1, t->mId2 + 1);
        // tVector pos0 = vertices_array[t->mId0]->mPos,
        //         pos1 = vertices_array[t->mId1]->mPos,
        //         pos2 = vertices_array[t->mId2]->mPos;

        // double diff0 = (pos0 - pos1).norm(), diff1 = (pos0 - pos2).norm(),
        //        diff2 = (pos1 - pos2).norm();
        // std::cout << " diff  - " << diff0 + diff1 + diff2 << std::endl;
        // if (diff0 < thre || diff1 < thre || diff2 < thre)
        // {
        //     std::cout << "pos0 = " << pos0.transpose() << std::endl;
        //     std::cout << "pos1 = " << pos1.transpose() << std::endl;
        //     std::cout << "pos2 = " << pos2.transpose() << std::endl;
        //     std::cout << "tri id = " << i << std::endl;
        //     exit(1);
        // }
        fout << cur_str;
    }
    if (silent == false)
        printf("[debug] export obj to %s\n", export_path.c_str());
    return true;
}

// bool cObjExporter::ExportObj(std::string export_path,
//                              const std::vector<tVertexPtr > &vertices_array,
//                              const std::vector<tTriangle *> &triangles_array)
// {

// bool cObjExporter::ExportObj(std::string export_path,
//                              const std::vector<tVertexPtr > &vertices_array,
//                              const std::vector<tTriangle *> &triangles_array)
// {

// }