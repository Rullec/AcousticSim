#pragma once
#include "Primitives.h"
#include <string>
#include <vector>
/**
 * \brief           export object (with no texture and texture coords)
 */
// class cObjExporter
// {
// public:
//     static bool ExportObj();
//     // static bool ExportObj;
// };

// class cObjExporter
// {
// public:
//     static bool ExportObj(std::string export_path,
//                           const std::vector<tVertex *> &vertices_array,
//                           const std::vector<tTriangle *> &triangles_array);
// };
class cObjExporter
{
public:
    static bool ExportObj(std::string export_path,
                          const std::vector<tVertex *> &vertices_array,
                          const std::vector<tTriangle *> &triangles_array, bool silent = false);
};