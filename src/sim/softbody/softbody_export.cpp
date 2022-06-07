#include "sim/softbody/SoftBodyImplicit.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(softbody, m)
{
    py::class_<cSoftBodyImplicit, std::shared_ptr<cSoftBodyImplicit>>(
        m, "softbody")
        .def(py::init<>())
        .def("InitFromFile", &cSoftBodyImplicit::InitFromFile)
        .def("GetGlobalStiffnessMatrixEntrys",
             &cSoftBodyImplicit::GetGlobalStiffnessMatrixEntrys)
        .def("GetGlobalRawMassMatrixEntrys",
             &cSoftBodyImplicit::GetGlobalRawMassMatrixEntrys)
        .def("GetMassMatrixDiag", &cSoftBodyImplicit::GetMassMatrixDiag)
        .def("GetVertexPos", &cSoftBodyImplicit::GetVertexPos)
        .def("GetTetVertexIdx", &cSoftBodyImplicit::GetTetVertexIdx)
        .def("GetYoungsModulus", &cSoftBodyImplicit::GetYoungsModulus)
        .def("GetPoissonRatio", &cSoftBodyImplicit::GetPoissonRatio)
        .def("GetRayleighDamping", &cSoftBodyImplicit::GetRayleighDamping);
}