#include "sim/softbody/FourOrderTensor.h"
#include <iostream>

void test_tensor()
{
    size_t i = 3, j = 3, k = 3, l = 3;
    cFourOrderTensor tensor(i, j, k, l);
    tMatrixXd res = tensor.ExpandToMatrix();
    res = tMatrixXd::Random(res.rows(), res.cols());
    tensor.LoadFromAMatrix(i,j,k,l,res);
    tMatrixXd new_res = tensor.ExpandToMatrix();
    std::cout << res << std::endl;
    std::cout << new_res << std::endl;
}