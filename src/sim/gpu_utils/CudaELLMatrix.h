#include "sim/gpu_utils/CudaArray.h"
#include "sim/gpu_utils/CudaMatrix.h"
#include <vector>

/*
    3 * 3 based ELL matrix: each element is a 3x3 matrix
*/
class cCudaELLMatrix
{
public:
    using tELLRow = std::vector<tCudaMatrix3f>;
    using tELLColumnIdPerRow = std::vector<int>;
    explicit cCudaELLMatrix();
    void Init(int rows_3block, int cols_3block,
              const std::vector<tELLColumnIdPerRow> &column_idx_per_row);
    void AddValue(int row, int col, const tCudaMatrix3f &result);
    void SetZero();
protected:
    int mRows, mCols;

    std::vector<tELLRow> mRowData; // data in each row

    std::vector<std::vector<int>> mColumnIdLst; //
};

void ProdELLMatrixVector(const cCudaELLMatrix & A, const cCudaArray<tCudaMatrix3f> & x, cCudaArray<tCudaMatrix3f> & result);