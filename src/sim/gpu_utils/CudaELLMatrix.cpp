#include "CudaELLMatrix.h"
#include <algorithm>
#include <cassert>
cCudaELLMatrix::cCudaELLMatrix()
{
    mRows = 0;
    mCols = 0;
    mRowData.clear();
    mColumnIdLst.clear();
}
void cCudaELLMatrix::Init(int rows, int cols,
                          const std::vector<tELLColumnIdPerRow> &column_id_lst)
{
    mRows = rows, mCols = cols;

    mRowData.resize(rows);
    assert(column_id_lst.size() == mRows);
    mColumnIdLst = column_id_lst;
    // sort
    for (int i = 0; i < mRows; i++)
    {
        auto &column_id_per_row = mColumnIdLst[i];
        int row_data_size = mColumnIdLst[i].size();
        mRowData[i].resize(row_data_size);

        std::sort(column_id_per_row.begin(), column_id_per_row.end());
    }
}

void cCudaELLMatrix::AddValue(int row, int col, const tCudaMatrix3f &result)
{
    auto &cur_row_data = mRowData[row];
    auto &column_id_lst = mColumnIdLst[row];
    
}

void cCudaELLMatrix::SetZero()
{
    for(auto & row : mRowData)
    {
        for(auto & col : row)
        {
            col.setZero();
        }
    }
}