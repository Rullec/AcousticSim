#pragma once
#include "sim/gpu_utils/CudaDef.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <type_traits>
// namespace HandyMatrix
// {
template <typename dtype, int N, int M> class tCudaMatrix
{
public:
    static constexpr int mRandSeed = 0;
    static constexpr int mRows = N;
    static constexpr int mCols = M;
    static constexpr int mElements = N * M;
    dtype mData[mElements];
    SIM_CUDA_CALLABLE explicit tCudaMatrix() { setZero(); }

    SIM_CUDA_CALLABLE explicit tCudaMatrix(
        const std::initializer_list<dtype> &list)
    {
        assert(list.size() == N * M);
        int i = 0;
        for (auto &res : list)
            mData[i++] = res;
        // for (int i = 0; i < N * M; i++)
    }
    // SIM_CUDA_CALLABLE operator dtype() const
    // {
    //     static_assert(N == 1 && M == 1);
    //     return mData[0];
    // }
    SIM_CUDA_CALLABLE size_t size() const { return mElements; }
    SIM_CUDA_CALLABLE void setValue(const dtype &val)
    {
        // std::cout << "setvalue begin, val = " << val << "bytes = " <<
        // sizeof(dtype) * mElements << std::endl;
        for (int i = 0; i < mElements; i++)
            mData[i] = val;
    }
    SIM_CUDA_CALLABLE void setRandom()
    {
        for (int i = 0; i < mElements; i++)
        {
            if (std::is_integral_v<dtype>)
            {
                mData[i] = static_cast<dtype>(rand());
            }
            else
            {
                // is floating point_v
                mData[i] = static_cast<dtype>(rand() * 1.0 / RAND_MAX);
            }
        }
    }
    SIM_CUDA_CALLABLE void setZero()
    {
        // printf("---begin to set zero---\n");
        std::memset(mData, 0, sizeof(dtype) * mElements);
        // for (int i = 0; i < mElements; i++)
        //     mData[i] = 0;
        // printf("---end to set zero---\n");
    }
    SIM_CUDA_CALLABLE void setIdentity()
    {
        // printf("------begin to set I------\n");
        static_assert(N == M);
        this->setZero();
        for (int i = 0; i < std::min(N, M); i++)
        {
            // printf("--!!!set %d %d = 1--\n", i, i);
            (*this)(i, i) = 1;
        }
        // printf("------end to set I------\n");
    }
    SIM_CUDA_CALLABLE tCudaMatrix(const dtype &val)
    {
        // std::cout << "begin to set val " << val << std::endl;
        setValue(val);
        // // std::cout << "after to set val " << mData[0] << std::endl;
    }
    static SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, M> Zero()
    {
        tCudaMatrix<dtype, N, M> mat;
        return mat;
    }
    static SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, M> Ones()
    {
        tCudaMatrix<dtype, N, M> mat;
        mat.setValue(1);
        return mat;
    }
    static SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, M> Identity()
    {
        tCudaMatrix<dtype, N, M> mat;
        mat.setIdentity();
        return mat;
    }
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, M> cwiseAbs()
    {
        tCudaMatrix<dtype, N, M> data(*this);
        for (int i = 0; i < mElements; i++)
            data.mData[i] = std::fabs(data.mData[i]);
        return data;
    }
    SIM_CUDA_CALLABLE dtype maxCoeff()
    {
        dtype cur_max = -1e16;
        for (int i = 0; i < mElements; i++)
        {
            cur_max = SIM_MAX(cur_max, mData[i]);
        }
        return cur_max;
    }
    SIM_CUDA_CALLABLE dtype minCoeff()
    {
        dtype cur_min = 1e16;
        for (int i = 0; i < mElements; i++)
        {
            cur_min = SIM_NIN(cur_min, mData[i]);
        }
        return cur_min;
    }
    SIM_CUDA_CALLABLE dtype operator()(int i, int j) const
    {

        if (false == (i >= 0 && i < N))
        {
            printf("visit (%d, %d) in mat size(%d, %d)\n", i, j, N, M);
            assert(false);
        }
        if (false == (j >= 0 && j < M))
        {
            printf("visit (%d, %d) in mat size(%d, %d)\n", i, j, N, M);
            assert(false);
        }
        return mData[j * mRows + i];
    }
    SIM_CUDA_CALLABLE dtype &operator()(int i, int j)
    {
        assert(i >= 0 && i < N);
        assert(j >= 0 && j < M);
        return mData[j * mRows + i];
    }
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, M, N> transpose() const
    {
        tCudaMatrix<dtype, M, N> res;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++)
            {
                res(j, i) = (*this)(i, j);
            }
        return res;
    }
    SIM_CUDA_CALLABLE dtype operator[](int idx) const { return mData[idx]; }
    SIM_CUDA_CALLABLE dtype &operator[](int idx) { return mData[idx]; }
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, 1> col(int idx) const
    {
        tCudaMatrix<dtype, N, 1> res;
        for (int i = 0; i < N; i++)
        {
            res[i] = (*this)(i, idx);
        }
        return res;
    }

    SIM_CUDA_CALLABLE tCudaMatrix<dtype, 1, M> row(int idx) const
    {
        tCudaMatrix<dtype, 1, M> res;
        for (int i = 0; i < M; i++)
        {
            res[i] = (*this)(idx, i);
        }
        return res;
    }
    SIM_CUDA_CALLABLE void setcol(const tCudaMatrix<dtype, N, 1> &target_col,
                                  int col_id)
    {
        for (int i = 0; i < N; i++)
        {
            mData[col_id * mRows + i] = target_col[i];
        }
    }
    SIM_CUDA_CALLABLE void setrow(const tCudaMatrix<dtype, 1, M> &target_col,
                                  int row_id)
    {
        for (int col_id = 0; col_id < M; col_id++)
        {
            mData[col_id * mRows + row_id] = target_col[i];
        }
    }
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, M> operator*(const dtype &val) const
    {
        tCudaMatrix<dtype, N, M> res = (*this);
        for (int i = 0; i < mElements; i++)
            res[i] *= val;
        return res;
    }
    SIM_CUDA_CALLABLE void operator*=(const dtype &val)
    {
        for (int i = 0; i < mElements; i++)
            mData[i] *= val;
    }
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, M> operator/(const dtype &val) const
    {
        tCudaMatrix<dtype, N, M> res = (*this);
        for (int i = 0; i < mElements; i++)
            res[i] /= val;
        return res;
    }
    // pre postioning negative
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, M> operator-()
    {
        tCudaMatrix<dtype, N, M> negative = (*this);
        for (int i = 0; i < mElements; i++)
        {
            negative.mData[i] *= -1;
        }
        return negative;
    }
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, M>
    operator-(const tCudaMatrix<dtype, N, M> &another_mat) const
    {
        tCudaMatrix<dtype, N, M> negative = (*this);
        for (int i = 0; i < mElements; i++)
        {
            negative.mData[i] -= another_mat.mData[i];
        }
        return negative;
    }
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, M>
    operator+(const tCudaMatrix<dtype, N, M> &another_mat) const
    {
        tCudaMatrix<dtype, N, M> negative = (*this);
        for (int i = 0; i < mElements; i++)
        {
            negative.mData[i] += another_mat.mData[i];
        }
        return negative;
    }
    SIM_CUDA_CALLABLE void
    operator+=(const tCudaMatrix<dtype, N, M> &another_mat)
    {

        for (int i = 0; i < mElements; i++)
        {
            mData[i] += another_mat.mData[i];
        }
    }
    SIM_CUDA_CALLABLE void
    operator-=(const tCudaMatrix<dtype, N, M> &another_mat)
    {

        for (int i = 0; i < mElements; i++)
        {
            mData[i] -= another_mat.mData[i];
        }
    }

    template <int K>
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, K>
    operator*(const tCudaMatrix<dtype, M, K> &mat) const
    {
        // printf("N %d M %d K %d\n", N, M, K);
        tCudaMatrix<dtype, N, K> res;
        for (int i = 0; i < N; i++)
            for (int j = 0; j < K; j++)
            {
                // printf("row %d and col %d\n",i, j);
                res(i, j) = (*this).row(i).dot((mat.col(j)));
            }
        return res;
    }
    friend SIM_CUDA_CALLABLE std::ostream &
    operator<<(std::ostream &fout, const tCudaMatrix<dtype, N, M> &mat)
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                fout << mat(i, j);
                if (j != M - 1)
                    fout << " ";
            }
            if (i != N - 1)
                fout << std::endl;
        }
        return fout;
    };
    SIM_CUDA_CALLABLE void print() const
    {
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                printf("%.5f", (*this)(i, j));
                if (j != M - 1)
                    printf(" ");
            }
            if (i != N - 1)
                printf("\n");
        }
    }
    friend SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, M>
    operator*(dtype val, const tCudaMatrix<dtype, N, M> &mat)
    {
        tCudaMatrix<dtype, N, M> new_mat(mat);
        for (int i = 0; i < new_mat.mElements; i++)
            new_mat.mData[i] *= val;
        return new_mat;
    };

    SIM_CUDA_CALLABLE dtype norm() const
    {
        dtype sum = 0;
        for (int i = 0; i < N * M; i++)
            sum += mData[i] * mData[i];
        return std::sqrt(sum);
    }

    SIM_CUDA_CALLABLE void normalize()
    {
        dtype cur_norm = norm();
        for (int i = 0; i < mElements; i++)
            mData[i] /= cur_norm;
        // val ;
    }
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, M> normalized() const
    {
        tCudaMatrix<dtype, N, M> new_value = (*this);
        new_value.normalize();
        return new_value;
    }
    template <int block_row, int block_col>
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, block_row, block_col>
    block(int st_row, int st_col) const
    {
        tCudaMatrix<dtype, block_row, block_col> mat;
        // printf("visit blocksize [%d, %d] st %d %d total size %d %d\n",
        // block_row, block_col, st_row, st_col, N, M);
        for (int i = st_row; i < st_row + block_row; i++)
            for (int j = st_col; j < st_col + block_col; j++)
            {
                mat(i - st_row, j - st_col) = (*this)(i, j);
            }
        return mat;
    }
    template <int K, int L>
    SIM_CUDA_CALLABLE void setBlock(int row_st, int col_st,
                                    const tCudaMatrix<dtype, K, L> &new_block)
    {
        assert(row_st + K <= N);
        assert(col_st + L <= M);
        for (int i = 0; i < K; i++)
            for (int j = 0; j < L; j++)
            {
                (*this)(row_st + i, col_st + j) = new_block(i, j);
            }
    }

    SIM_CUDA_CALLABLE tCudaMatrix<dtype, N, N> inverse() const
    {
        tCudaMatrix<dtype, N, N> res(*this);
        if (N == 2 && M == 2)
        {
            dtype a = res[0], c = res[1], b = res[2], d = res[3];
            dtype deno = 1.0 / (a * d - b * c);
            res[0] = d * deno;
            res[1] = -c * deno;
            res[2] = -b * deno;
            res[3] = a * deno;
        }
        else if (N == 3 && M == 3)
        {
            dtype a = res(0, 0), b = res(0, 1), c = res(0, 2), d = res(1, 0),
                  e = res(1, 1), f = res(1, 2), g = res(2, 0), h = res(2, 1),
                  i = res(2, 2);
            res(0, 0) = e * i - f * h;
            res(0, 1) = c * h - b * i;
            res(0, 2) = b * f - c * e;
            res(1, 0) = f * g - d * i;
            res(1, 1) = a * i - c * g;
            res(1, 2) = c * d - a * f;
            res(2, 0) = d * h - e * g;
            res(2, 1) = b * g - a * h;
            res(2, 2) = a * e - b * d;
            float prefix = 1.0 / (a * (e * i - f * h) - b * (d * i - f * g) +
                                  c * (d * h - e * g));
            res *= prefix;
        }
        else
        {
            assert(false && "unsupported inverse");
            exit(1);
        }
        return res;
    }

    template <int A, int B>
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, A * N, B * M>
    KroneckerProduct(const tCudaMatrix<dtype, A, B> &mat)
    {
        tCudaMatrix<dtype, A * N, B * M> result;
        // result.setBlock(0, 0, mat);
        for (int i = 0; i < N; i++)
            for (int j = 0; j < M; j++)
            {
                // printf("set from (%d, %d) block size (%d %d)\n", i * A, j *
                // B,
                //    A, B);
                result.setBlock(i * A, j * B, mat * (*this)(i, j));
            }
        return result;
    }
#if defined(__CUDACC__)
    __device__ bool hasNan() const
    {
        for (int i = 0; i < mElements; i++)
        {
            if (isnan(mData[i]))
                return true;
        }
        return false;
    }
#else

    SIM_CUDA_CALLABLE bool hasNan() const
    {
        for (int i = 0; i < mElements; i++)
        {
            if (cCudaMath::IsNan(mData[i]))
                return true;
        }
        return false;
    }
#endif
    // ------------vector method--------------
    template <int K>
    SIM_CUDA_CALLABLE dtype dot(const tCudaMatrix<dtype, K, 1> &vec) const
    {
        static_assert(N == 1 || M == 1);
        assert(N == K || M == K);
        dtype sum = 0;
        for (int i = 0; i < K; i++)
        {
            sum += vec[i] * mData[i];
        }
        return sum;
    }
    template <int K>
    SIM_CUDA_CALLABLE tCudaMatrix<dtype, K, 1> segment(int st) const
    {
        assert(K <= N * M);
        static_assert(N == 1 || M == 1); // a vecotr method
        tCudaMatrix<dtype, K, 1> res;
        for (int i = 0; i < K; i++)
            res[i] = mData[st + i];
        return res;
    }

    template <int K>
    SIM_CUDA_CALLABLE void setSegment(int st,
                                      const tCudaMatrix<dtype, K, 1> &new_seg)
    {
        assert(st + K < N * M);
        for (int i = 0; i < K; i++)
            mData[st + i] = new_seg[i];
    }
};

/*
1. define square matrix: Matrix2i, Matrix2f, Natrix3i, Matrix3f
2. define vector: Vector2i, Vector2f, Vector3f, Vector3f
3. define row vecotr: RowVector2i, RowVector2f...

So we need an macro to define them all: Given
1. data type (float), data type suffix (f)
2. size
*/

#define CUDA_DECLARE_MATRIX(type, type_suffix, rows, cols)                     \
    using tCudaMatrix##rows##cols##type_suffix = tCudaMatrix<type, rows, cols>;

#define CUDA_DECLARE_SQR_MATRIX(type, type_suffix, size)                       \
    using tCudaMatrix##size##type_suffix = tCudaMatrix<type, size, size>;

#define CUDA_DECLARE_VECTOR(type, type_suffix, size)                           \
    using tCudaVector##size##type_suffix = tCudaMatrix<type, size, 1>;

#define CUDA_DECLARE_MATRIX_AND_VECTOR(type, type_suffix, size)                \
    CUDA_DECLARE_VECTOR(type, type_suffix, size)                               \
    CUDA_DECLARE_SQR_MATRIX(type, type_suffix, size)

CUDA_DECLARE_MATRIX_AND_VECTOR(int, i, 2);
CUDA_DECLARE_MATRIX_AND_VECTOR(int, i, 3);
CUDA_DECLARE_MATRIX_AND_VECTOR(int, i, 4);
CUDA_DECLARE_MATRIX_AND_VECTOR(float, f, 2);
CUDA_DECLARE_MATRIX_AND_VECTOR(float, f, 3);
CUDA_DECLARE_MATRIX_AND_VECTOR(float, f, 4);
CUDA_DECLARE_MATRIX_AND_VECTOR(float, f, 9);
CUDA_DECLARE_MATRIX_AND_VECTOR(float, f, 12);
// CUDA_DECLARE_MATRIX_AND_VECTOR(float, f, 9);
// CUDA_DECLARE_VECTOR(int, i, 12);
CUDA_DECLARE_VECTOR(int, i, 32);
CUDA_DECLARE_VECTOR(int, i, 8);
CUDA_DECLARE_MATRIX(float, f, 3, 2);
CUDA_DECLARE_MATRIX(float, f, 9, 2);
using tCudaMatrix9f = tCudaMatrix<float, 9, 9>;
using tCudaMatrix3f = tCudaMatrix<float, 3, 3>;

// typedef tCudaMatrix<float, 1, 3> tCudaRowVector3f;
// typedef tCudaMatrix<float, 3, 3> tCudaMatrix3f;
// typedef tCudaMatrix<float, 1, 2> tCudaRowVector2f;
// typedef tCudaMatrix<float, 2, 2> tCudaMatrix2f;
// typedef tCudaMatrix<float, 9, 9> tCudaMatrix9f;
// typedef tCudaMatrix<float, 3, 2> tCudaMatrix32f;

// typedef tCudaMatrix<float, 3, 1> tCudaVector3f;
// typedef tCudaMatrix<float, 2, 1> tCudaVector2f;
