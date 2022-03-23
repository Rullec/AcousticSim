#include "RenderUtil.h"
#include "geometries/Primitives.h"
#include "utils/ColorUtil.h"

void cRenderUtil::CalcTriangleDrawBufferSingle(tVertexPtr v0, tVertexPtr v1,
                                               tVertexPtr v2,

                                               const tVector &color,
                                               Eigen::Map<tVectorXf> &buffer,
                                               int &st_pos)
{
    buffer.segment(st_pos, 3) = v0->mPos.segment(0, 3).cast<float>();
    buffer.segment(st_pos + 3, 4) = color.cast<float>();
    // buffer[st_pos + 6] = 0.5;
    buffer.segment(st_pos + 7, 3) = v0->mNormal.segment(0, 3).cast<float>();

    st_pos += RENDERING_SIZE_PER_VERTICE;
    buffer.segment(st_pos, 3) = v1->mPos.segment(0, 3).cast<float>();
    buffer.segment(st_pos + 3, 4) = color.cast<float>();
    // buffer[st_pos + 6] = 0.5;
    buffer.segment(st_pos + 7, 3) = v1->mNormal.segment(0, 3).cast<float>();

    st_pos += RENDERING_SIZE_PER_VERTICE;
    buffer.segment(st_pos, 3) = v2->mPos.segment(0, 3).cast<float>();
    buffer.segment(st_pos + 3, 4) = color.cast<float>();
    buffer.segment(st_pos + 7, 3) = v2->mNormal.segment(0, 3).cast<float>();
    // buffer[st_pos + 6] = 0.5;
    st_pos += RENDERING_SIZE_PER_VERTICE;
}

void cRenderUtil::CalcEdgeDrawBufferSingle(tVertexPtr v0, tVertexPtr v1,
                                           const tVector &edge_normal,
                                           Eigen::Map<tVectorXf> &buffer,
                                           int &st_pos, const tVector &color)
{

    tVector3f bias_amp =
        1e-4f * edge_normal.cast<float>().segment(0, 3); // 0.1 mm

    // pos, color, normal
    buffer.segment(st_pos, 3) = v0->mPos.segment(0, 3).cast<float>() + bias_amp;
    buffer.segment(st_pos + 3, 4) = color.cast<float>();
    buffer.segment(st_pos + 7, 3) = tVector3f(0, 0, 0);
    st_pos += RENDERING_SIZE_PER_VERTICE;

    buffer.segment(st_pos, 3) = v1->mPos.segment(0, 3).cast<float>() + bias_amp;
    buffer.segment(st_pos + 3, 4) = color.cast<float>();
    buffer.segment(st_pos + 7, 3) = tVector3f(0, 0, 0);
    st_pos += RENDERING_SIZE_PER_VERTICE;
}

void cRenderUtil::CalcPointDrawBufferSingle(const tVector &v_pos,
                                            const tVector &v_color,
                                            Eigen::Map<tVectorXf> &buffer,
                                            int &st_pos)
{
    buffer.segment(st_pos, 3) = v_pos.segment(0, 3).cast<float>();
    buffer.segment(st_pos + 3, 4) = v_color.cast<float>();
    st_pos += RENDERING_SIZE_PER_VERTICE;
}