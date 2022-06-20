#include "ViscoMassSpring.h"
#include "geometries/Primitives.h"
#include "imgui.h"
#include "utils/ColorUtil.h"
#include "utils/ObjUtil.h"
cViscoMassSpring::cViscoMassSpring(int id_)
    : cBaseObject(eObjectType::VISCOSITY_MASS_SPRING_TYPE, id_)
{
}

cViscoMassSpring::~cViscoMassSpring() {}

void cViscoMassSpring::Init(const Json::Value &conf)
{
    mInitProp.mCurX = 1;
    mInitProp.mCurV = 0; // pos and vel
    mMass = 1;
    mK = 1;
    mEta = 1; // spring stiffness and viscosity coef
    mDt = 1e-2;
    // load a ball
    cObjUtil::LoadObj("data/sphere.obj", this->mVertexArray, this->mEdgeArray,
                      this->mTriangleArray);
    for (auto &x : mVertexArray)
    {
        x->mPos *= 0.1;
        x->mPos[1] += 0.1;
    }
    for (auto &t : mTriangleArray)
        t->mColor = ColorBlue;
    for (auto &t : mEdgeArray)
        t->mColor = ColorBlack;
    mInitPos.resize(this->mVertexArray.size());

    CalcTriangleInitArea();
    UpdateTriangleNormal();
    UpdateVertexNormalFromTriangleNormal();

    for (int i = 0; i < mInitPos.size(); i++)
    {
        mInitPos[i] = mVertexArray[i]->mPos;
    }
    Reset();
}

void cViscoMassSpring::Reset()
{
    mCurX = mInitProp.mCurX;
    mCurV = mInitProp.mCurV; // pos and vel
    // mMass = mInitProp.mMass;
    // mK = mInitProp.mK;
    // mEta = mInitProp.mEta; // spring stiffness and viscosity coef
    UpdateMesh();
}
/*

v_next = (dt * (fext - k * x) + m * v) / (m + dt * (k * dt + eta))
*/
void cViscoMassSpring::Update(float dt)
{
    dt = mDt;
    // 1. get ext force
    double fext = 0;
    double v_next = (dt * (fext - mK * mCurX) + mMass * mCurV) /
                    (mMass + dt * (mK * dt + mEta));
    double x_next = dt * v_next + mCurX;

    // 2. calculate next pos
    mCurV = v_next;
    mCurX = x_next;
    printf("cur x %.3f cur v %.3f\n", mCurX, mCurV);

    // 3. update buf
    UpdateMesh();
}

void cViscoMassSpring::UpdateImGui()
{
    ImGui::DragFloat("k", &mK, 1.0, 0.0, 100);
    ImGui::DragFloat("eta", &mEta, 1, 0, 100);
    ImGui::DragFloat("mass", &mMass, 1, 1, 100);
    ImGui::DragFloat("dt", &mDt, 1e-3, 1e-3, 1e-1);
    if (ImGui::Button("reset"))
    {
        Reset();
    }
}

void cViscoMassSpring::UpdateMesh()
{

    for (int i = 0; i < mInitPos.size(); i++)
    {
        mVertexArray[i]->mPos = mInitPos[i] + tVector(mCurX, 0, 0, 0);
    }
}

void cViscoMassSpring::ApplyUserPerturbForceOnce(tPerturb *) { return; }