#include "sim/AcousticSoftBody.h"
#include "sim/AudioOutput.h"
extern cAudioOutputPtr gAudioOutput;
cAcousticSoftBody::cAcousticSoftBody(int id_) : cSoftBodyImplicit(id_)
{
}
void cAcousticSoftBody::Init(const Json::Value &conf)
{
    cSoftBodyImplicit::Init(conf);

    // Init noise
    UpdateDeformationGradient();
    mGlobalStiffnessSparseMatrix = CalcGlobalStiffnessSparseMatrix();

    SolveVibration(mInvLumpedMassMatrixDiag.cwiseInverse(), mGlobalStiffnessSparseMatrix, GetRayleightDamping(), mXcur, mXprev);
    
}
void cAcousticSoftBody::Update(float dt)
{
    cSoftBodyImplicit::Update(dt);

    // update second
}

