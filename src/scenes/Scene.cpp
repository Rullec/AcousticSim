#include "Scene.h"
cScene::cScene() { Reset(); }
cScene::~cScene() {}

void cScene::Update(double dt)
{
    mCurdt = dt;
    mCurTime += dt;
}

void cScene::Reset() { mCurTime = 0; }