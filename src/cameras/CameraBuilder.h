#pragma once
#include "CameraBase.h"
#include "utils/DefUtil.h"
#include "utils/JsonUtil.h"
#include <memory>

class cCameraBuilder
{
public:
    static CameraBase *BuildCamera(const Json::Value &camera_json);
};