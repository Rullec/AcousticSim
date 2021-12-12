 #pragma once
#include <memory>
#include <string>
class cScene : public std::enable_shared_from_this<cScene>
{
public:
    explicit cScene();
    virtual ~cScene();
    virtual void Init(const std::string &conf_path) = 0;
    virtual void Update(double dt);
    virtual void Reset();

protected:
    double mCurdt;
    double mCurTime;
};
