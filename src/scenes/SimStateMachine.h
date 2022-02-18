#pragma once
#include <string>
enum eSimState
{
    SIMSTATE_BEFORE_INIT = 0,
    SIMSTATE_PAUSE,
    SIMSTATE_RUN,
    SIMSTATE_RUN_SINGLE_FRAME,
    NUM_OF_SIMSTATE
};

struct tSimStateMachine
{
public:
    explicit tSimStateMachine();
    void SimulatorInitDone(eSimState start_state);
    void SimulationForwardOneFrameDone();
    eSimState GetCurrentState() const;
    void Key(int key, int scancode, int action, int mods);
    static std::string BuildStateStr(eSimState state);
    bool IsRunning();

protected:
    eSimState mCurState;

    void KeyFromPause(int key, int scancode, int action, int mods);
    void KeyFromRun(int key, int scancode, int action, int mods);
    void KeyFromRunSingleFrame(int key, int scancode, int action, int mods);
};