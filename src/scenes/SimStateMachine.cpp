#include "SimStateMachine.h"
#include "GLFW/glfw3.h"

static const std::string gSimState[eSimState::NUM_OF_SIMSTATE] = {
    "BeforeInit",
    "Pause",
    "Run",
    "Run(single_frame)"};
std::string tSimStateMachine::BuildStateStr(eSimState state)
{
    return gSimState[state];
}

tSimStateMachine::tSimStateMachine()
{
    mCurState = SIMSTATE_BEFORE_INIT;
    // if (start_run)
    //     mCurState = eSimState::SIMSTATE_RUN;
    // else
    //     mCurState = eSimState::SIMSTATE_PAUSE;
}

void tSimStateMachine::SimulatorInitDone(eSimState start_state)
{
    mCurState = start_state;
}

void tSimStateMachine::SimulationForwardOneFrameDone()
{
    if (mCurState == eSimState::SIMSTATE_RUN_SINGLE_FRAME)
    {
        mCurState = eSimState::SIMSTATE_PAUSE;
    }
}
eSimState tSimStateMachine::GetCurrentState() const
{
    return mCurState;
}
void tSimStateMachine::Key(int key, int scancode, int action, int mods)
{
    if (action == GLFW_PRESS)
    {
        switch (key)
        {
        case GLFW_KEY_O:
            mCurState = eSimState::SIMSTATE_RUN_SINGLE_FRAME;
            break;
        case GLFW_KEY_I:
        {
            if (mCurState == eSimState::SIMSTATE_PAUSE)
                mCurState = eSimState::SIMSTATE_RUN;
            else if (mCurState == eSimState::SIMSTATE_RUN)
                mCurState = eSimState::SIMSTATE_PAUSE;
            break;
        }
        default:
            break;
        }
    }
}
bool tSimStateMachine::IsRunning()
{
    return mCurState == SIMSTATE_RUN || mCurState == SIMSTATE_RUN_SINGLE_FRAME;
}