class cMPIUtil
{
public:
    static bool Init();
    static int GetWorldSize();
    static int GetWorldRank();
    static void End();
};