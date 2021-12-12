#pragma once
#include "spdlog/sinks/stdout_color_sinks.h"
#include <memory>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>
typedef std::shared_ptr<spdlog::logger> tLogger;
typedef spdlog::level::level_enum eLogLevel;

class cLogUtil
{
public:
    static void SetLoggingLevel(const std::string &level_name);
    static tLogger CreateLogger(const std::string &loggername);
    static void DropLogger(const std::string &loggername);
    static tLogger GetLogger(const std::string &loggername);
    static void Printf(const tLogger &logger, eLogLevel level, const char *fmt,
                       va_list args);

    // const static tLogger mGlobalLogger;
    inline const static tLogger mGlobalLogger =
        cLogUtil::CreateLogger("Global");

private:
    inline const static size_t buf_size = 1000;
    inline static char buf[buf_size];
    static int GetLevelEnumFromString(const std::string &level_name);
};

#if defined(_WIN32)
#define SIM_UNREACHABLE __assume(0);
#else
#define SIM_UNREACHABLE __builtin_unreachable();
#endif

#define __FILENAME__                                                           \
    (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define SPD_AUGMENTED_LOG(X, ...)                                              \
    cLogUtil::mGlobalLogger->X(                                                \
        fmt::format("[{}:{} @ {}] ", __FILENAME__, __LINE__, __FUNCTION__) +   \
        fmt::format(__VA_ARGS__))

#define SIM_OUTPUT(...) cLogUtil::mGlobalLogger->info(fmt::format(__VA_ARGS__))

#define SIM_TRACE(...) SPD_AUGMENTED_LOG(trace, __VA_ARGS__)
#define SIM_DEBUG(...) SPD_AUGMENTED_LOG(debug, __VA_ARGS__)
#define SIM_INFO(...) SPD_AUGMENTED_LOG(info, __VA_ARGS__)
#define SIM_WARN(...) SPD_AUGMENTED_LOG(warn, __VA_ARGS__)
#define SIM_ERROR(...)                                                         \
    {                                                                          \
        SPD_AUGMENTED_LOG(error, __VA_ARGS__);                                 \
        SIM_UNREACHABLE;                                                       \
    }

#define SIM_ASSERT_INFO(x, ...)                                                \
    {                                                                          \
        bool ___ret___ = static_cast<bool>(x);                                 \
        if (!___ret___)                                                        \
        {                                                                      \
            SIM_ERROR(__VA_ARGS__);                                            \
        }                                                                      \
    }

#define SIM_ASSERT(x) SIM_ASSERT_INFO((x), #x)