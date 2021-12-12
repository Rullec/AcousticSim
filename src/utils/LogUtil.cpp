#include "LogUtil.h"
#include <cstdarg>

void cLogUtil::SetLoggingLevel(const std::string &level_name)
{
    spdlog::set_level(static_cast<spdlog::level::level_enum>(
        GetLevelEnumFromString(level_name)));
}

int cLogUtil::GetLevelEnumFromString(const std::string &level_name)
{
    if (level_name == "trace")
    {
        return spdlog::level::trace;
    }
    else if (level_name == "debug")
    {
        return spdlog::level::debug;
    }
    else if (level_name == "info")
    {
        return spdlog::level::info;
    }
    else if (level_name == "warn")
    {
        return spdlog::level::warn;
    }
    else if (level_name == "error")
    {
        return spdlog::level::err;
    }
    else if (level_name == "critical")
    {
        return spdlog::level::critical;
    }
    else if (level_name == "off")
    {
        return spdlog::level::off;
    }
    else
    {
        SIM_ERROR("Unknown logging level [{}]. Levels = trace, debug, info, "
                  "warn, error, "
                  "critical, off",
                  level_name);
    }
}

tLogger cLogUtil::CreateLogger(const std::string &loggername)
{
    tLogger logger = spdlog::stdout_color_mt(loggername);
    logger->set_level(spdlog::level::trace);
    logger->set_pattern("[%^%l%$]%v");
    // if you want to add time
    // logger->set_pattern("[%n]%v");
    return logger;
}

void cLogUtil::DropLogger(const std::string &loggername)
{
    spdlog::drop(loggername);
}

tLogger cLogUtil::GetLogger(const std::string &loggername)
{
    return spdlog::get(loggername);
}

void InfoPrintf(const tLogger &logger, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    // size_t real_size = vsnprintf(buf, buf_size, fmt, ap);

    cLogUtil::Printf(logger, eLogLevel::info, fmt, ap);
    va_end(ap);
}

void WarnPrintf(const tLogger &logger, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    // size_t real_size = vsnprintf(buf, buf_size, fmt, ap);

    cLogUtil::Printf(logger, eLogLevel::warn, fmt, ap);
    va_end(ap);
}

void ErrorPrintf(const tLogger &logger, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    // size_t real_size = vsnprintf(buf, buf_size, fmt, ap);

    cLogUtil::Printf(logger, eLogLevel::err, fmt, ap);
    va_end(ap);
}

void DebugPrintf(const tLogger &logger, const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    // size_t real_size = vsnprintf(buf, buf_size, fmt, ap);

    cLogUtil::Printf(logger, eLogLevel::debug, fmt, ap);
    va_end(ap);
}

void cLogUtil::Printf(const tLogger &logger, eLogLevel level, const char *fmt,
                      va_list args)
{
    size_t real_size = vsnprintf(buf, buf_size, fmt, args);

    const std::string &log = std::string(buf, std::min(real_size, buf_size));
    switch (level)
    {
    case eLogLevel::info:
        logger->info(log);
        break;
    case eLogLevel::debug:
        logger->debug(log);
        break;
    default:
        logger->error(log);
        break;
    }
}

// // offer a template specialization for Eigen variables
// #include "MathUtil.h"
// template <> struct fmt::formatter<Eigen::Vector4d>
// {
//     // Presentation format: 'f' - fixed, 'e' - exponential.
//     char presentation = 'f';

//     // Parses format specifications of the form ['f' | 'e'].
//     constexpr auto parse(format_parse_context &ctx)
//     {
//         // auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) //
//         // c++11 [ctx.begin(), ctx.end()) is a character range that contains
//         a
//         // part of the format string starting from the format specifications
//         to
//         // be parsed, e.g. in
//         //
//         //   fmt::format("{:f} - point of interest", point{1, 2});
//         //
//         // the range will contain "f} - point of interest". The formatter
//         should
//         // parse specifiers until '}' or the end of the range. In this
//         example
//         // the formatter should parse the 'f' specifier and return an
//         iterator
//         // pointing to '}'.

//         // Parse the presentation format and store it in the formatter:
//         auto it = ctx.begin(), end = ctx.end();
//         if (it != end && (*it == 'f' || *it == 'e'))
//             presentation = *it++;

//         // Check if reached the end of the range:
//         if (it != end && *it != '}')
//             throw format_error("invalid format");

//         // Return an iterator past the end of the parsed range:
//         return it;
//     }

//     // Formats the point p using the parsed format specification
//     (presentation)
//     // stored in this formatter.
//     template <typename FormatContext>
//     auto format(const Eigen::Vector4d &p, FormatContext &ctx)
//     {
//         // auto format(const point &p, FormatContext &ctx) ->
//         // decltype(ctx.out()) // c++11 ctx.out() is an output iterator to
//         write
//         // to.
//         return format_to(ctx.out(),
//                          presentation == 'f' ? "({:.4f}, {:.4f}, {:.4f},
//                          {:.4f})"
//                                              : "({:.4e}, {:.4e}, {:.4e},
//                                              {:.4e})",
//                          p[0], p[1], p[2], p[3]);
//     }
// };