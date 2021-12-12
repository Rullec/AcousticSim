#pragma once
#include <vector>
#include <string>

class cStringUtil
{
public:
    static std::vector<std::string> SplitString(const std::string &raw_string, const std::string &delimiter);
    static void RemoveEmptyLine(std::vector<std::string> &lines);
};