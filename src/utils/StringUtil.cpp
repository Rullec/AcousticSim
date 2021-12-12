#include "StringUtil.h"

std::vector<std::string> cStringUtil::SplitString(const std::string &s_, const std::string &delimiter)
{
    size_t pos = 0;
    std::string token;
    std::string s = s_;
    std::vector<std::string> split_string_array(0);
    while ((pos = s.find(delimiter)) != std::string::npos)
    {
        token = s.substr(0, pos);
        split_string_array.push_back(token);
        // std::cout << token << std::endl;
        s.erase(0, pos + delimiter.length());
    }
    split_string_array.push_back(s);
    return split_string_array;
}

void cStringUtil::RemoveEmptyLine(std::vector<std::string> &cont)
{
    std::vector<std::string>::iterator it = cont.begin();
    while (it != cont.end())
    {
        if ((*it).size() == 0)
        {
            it = cont.erase(it);
            it = cont.begin();
        }
        else
        {
            it++;
        }
    }
}