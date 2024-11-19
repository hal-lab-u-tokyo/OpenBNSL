#pragma once
#include <string>

class Template {
public:
    Template(const std::string& str);
private:
    std::string m_str;
public:
    void print();
};