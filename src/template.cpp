#include "template.h"
#include <iostream>

Template::Template(const std::string& str) : m_str(str) {}
void Template::print() {
    std::cout << m_str << std::endl;
}