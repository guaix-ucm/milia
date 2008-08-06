#include "milia/metric.h"
#include "milia/exception.h"

#include <stdexcept>
#include <iostream>

int main() {
std::cout << "Testing exceptions\n";
try {
milia::metric a(0, -10, 30);
}
catch(milia::exception& ex) {
std::cout << ex.what() << '\n';
return 0;
}
catch(...) {
std::cout << "Ooops\n";
return 1;
}
return 1;
}
