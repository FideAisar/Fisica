#include <iostream>

int main() {
  int pippo = 12;
  int &pasticcio = pippo;

  std::cout << pippo << " equal to " << pasticcio << std::endl;
  pippo = 15;
  std::cout << pippo << " equal to " << pasticcio << std::endl;

  return 0;
}
