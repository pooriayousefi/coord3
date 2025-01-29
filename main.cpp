import coord3;

#include <stdexcept>
#include <exception>
#include <iostream>

using namespace R3;

// entry point
int main()
{
	try
	{
		cartesian_coordinate<float> c1{ 1.f,2.f,3.f }, c2{ 4.f,5.f,6.f };
		auto c3{ c1 + c2 };
		std::cout 
			<< "\nc1: " << c1
			<< "\nc2: " << c2
			<< "\nc3: " << c3
			<< std::endl;
		return 0;
	}
	catch (const std::exception& ex)
	{
		std::cerr << ex.what() << std::endl;
		return 1;
	}
}