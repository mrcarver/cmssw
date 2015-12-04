#include <fstream>
#include <iostream>
#include <stdio.h>
#include <ctime>
#include <time.h>
int main(){

	std::ofstream file2("LUT2.dat",std::ofstream::out);
	
	std::time_t timer;// = std::time();
	
	time ( &timer);
	
	std::cout<<std::asctime(std::localtime(&timer))<<"\n";
	
	//for(int i=0;i<100;i++){//100
	//for(int i=0;i<1000;i++){//1k
	//for(int i=0;i<1000000;i++){//1M
	for(int i=0;i<1000000000;i++){//1G
	
		file2 << i%9 << std::endl;
	
	
	}
	
	time ( &timer);
	
	std::cout<<std::asctime(std::localtime(&timer))<<"\n";
	
	return 0;

}
