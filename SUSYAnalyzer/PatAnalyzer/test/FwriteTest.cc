#include <iostream>
#include <stdio.h>
#include <time.h>
#include <ctime>
int main(){

	FILE *output;
	output = fopen("FwriteTest.txt","wb");
	
	std::time_t timer;// = std::time();
	
	time ( &timer);
	
	std::cout<<std::asctime(std::localtime(&timer))<<"\n";

	for(int x=0;x<1000;x++){

		unsigned short word[35651584] = {0};

		for(int i=0;i<35651584;i++){


			word[i] = x+i;


		}
		
		fwrite(word,sizeof(word),sizeof(word),output);

	}
	
	time ( &timer );
	std::cout<<std::asctime(std::localtime(&timer))<<"\n";
	
	fclose(output);

	return 0;
}
