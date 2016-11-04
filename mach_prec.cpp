#include <iostream>
#include <cmath>
using namespace std;

int main(){
	double unit(1.0), subject(1.0), combo(2.0);
	bool loopcheck(true);
	while(combo>unit){
		subject *= 0.5;
		combo = unit+subject;
	}
	subject *= 2;
	cout << subject << endl;
	return 0;
}
	
		
		
