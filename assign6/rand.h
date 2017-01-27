#pragma once
// This is the random number generator used for the Computational Methods 
// Assignment 3, obtained from Numerical Recipes 3rd Edition. The parameters
// used are found to be the ones which provide the most "random" numbers by
// the people who wrote the code.


struct Ran {
        public:
        unsigned long long int u, v, w;
	// The constructor for this struct initialises three numbers that are
	// independently updated and then combined to produce a new psuedorandom
	// number. The seed "j" provided is used to initialise this chain.
        Ran(){
	}
	
	inline void seed(unsigned long long int j) {
                v = 4101842887655102017LL;
		w = 1;
		u = j^v; int64();
                v = u; int64();
                w = v; int64();
        }

	// The variable "v" is updated via a 64bit XOR shift method, whilst the
	// variable "w" is updated via a multiply with carry method. A new variable
	// "x" is produced essentially as a copy of "u", and another 64bit XOR shift
	// is carried out and then all three numbers are combined in a bitwise XOR,
	// to provide a random integer.
        inline unsigned long long int int64() {
                u= u*2862933555777941757LL + 7046029254386353087LL;
                u = u*2862933555777941757LL + 7046029254386353087LL;
                v ^= v >> 17;
                v ^= v << 31;
                v ^= v >> 8;
                w = 4294957665U*(w & 0xffffffff)+ (w>>32);
                unsigned long long int x = u^(u<<21);
                x ^= x >>35;
                x ^= x << 4;
                return (x+v)^w;
        }

	// The same method can be used to provide a double-precision floating point
	// variable between 0 and 1, by dividing by the maximum number obtainable from
	// the int64() method. This can also be done for unsigned integers by providing
	// a variable type conversion of the output of int64().
        inline double doub() { return 5.42101086242552217E-20*int64();}
        inline unsigned int int32() { return (unsigned int)int64(); }

};

