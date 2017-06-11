#include <iostream>
using namespace std;
#include <vector>
 
   int main()
   {
    cout<<"life is lonly";
    
    struct Particle {
	    int id;
    };
    
    vector<Particle> x;
    
    
    
    x.push_back(Particle{0});
    x.push_back(Particle{1});
    x.push_back(Particle{2});
    
    for(Particle &p : x) {
        p = Particle{9};
    }
    
    for(Particle p : x) {
        cout << p << endl;
    }
    
    
    
    return 0;
    }