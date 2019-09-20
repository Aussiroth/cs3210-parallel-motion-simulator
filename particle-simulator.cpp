#include <iostream>
#include <cmath>
#include <vector> 
using namespace std;

class Particle 
{ 
    // Access specifier 
    public: 
        operator string() const { 
            char buffer [100];;
            sprintf(buffer, "Particle %d: %.8lf %.8lf %.8lf %.8lf", i, x, y, vX, vY); 
            return buffer;
        }

        Particle(int i, double x, double y, double vX, double vY, int l) 
        {
            this -> i = i;
            this -> x = x;
            this -> y = y;
            this -> vX = vX;
            this -> vY = vY;
	    this -> l = l;
        }
  
        // Data Members
        int i;
	int l;
        double x;
        double y;
        double vX;
        double vY; 

        void move() {
	    //check for x wall collisions
	    //check for y wall collisions
	    double xCollide = vX < 0 ? x/(0-vX) : ((double)l-x)/vX;
	    double yCollide = vY < 0 ? y/(0-vY) : ((double)l-y)/vY;
	    if (xCollide >= 1 && yCollide >= 1) {
            	x += vX;
            	y += vY;
	    }
	    else {
	    	if (xCollide < yCollide) {
		    x += xCollide * vX;
		    y += xCollide * vY;
		    vX = -vX;
		    //after handling x collision, need to stop the ball at the edge of box if it collides with y too
		    if (yCollide < 1) {
		    	x += (yCollide-xCollide) * vX;
			y += (yCollide-xCollide) * vY;
		    }
		    else {
		    	x += (1-xCollide) * vX;
			y += (1-xCollide) * vY;
		    }
		}
		else {
		    x += yCollide * vX;
		    y += yCollide * vY;
		    vY = -vY;
		    if (xCollide < 1) {
		    	x += (xCollide-yCollide) * vX;
			y += (xCollide-yCollide) * vY;
		    }
		    else {
		        x += (1-yCollide) * vX;
			y += (1-yCollide) * vY;
		    }
		}
	    }
        }
}; 

int main ()
{
    int n; // no. of particles
    int l; // size of square
    int r; // radius of particle
    int s; // no. of steps
    string command; // simulator command
    cin >> n >> l >> r >> s >> command;

    //
    vector<Particle> particles; 
    for (int i = 0; i < n; ++i)
    {
        int index; 
        double x;
        double y;
        double vX;
        double vY; 
        int count;
        count = scanf("%d %lf %lf %lf %lf", &index, &x, &y, &vX, &vY);
        if (count == EOF || count <= 0) break;

        particles.push_back(Particle (index, x, y, vX, vY, l));
    }
    
    if (!command.compare("print"))
    {
        cout << "Input read: " << endl;
        for (int i = 0; i < particles.size(); ++i)
        {
            
            cout << (string) particles[i] << endl;
        }

        for (int i = 0; i < s; ++i)
        {
            cout << "Timestep " << i << endl;
            for (int j = 0; j < particles.size(); ++j)
            {
                particles[j].move();
                cout << (string) particles[j] << endl;
            }
        }
    }
    return 0;
}
