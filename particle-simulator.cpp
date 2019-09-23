#include <iostream>
#include <cmath>
#include <vector> 
using namespace std;

int n, l, r, s;

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
				//same as above but for y collision before x
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
				}
						y += (1-yCollide) * vY;
					}
			}
        }
}; 

double timeCollision(Particle, Particle);

int main ()
{
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
	
	/*try to check collision
	in case to check function works LOL
	cout << "Collision check between particle 0 and 1" << endl;
	cout << timeCollision(particles[0], particles[1]) << endl;
    */
	
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

//Input: 2 Particles
//Output: Returns time taken before collision occurs if they collide, negative value otherwise.
double timeCollision(Particle first, Particle second) {
	
	double c = pow((first.x-second.x), 2) + pow((first.y - second.y), 2) - 4*r*r;
	double b = 2*((first.x - second.x)*(first.vX - second.vX) + (first.y - second.y)*(first.vY-second.vY));
	double a = pow((first.vX-second.vX), 2) + pow((first.vY - second.vY), 2);
	
	//check for solution
	if (b*b-4*a*c < 0) {
		return -100.0;
	}
	
	//else if there is a solution
	double solfirst = (-sqrt(b*b-4*a*c)-b)/(2*a);
	double solsecond = (-b+sqrt(b*b-4*a*c))/(2*a);
	return solfirst;
}