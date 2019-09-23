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

class JaggedMatrix
{
    public:
        int length;
        double **matrix;

        JaggedMatrix(int i) 
        {
            this->length = i;
            matrix = (double**) calloc(i, sizeof(double *));
            for (int k = 0; k < i; ++k) 
            {
                matrix[k] = (double *) calloc(k+1, sizeof(double));
            }
        }

        int get(int i, int j)
        {
            if (i < j) 
            {
                return matrix[j][i];
            }
            return matrix[i][j];
        }

        void put(int i, int j, double value)
        {
            if (i < j) 
            {
                matrix[j][i] = value;
            } else {
                matrix[i][j] = value;
            }
           
        }

        void destroy()
        {
            for (int k = 0; k < length; ++k)
            {
                free(matrix[k]);
            }
            free(matrix);
        }
};


double timeParticleCollision(Particle, Particle);
double timeWallCollision(Particle);
void moveParticlesParallel(vector<Particle> particles);

int main ()
{
    string command; // simulator command
    cin >> n >> l >> r >> s >> command;

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
	in case to check function works
	cout << "Collision check between particle 0 and 1" << endl;
	cout << timeParticleCollision(particles[0], particles[1]) << endl;
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
            moveParticlesParallel(particles);
            for (int j = 0; j < particles.size(); ++j)
            {
                // particles[j].move();
                cout << (string) particles[j] << endl;
            }
        }
    }
    return 0;
}


void moveParticlesParallel(vector<Particle> particles) 
{
    int n = particles.size();
    JaggedMatrix particleCollisionTimes = JaggedMatrix(n);
    double wallCollisionTimes[n] = {};
    // int minTimesIndex[n] = {};
    // for (int i = 0; i < n; ++i) minTimesIndex[i] = -1;
    // double taskCount[n] = {0};
    
    # pragma omp parallel for
    for (int i = 0; i < n-1; ++i)
    {
        wallCollisionTimes[i] = timeWallCollision(particles[i]);
        # pragma omp parallel for
        for (int j = i+1; j < n-1; ++j)
        {
            double t = timeParticleCollision(particles[i], particles[j]);
            particleCollisionTimes.put(i, j, t);
            // getCollisions();
        }
    }
    // free memory
    particleCollisionTimes.destroy();
}

//Input: 2 Particles
//Output: Returns time taken before collision occurs if they collide, negative value otherwise.
double timeParticleCollision(Particle first, Particle second) {
	
	//a, b and c are as in the quadratic formula representation.
	//t, the time taken for the 2 circles to touch, is the unknown variable we are solving for
	//by taking difference in circle centres, setting an unknown t for collision time, and then taking distance moved in time t,
	//we can solve for t such that the circle centers are <= 2r and therefore collide. 4r^2 is to solve for radius distance.
	double c = pow((first.x-second.x), 2) + pow((first.y - second.y), 2) - 4*r*r;
	double b = 2*((first.x - second.x)*(first.vX - second.vX) + (first.y - second.y)*(first.vY-second.vY));
	double a = pow((first.vX-second.vX), 2) + pow((first.vY - second.vY), 2);
	
	//check for solution
	if (b*b-4*a*c < 0) {
		return -100.0;
	}
	
	//else if there is a solution, the one with smaller value should be the main collision. Second value is after the 2 circles phase through each other
	double solfirst = (-sqrt(b*b-4*a*c)-b)/(2*a);
	double solsecond = (-b+sqrt(b*b-4*a*c))/(2*a);
	return solfirst;
}

//Input: 1 Particle
//Output: Returns time taken before collision occurs if it collides with wall, negative value otherwise.
double timeWallCollision(Particle particle) {
	//check for x wall collisions
    //check for y wall collisions
    double xCollide = particle.vX < 0 ? particle.x/(0-particle.vX) : ((double)l-particle.x)/particle.vX;
    double yCollide = particle.vY < 0 ? particle.y/(0-particle.vY) : ((double)l-particle.y)/particle.vY;
}