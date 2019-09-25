#include <iostream>
#include <cmath>
#include <vector> 
#include <algorithm>
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

        Particle();

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
			cout << x << " " << y << " " << vX << " " << vY << " " << xCollide << " " << yCollide << endl;
			if (xCollide >= 1 && yCollide >= 1) {
				x += vX;
				y += vY;
			}
			else 
            {
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

        int getIndex()
        {
            return this->i;
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

class CollisionEvent 
{
    public:
        Particle first;
        double time;

        CollisionEvent(Particle first, double time)
        {
            this->first = first;
            this->time = time;
        }

        virtual void execute();

        double getTime()
        {
            return this->time;
        }

        double getSmallestIndex()
        {
            return first.getIndex();
        }
};

class ParticleCollisionEvent: public CollisionEvent
{
    public:
        Particle second;

        ParticleCollisionEvent(Particle first, Particle second, double time)
        : CollisionEvent(first, time)
        {
            this->second = second;

        }
    
        void execute() 
        {
            //move them to proper position first
			first.x += time * first.vX;
			first.y += time * first.vY;
			second.x += time * second.vX;
			second.y += time * second.vY;
			
			//perform collision here
            //find normal vector
			double normalX = first.x - second.x;
			double normalY = first.y - second.y;
			double normalMag = sqrt(pow(normalX, 2) + pow(normalY, 2));
			normalX = normalX/normalMag; normalY = normalY/normalMag;
			double tangentX = -normalY;
			double tangentY = normalX;
			
			//compute velocity vectors wrt to normal and tangent
			double vFirstNormal = normalX * first.vX + normalY * first.vY;
			double vFirstTangent = tangentX * first.vX + tangentY * first.vY;
			double vSecondNormal = normalX * second.vX + normalY * second.vY;
			double vSecondTangent = tangentX * second.vX + tangentY * second.vY;
			
			//collision simply swaps velocities
			vFirstNormal = vSecondNormal;
			vSecondNormal = vFirstNormal;
			
			first.vX = vFirstNormal * normalX + vFirstTangent * tangentX;
			first.vY = vFirstNormal * normalY + vFirstTangent * tangentY;
			second.vX = vSecondNormal * normalX + vSecondTangent * tangentX;
			second.vY = vSecondNormal * normalY + vSecondTangent * tangentY;
			
			//Continue to move them here
			//Check for wall collisions and stop the particle at wall if so
			double timeToMove;
			double xCollide = first.vX < 0 ? (first.x-r)/(0-first.vX) : ((double)l-r-first.x)/first.vX;
			double yCollide = first.vY < 0 ? (first.y-r)/(0-first.vY) : ((double)l-r-first.y)/first.vY;
			if (xCollide >= 1-time && yCollide >= 1-time) 
			{
				timeToMove = 1-time;
			}
			else
			{
				timeToMove = min(xCollide, yCollide);
			}
			first.x += timeToMove * first.vX;
			first.y += timeToMove * first.vY;
			
			xCollide = first.vX < 0 ? (first.x-r)/(0-first.vX) : ((double)l-r-first.x)/first.vX;
			yCollide = first.vY < 0 ? (first.y-r)/(0-first.vY) : ((double)l-r-first.y)/first.vY;
			if (xCollide > 1-time && yCollide > 1-time) 
			{
				timeToMove = 1-time;
			}
			else 
			{
				timeToMove = min(xCollide, yCollide);
			}
			second.x += timeToMove * second.vX;
			second.y += timeToMove * second.vY;
        }

        double getSmallestIndex()
        {
            return first.getIndex() < second.getIndex() ? first.getIndex() : second.getIndex();
        }
};

class WallCollisionEvent: public CollisionEvent
{
    public:

        WallCollisionEvent(Particle first, double time)
        : CollisionEvent(first, time){}

        void execute() {
            //check for x wall collisions
			//check for y wall collisions
			double xCollide = first.vX < 0 ? first.x/(0-first.vX) : ((double)first.l-first.x)/first.vX;
			double yCollide = first.vY < 0 ? first.y/(0-first.vY) : ((double)first.l-first.y)/first.vY;
			if (xCollide < yCollide) {
				first.x += xCollide * first.vX;
				first.y += xCollide * first.vY;
				first.vX = -first.vX;
				//after handling x collision, need to stop the ball at the edge of box if it collides with y too
				if (yCollide < 1) {
					first.x += (yCollide-xCollide) * first.vX;
					first.y += (yCollide-xCollide) * first.vY;
				}
				else {
					first.x += (1-xCollide) * first.vX;
					first.y += (1-xCollide) * first.vY;
				}
			}
			//same as above but for y collision before x
			else {
				first.x += yCollide * first.vX;
				first.y += yCollide * first.vY;
				first.vY = -first.vY;
				if (xCollide < 1) {
					first.x += (xCollide-yCollide) * first.vX;
					first.y += (xCollide-yCollide) * first.vY;
				}
				else {
					first.x += (1-yCollide) * first.vX;
				}
				first.y += (1-yCollide) * first.vY;
			}
        }
};

// class NoCollisionEvent: public CollisionEvent
// {
//     public:
//         void execute() {
//             //simply move the particle
// 			first.x += first.vX;
// 			first.y += first.vY;
//         }
// };

struct collisionEventGreaterThan {
    bool operator()(CollisionEvent a, CollisionEvent b) const {
        if (a.getTime() == b.getTime()) return a.getSmallestIndex() > b.getSmallestIndex(); 
        return a.getTime() > b.getTime();
    }
};

void moveParticlesParallel(vector<Particle> particles);
double timeParticleCollision(Particle, Particle);
double timeWallCollision(Particle);
void performCollision(Particle&, Particle&);

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
	in case to check function works*/
	/*
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
		performCollision(particles[0], particles[1]);
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
    // time of particle-particle collisions
    // JaggedMatrix particleCollisionTimes = JaggedMatrix(n);
    // time of particle-wall collisions
    // double wallCollisionTimes[n] = {};

    vector<CollisionEvent> events;
    
    // calculate collision times
    # pragma omp parallel for
    for (int i = 0; i < n-1; ++i)
    {
        double wallCollisionTime = timeWallCollision(particles[i]);
        if (wallCollisionTime < 1)
        {
            WallCollisionEvent event = WallCollisionEvent(particles[i], wallCollisionTime);
            events.push_back(event);
        }

        # pragma omp parallel for
        for (int j = i+1; j < n-1; ++j)
        {
            double particleCollisionTime = timeParticleCollision(particles[i], particles[j]);
            if (particleCollisionTime >= 0 && particleCollisionTime < 1)
            {
                ParticleCollisionEvent event = ParticleCollisionEvent(particles[i],particles[j], particleCollisionTime);
                events.push_back(event);
            }
        }
    }

    // TODO: create min heap for collisions
    make_heap(events.begin(), events.end(), collisionEventGreaterThan());

    // TODO: process min heap of collisions
    bool map[n] = { false };
    #pragma omp parallel

    // TODO: free all memory
    // particleCollisionTimes.destroy();
}

//Input: 2 Particles
//Output: Returns time taken before collision occurs if they collide, negative value otherwise.
double timeParticleCollision(Particle first, Particle second)
{
	
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
double timeWallCollision(Particle particle)
{
	//check for x wall, y wall collisions
    double xCollide = particle.vX < 0 ? particle.x/(0-particle.vX) : ((double)l-particle.x)/particle.vX;
    double yCollide = particle.vY < 0 ? particle.y/(0-particle.vY) : ((double)l-particle.y)/particle.vY;
	return min(xCollide, yCollide);
}

//Takes in 2 particles by reference, assumes they have been moved to a position of collision
//Changes the velocities of the 2 particles as if it was a perfectly elastic collision of particles of equal mass.
void performCollision(Particle& first, Particle& second)
{
	//find normal vector
	double normalX = first.x - second.x;
	double normalY = first.y - second.y;
	double normalMag = sqrt(pow(normalX, 2) + pow(normalY, 2));
	normalX = normalX/normalMag; normalY = normalY/normalMag;
	double tangentX = -normalY;
	double tangentY = normalX;
	
	//compute velocity vectors wrt to normal and tangent
	double vFirstNormal = normalX * first.vX + normalY * first.vY;
	double vFirstTangent = tangentX * first.vX + tangentY * first.vY;
	double vSecondNormal = normalX * second.vX + normalY * second.vY;
	double vSecondTangent = tangentX * second.vX + tangentY * second.vY;
	//collision simply changes normal velocities to be in opposite direction
	
	vFirstNormal = -vFirstNormal;
	vSecondNormal = -vSecondNormal;
	
	first.vX = vFirstNormal * normalX + vFirstTangent * tangentX;
	first.vY = vFirstNormal * normalY + vFirstTangent * tangentY;
	second.vX = vSecondNormal * normalX + vSecondTangent * tangentX;
	second.vY = vSecondNormal * normalY + vSecondTangent * tangentY;
}