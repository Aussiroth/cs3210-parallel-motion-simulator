#include <iostream>
#include <cmath>
#include <vector> 
#include <algorithm>
#include <chrono>
#include <random>
#include <atomic>
using namespace std;

int n, l, r, s;

mt19937 rng;
random_device rd;

class Particle
{ 
	public: 
		operator string() const { 
			char buffer [200];
			snprintf(buffer, 200, "%d %.8lf %.8lf %.8lf %.8lf", i, x, y, vX, vY); 
			return buffer;
		}

		int i;
		double x;
		double y;
		double vX;
		double vY; 

		int pColl;
		int wColl;

		Particle() {};

		Particle(int i, double x, double y, double vX, double vY, int l) 
		{
			this -> i = i;
			this -> x = x;
			this -> y = y;
			this -> vX = vX;
			this -> vY = vY;
			this -> pColl = 0;
			this -> wColl = 0;
		}

		int getIndex()
		{
			return this->i;
		}
		
		string getFullRepresentation()
		{
			
			char buffer[200];
			snprintf(buffer, 200, "%d %.8lf %.8lf %.8lf %.8lf %d %d", i, x, y, vX, vY, pColl, wColl); 
			string res(buffer);
			return res;
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

		double get(int i, int j)
		{
			if (i < j) 
			{
				return matrix[j][i];
			}
			return matrix[i][j];
		}

		void set(int i, int j, double value)
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
		void print()
		{
			for (int i = 0; i < length; i++)
			{
				for (int j = 0; j <= i; j++)
					printf("%f ", matrix[i][j]);
				printf("\n");
			}
			
		}
};

class CollisionEvent 
{
	bool operator < (CollisionEvent other)
	{
		if (this->time == other.getTime()) return this->getSmallestIndex() > other.getSmallestIndex(); 
		return this->time > other.getTime();
	}

	public:
	Particle* first;
	double time;

	CollisionEvent() {}

	virtual ~CollisionEvent() {}

	CollisionEvent(Particle* first, double time)
	{
		this->first = first;
		this->time = time;
	}

	virtual void execute() {};

	double getTime()
	{
		return this->time;
	}

	double getSmallestIndex()
	{
		return (*first).getIndex();
	}
};

class ParticleCollisionEvent: public CollisionEvent
{
	public:
		bool operator == (ParticleCollisionEvent other)
		{
			int firstIndex = (*this->first).getIndex();
			int secondIndex = (*this->second).getIndex();
			int otherFirstIndex = (*other.first).getIndex();
			int otherSecondIndex = (*other.second).getIndex();
			return (firstIndex == otherSecondIndex && secondIndex == otherFirstIndex) ||
				(firstIndex == otherFirstIndex && secondIndex == otherSecondIndex);
		}
		Particle* second;


		ParticleCollisionEvent(Particle* first, Particle* second, double time)
			: CollisionEvent(first, time)
		{
			this->second = second;
		}

		void execute() 
		{
			if (first->getIndex() >= second->getIndex())
				return;
			
			//move them to proper position first
			first->x += this->time * first->vX;
			first->y += this->time * first->vY;
			second->x += this->time * second->vX;
			second->y += this->time * second->vY;

			//perform collision here
			//find normal vector
			double normalX = first->x - second->x;
			double normalY = first->y - second->y;
			double normalMag = sqrt(pow(normalX, 2) + pow(normalY, 2));
			normalX = normalX/normalMag; normalY = normalY/normalMag;
			double tangentX = -normalY;
			double tangentY = normalX;

			//compute velocity vectors wrt to normal and tangent
			double vFirstNormal = normalX * first->vX + normalY * first->vY;
			double vFirstTangent = tangentX * first->vX + tangentY * first->vY;
			double vSecondNormal = normalX * second->vX + normalY * second->vY;
			double vSecondTangent = tangentX * second->vX + tangentY * second->vY;

			//collision simply swaps velocities
			double temp = vFirstNormal;
			vFirstNormal = vSecondNormal;
			vSecondNormal = temp;

			first->vX = vFirstNormal * normalX + vFirstTangent * tangentX;
			first->vY = vFirstNormal * normalY + vFirstTangent * tangentY;
			second->vX = vSecondNormal * normalX + vSecondTangent * tangentX;
			second->vY = vSecondNormal * normalY + vSecondTangent * tangentY;

			//eliminate negative 0s
			if (first->vX == -0.0) first->vX = 0.0;
			if (first->vY == -0.0) first->vY = 0.0;
			if (second->vX == -0.0) second->vX = 0.0;
			if (second->vY == -0.0) second->vY = 0.0;

			//Continue to move them here
			//Check for wall collisions and stop the particle at wall if so
			double timeToMove;
			double xCollide = first->vX < 0 ? (first->x-r)/(0-first->vX) : ((double)l-r-first->x)/first->vX;
			double yCollide = first->vY < 0 ? (first->y-r)/(0-first->vY) : ((double)l-r-first->y)/first->vY;
			if (xCollide >= 1-this->time && yCollide >= 1-this->time) 
			{
				timeToMove = 1-this->time;
			}
			else
			{
				timeToMove = min(xCollide, yCollide);
			}
			first->x += timeToMove * first->vX;
			first->y += timeToMove * first->vY;

			xCollide = second->vX < 0 ? (second->x-r)/(0-second->vX) : ((double)l-r-second->x)/second->vX;
			yCollide = second->vY < 0 ? (second->y-r)/(0-second->vY) : ((double)l-r-second->y)/second->vY;
			if (xCollide >= 1-this->time && yCollide >= 1-this->time) 
			{
				timeToMove = 1-this->time;
			}
			else
			{
				timeToMove = min(xCollide, yCollide);
			}
			second->x += timeToMove * second->vX;
			second->y += timeToMove * second->vY;
			
			first->pColl++;
			second->pColl++;
		}

		double getSmallestIndex()
		{
			return (*first).getIndex() < (*second).getIndex() ? (*first).getIndex() : (*second).getIndex();
		}
};

class WallCollisionEvent: public CollisionEvent
{
	public:

		WallCollisionEvent(Particle* first, double time)
			: CollisionEvent(first, time){}

		void execute() {
			//check for x wall collisions
			//check for y wall collisions
			double xCollide = first->vX < 0 ? (first->x-r)/(0-first->vX) : ((double)l-first->x-r)/first->vX;
			double yCollide = first->vY < 0 ? (first->y-r)/(0-first->vY) : ((double)l-first->y-r)/first->vY;
			if (xCollide < yCollide) {
				first->x += xCollide * first->vX;
				first->y += xCollide * first->vY;
				first->vX = -first->vX;
				//after handling x collision, need to stop the ball at the edge of box if it collides with y too
				if (yCollide < 1) {
					first->x += (yCollide-xCollide) * first->vX;
					first->y += (yCollide-xCollide) * first->vY;
				}
				else {
					first->x += (1-xCollide) * first->vX;
					first->y += (1-xCollide) * first->vY;
				}
			}
			//collision with corner of box reverses both
			else if (xCollide == yCollide) {
				first->x += xCollide * first->vX;
				first->y += xCollide * first->vY;
				first->vX = -first->vX;
				first->vY = -first->vY;
				first->x += (1-xCollide) * first->vX;
				first->y += (1-xCollide) * first->vY;
			}
			//same as x collision but for y wall collision happening first
			else {
				first->x += yCollide * first->vX;
				first->y += yCollide * first->vY;
				first->vY = -first->vY;
				if (xCollide < 1) {
					first->x += (xCollide-yCollide) * first->vX;
					first->y += (xCollide-yCollide) * first->vY;
				}
				else {
					first->x += (1-yCollide) * first->vX;
					first->y += (1-yCollide) * first->vY;
				}
			}
			first->wColl++;
		}
};

class NoCollisionEvent: public CollisionEvent
{
	public:

		NoCollisionEvent(Particle* first)
			: CollisionEvent(first, 1.0)
		{}

		void execute() {
			//simply move the particle
			first->x += first->vX;
			first->y += first->vY;
		}
};

void moveParticles(vector<Particle*> particles);
double timeParticleCollision(Particle&, Particle&);
double timeWallCollision(Particle&);

int main ()
{
	string command; // simulator command
	cin >> n >> l >> r >> s >> command;

	rng.seed(rd());
	uniform_real_distribution<double> pos(r, l-r);
	uniform_real_distribution<double> velocity((double)l/(8*r), (double)l/4);
	vector<Particle*> particles; 
	int scanned;
	for (scanned = 0; scanned < n; ++scanned)
	{
		int index; 
		double x;
		double y;
		double vX;
		double vY; 
		int count;
		count = scanf("%d %lf %lf %lf %lf", &index, &x, &y, &vX, &vY);
		if (count == EOF || count <= 0) break;

		particles.push_back(new Particle(index, x, y, vX, vY, l));
	}
	for (int j = scanned; j < n; j++)
	{
		double x = pos(rng);
		double y = pos(rng);
		double vX = velocity(rng);
		double vY = velocity(rng);
		particles.push_back(new Particle(j, x, y, vX, vY, l));
	}

	auto start = chrono::high_resolution_clock::now();

	for (int i = 0; i < s; ++i)
	{	
		moveParticles(particles);
		if (!command.compare("print"))
		{
			for (int j = 0; j < particles.size(); ++j)
			{
				cout << i << " " << (string) *particles[j] << endl;
			}
		}
	}

	auto finish = std::chrono::high_resolution_clock::now();

	for (int j = 0; j < n; ++j)
	{
		cout << particles[j]->getFullRepresentation() << endl;
	}
	double timeTaken = (double)chrono::duration_cast<chrono::nanoseconds>(finish-start).count()/1000000000;
	 printf("Time taken: %.5f s\n", timeTaken);

	return 0;
}


void moveParticles(vector<Particle*> particles) 
{
	int n = particles.size();
	// time of particle-particle collisions
	JaggedMatrix particleCollisionTimes = JaggedMatrix(n);
	// time of particle-wall collisions
	double wallCollisionTimes[n] = {};

	vector<CollisionEvent> events;

	// calculate collision times
#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
		wallCollisionTimes[i] = timeWallCollision(*particles[i]);

#pragma omp parallel for
		for (int j = i+1; j < n; ++j)
		{
			double particleCollisionTime = timeParticleCollision(*particles[i], *particles[j]);
			particleCollisionTimes.set(i, j, particleCollisionTime);
		}
	}
	CollisionEvent* found[n] = { nullptr };
	int foundCount = 0;
	while (foundCount != n)
	{   
		CollisionEvent* temp[n];
#pragma omp parallel for
		for (int i = 0; i < n; ++i)
		{   
			// first assume no collision
			temp[i] = new NoCollisionEvent(particles[i]);

			// check for particle-wall collision
			if (wallCollisionTimes[i] < (*temp[i]).getTime() && wallCollisionTimes[i] < 1)
			{
				delete(temp[i]);
				temp[i] = new WallCollisionEvent(particles[i], wallCollisionTimes[i]);
			}

			// check for particle-particle collision
			for (int j = 0; j < n; ++j)
			{
				if (i == j) continue;

				double time = particleCollisionTimes.get(i, j);
				if (time > -1 && time < (*temp[i]).getTime() && time < 1 && found[j] == NULL) {
					delete(temp[i]);
					temp[i] = new ParticleCollisionEvent(particles[i], particles[j], time);
				}
			}
		}

		for (int i = 0; i < n; ++i)
		{
			if (found[i] != NULL) continue;

			CollisionEvent* e = temp[i];

			// particle-particle collision
			if(ParticleCollisionEvent* v = dynamic_cast<ParticleCollisionEvent*>(e))
			{

				int otherIndex = (*(*v).second).getIndex();
				if (ParticleCollisionEvent* v2 = dynamic_cast<ParticleCollisionEvent*>(temp[otherIndex]))
				{
					if (*v == *v2) 
					{
						found[i] = temp[i];
						++foundCount;
					}
				}

			}
			// particle-wall collision or no collision
			else
			{
				found[i] = temp[i];
				++foundCount;
			}
		}
	}

#pragma omp parallel for
	for (int i = 0; i < n; ++i)
	{
		(*found[i]).execute();
	}
	for (int i = 0; i < n; ++i) delete(found[i]);
	particleCollisionTimes.destroy();
}

// input: 2 Particles
// output: Returns time taken before collision occurs if they collide, negative value otherwise.
double timeParticleCollision(Particle& first, Particle& second)
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
		return 100000.0;
	}
	//else if there is a solution, the one with smaller value should be the main collision. Second value is after the 2 circles phase through each other
	double solfirst = (-sqrt(b*b-4*a*c)-b)/(2*a);
	if (solfirst > 0)
	{
		return solfirst;
	}
	solfirst = (sqrt(b*b-4*a*c)-b)/(2*a);
	if (solfirst > 0)
	{
		return 0;
	}
	return 100000.0;
}

// input: 1 Particle
// output: Returns time taken before collision occurs if it collides with wall, negative value otherwise.
double timeWallCollision(Particle& particle)
{
	//check for x wall, y wall collisions
	double xCollide = particle.vX < 0 ? (particle.x-r)/(0-particle.vX) : ((double)l-particle.x-r)/particle.vX;
	double yCollide = particle.vY < 0 ? (particle.y-r)/(0-particle.vY) : ((double)l-particle.y-r)/particle.vY;
	return min(xCollide, yCollide);
}
