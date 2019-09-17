#include <iostream>
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

        Particle(int i, double x, double y, double vX, double vY) 
        {
            this -> i = i;
            this -> x = x;
            this -> y = y;
            this -> vX = vX;
            this -> vY = vY;
        }
  
        // Data Members
        int i; 
        double x;
        double y;
        double vX;
        double vY; 

        void move() {
            x += vX;
            y += vY;
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

        particles.push_back(Particle (index, x, y, vX, vY));
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