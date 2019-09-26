import random
import math

num_particles = int(input().rstrip())
box_size = int(input().rstrip())
radius = int(input().rstrip())
timesteps = int(input().rstrip())

print(num_particles)
print(box_size)
print(radius)
print(timesteps)
print("perf")

random.seed(114514)
pos = []
for i in range(0, num_particles):
    while True:
        acceptable = True
        #prevent particle spawning inside a wall...
        x = random.random() * (box_size - 2*radius) + radius
        y = random.random() * (box_size - 2*radius) + radius
        #retry to make particle if they spawn inside another particle...
        for j in range(0, len(pos)):
            dist = math.pow(x-pos[j][0], 2) + math.pow(x-pos[j][1], 2)
            if dist < 4 * radius * radius:
                acceptable = False
        if acceptable:
            break
    vX = random.random() * (box_size/(4*radius) - box_size/(32*radius)) + box_size/(32*radius)
    vY = random.random() * (box_size/(4*radius) - box_size/(32*radius)) + box_size/(32*radius)
    print(str(i) + " " + str(x) + " " + str(y) + " " + str(vX) + " " + str(vY))
    pos.append([x, y])
    