import numpy as np

# Objective function to be optimized (you need to define your own)
def objective_function(properties):
    # Replace this with your actual objective function
    return np.sum(np.square(properties))

class Particle:
    def __init__(self, dimension):
        self.position = np.random.rand(dimension)  # Initial random position
        self.velocity = np.random.rand(dimension)  # Initial random velocity
        self.best_position = self.position.copy()
        self.best_fitness = float('inf')

def particle_swarm_optimization(objective_function, num_particles, num_dimensions, max_iterations, inertia_weight, c1, c2):
    particles = [Particle(num_dimensions) for _ in range(num_particles)]
    global_best_position = None
    global_best_fitness = float('inf')

    for _ in range(max_iterations):
        for particle in particles:
            fitness = objective_function(particle.position)

            # Update personal best
            if fitness < particle.best_fitness:
                particle.best_fitness = fitness
                particle.best_position = particle.position.copy()

            # Update global best
            if fitness < global_best_fitness:
                global_best_fitness = fitness
                global_best_position = particle.position.copy()

        for particle in particles:
            # Update velocity and position
            inertia_term = inertia_weight * particle.velocity
            cognitive_term = c1 * np.random.rand(num_dimensions) * (particle.best_position - particle.position)
            social_term = c2 * np.random.rand(num_dimensions) * (global_best_position - particle.position)

            particle.velocity = inertia_term + cognitive_term + social_term
            particle.position += particle.velocity

    return global_best_position, global_best_fitness

# Parameters
num_particles = 20
num_dimensions = 8
max_iterations = 100
inertia_weight = 0.5
c1 = 1.5  # Cognitive coefficient
c2 = 1.5  # Social coefficient

# Run PSO
best_solution, best_fitness = particle_swarm_optimization(objective_function, num_particles, num_dimensions, max_iterations, inertia_weight, c1, c2)

print("Best solution:", best_solution)
print("Best fitness:", best_fitness)

