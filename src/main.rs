// Combi
// 20/02/24
// Neil Crago <n.j.crago@gmail.com>
//
// A program to explore Combinatorics using a Genetic Algorithm
// to solve the TSP (Travelling Salesman Problem)
// Every time i run this thing I get different answers,
// so must be a problem, though probabalistic it's too weird
// possibly unusable, bruteforce may be better
//

use rand::prelude::SliceRandom;
use rand::{thread_rng, Rng}; 
use std::f64;

const TOURNAMENT_SIZE: i32 = 2;
const RESTART_FREQUENCY: usize = 20;

// This struct represents a chromosome (city order) in the population
#[derive(Debug, Clone)]
struct Chromosome {
    genes: Vec<usize>,
    fitness: f64,
}

// This struct represents a city with its coordinates
#[derive(Clone, Debug)]
struct City {
    x: f64,
    y: f64,
}

// Function to calculate the distance between two cities
fn distance(city1: &City, city2: &City) -> f64 {
    let (x1, y1) = (city1.x, city1.y);
    let (x2, y2) = (city2.x, city2.y);
    f64::sqrt((x2 - x1).powi(2) + (y2 - y1).powi(2))
}

// This function calculates the total distance of a given tour
fn tour_length(cities: &[City], order: &[usize]) -> f64 {
    let mut total_distance = 0.0;
    for i in 0..cities.len() - 1 {
        let current_city = &cities[order[i]];
        let next_city = &cities[order[i + 1]];
        total_distance += distance(current_city, next_city);
    }
    // Add distance from last city back to first city to complete the loop
    total_distance += distance(&cities[order[cities.len() - 1]], &cities[order[0]]);
    total_distance
}

fn generate_initial_population(cities: &[City], population_size: usize) -> Vec<Chromosome> {
    let mut rng = thread_rng();
    let mut population = Vec::new();

    for _ in 0..population_size {
        let mut chromosome = Chromosome {
            genes: Vec::new(),
            fitness: 0.0,
        };

        // Create a vector of unique indices (0..cities.len())
        let mut indices: Vec<usize> = (0..cities.len()).collect();
        indices.shuffle(&mut rng);

        // Assign shuffled indices to the chromosome's genes
        chromosome.genes = indices;
        population.push(chromosome);
    }

    // Print the initial population for debugging
    // println!("Initial population: {:?}", population);

    population
}

fn calculate_fitness(population: &mut [Chromosome], cities: &[City]) {
    for chromosome in population.iter_mut() {
        chromosome.fitness = 1.0 / tour_length(cities, &chromosome.genes);
    }
}

fn selection(population: &[Chromosome], selection_size: usize) -> Vec<Chromosome> {
    let mut rng = thread_rng();
    let mut selected_chromosomes = Vec::new();
    for _ in 0..selection_size {
        let mut tournament = Vec::new();
        for _ in 0..TOURNAMENT_SIZE {
            // Adjust TOURNAMENT_SIZE (e.g., 2 or 3)
            let random_index = rng.gen_range(0..population.len());
            tournament.push(&population[random_index]);
        }
        let fittest_index = tournament
            .iter()
            .enumerate()
            .max_by_key(|(_, c)| c.fitness as i32)
            .unwrap()
            .0;
        selected_chromosomes.push(tournament[fittest_index].clone());
    }
    selected_chromosomes
}

fn crossover(parent1: &Chromosome, parent2: &Chromosome) -> (Chromosome, Chromosome) {
    let mut rng = thread_rng();
    let crossover_point = rng.gen_range(1..parent1.genes.len() - 1);
    let mut offspring1 = Chromosome {
        genes: Vec::new(),
        fitness: 0.0,
    };
    let mut offspring2 = Chromosome {
        genes: Vec::with_capacity(parent2.genes.len()),
        fitness: 0.0,
    };

    // Ordered crossover for offspring1
    offspring1
        .genes
        .extend(parent1.genes[0..crossover_point].iter().cloned());
    for gene in parent2.genes.iter() {
        if !offspring1.genes.contains(gene) {
            offspring1.genes.push(*gene);
        }
    }

    // Fill offspring2 with remaining genes from parent2 in order
    offspring2
        .genes
        .extend(offspring1.genes[0..crossover_point].iter().cloned());
    for gene in parent2.genes.iter() {
        if !offspring2.genes.contains(gene) {
            offspring2.genes.push(*gene);
        }
    }

    (offspring1, offspring2)
}

fn mutation(chromosome: &mut Chromosome, mutation_rate: f64) {
    let mut rng = thread_rng();
    for i in 0..chromosome.genes.len() {
        if rng.gen_bool(mutation_rate) {
            let random_index = rng.gen_range(0..chromosome.genes.len());
            chromosome.genes.swap(i, random_index);
        }
    }
}

fn genetic_algorithm_tsp(
    cities: &[City],
    population_size: usize,
    iterations: usize,
    mutation_rate: f64,
) -> Chromosome {
    let mut population = generate_initial_population(cities, population_size);
    for generation in 0..iterations {
        if generation % RESTART_FREQUENCY == 0 {
            // Restart every few generations
            population = generate_initial_population(cities, population_size);
        }
        calculate_fitness(&mut population, cities);
        let selected_parents = selection(&population, population_size / 2);
        let fittest_parent = population
            .iter()
            .max_by_key(|c| c.fitness as i32)
            .unwrap()
            .clone();
        let mut new_population = Vec::new();
        new_population.push(fittest_parent); // Elitism - include fittest parent
        for i in 0..population_size / 2 {
            let parent1 = &selected_parents[i];
            let parent2 = &selected_parents[(i + 1) % selected_parents.len()];
            let (mut child1, mut child2) = crossover(parent1, parent2);

            mutation(&mut child1, mutation_rate);
            mutation(&mut child2, mutation_rate);

            new_population.push(child1);
            new_population.push(child2);
        }
        population.append(&mut new_population);
        population.sort_by_key(|chromosome| chromosome.fitness as i32);
        population.truncate(population_size);
    }
    // Return the fittest chromosome
    population[0].clone()
}

fn main() {
    // Sample set of city coordinates
    let cities = vec![
        City { x: 1.0, y: 2.0 },
        City { x: 7.0, y: 11.0 },
        City { x: 3.0, y: 23.0 },
        City { x: 2.0, y: 80.0 },
        City { x: 9.0, y: 2.0 },
        City { x: 19.0, y: 11.0 },
        City { x: 29.0, y: 23.0 },
        City { x: 1.0, y: 80.0 },
    ];

    let population_size = 100;
    let iterations = 10;
    let mutation_rate = 0.001; // a reasonable range (e.g., 0.001 to 0.1)

    // let best_chromosome = genetic_algorithm_tsp(&cities, population_size, iterations, mutation_rate);
    let best_chromosome =
        genetic_algorithm_tsp(&cities, population_size, iterations, mutation_rate);

    let best_tour = &best_chromosome.genes;
    // println!("order = {:?}", &best_chromosome.genes);
    let best_distance = tour_length(&cities, best_tour);

    println!("Shortest tour distance (GA): {}", best_distance);
    println!("Tour order (GA): {:?}", best_tour);
}
