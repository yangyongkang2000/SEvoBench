from evox import algorithms, problems, workflows, monitors, use_state
import numpy as np
import jax.numpy as jnp
from jax import random
import time
nruns=10
dim=30
best_fitness_list = np.zeros(nruns)
run_times = np.zeros(nruns)
for _ in range(nruns):
    shade= algorithms.SHADE(
    lb=jnp.full(shape=(dim,), fill_value=-100),
    ub=jnp.full(shape=(dim,), fill_value=100),
    pop_size=100,
    )
    problem = problems.numerical.Rosenbrock()
    monitor = monitors.EvalMonitor(full_fit_history=False)
    workflow=workflows.StdWorkflow(shade,problem,monitors=[monitor])
    key = random.PRNGKey(np.random.randint(1e5))
    start_time = time.time()
    state = workflow.init(key)
    for i in range(299):
        state = workflow.step(state)
    run_times[_] = time.time() - start_time
    best_fitness, state = use_state(monitor.get_best_fitness)(state)
    best_fitness_list[_]=best_fitness
mean_fitness = np.mean(best_fitness_list)
var_fitness = np.var(best_fitness_list)
mean_time = np.mean(run_times)
print(f"\nResults after {nruns} runs:")
print(f"Mean best fitness: {mean_fitness:.3e}")
print(f"Fitness variance: {var_fitness:.3e}") 
print(f"Average runtime: {mean_time:.2f} seconds")
