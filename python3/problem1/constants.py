import os
vectorized_computation = True
results_folder = "./results"

max_iterations = 100000
number_of_children = 200

if not os.path.exists(results_folder):
    os.makedirs(results_folder, exist_ok=True)

