def generate_values(self,N): 
    rand_wait = np.random.randint(5) 
    rand_size = np.random.randint(1,20,N) 
    rand_dl = np.random.randint(1,10*N,N) 
    joblist = [[rand_size[c],rand_dl[c]] for c in range(N)] 
    jobs = {i:job for i,job in enumerate(joblist)} 
    rand_pen = np.random.rand(1) 
    return(rand_wait, jobs, rand_pen ) 
    
wait, jobs, pen = generate_values(1, 5) 

def evaluate(data , wait, jobs ,rand_pen): 
    N = len(data) 
    data = data 
    wait = wait 
    # Full evaluation where needed 
    values = np.zeros(len(data)) 
    for solution in range(len(data)): 
        time = 0 
        score = 0 
        for task in data[solution]: 
            size = jobs[task][0] 
            deadline = jobs[task][1] 
            penalty = (time+size) - deadline 
            if penalty > 0: 
                score += rand_pen*penalty 
                time += size + wait 
                values[solution] = score 
    return(values) 

data = [[0,1,2,3,4],[4,3,2,1,0]] 

evaluate(data, wait, jobs, pen)