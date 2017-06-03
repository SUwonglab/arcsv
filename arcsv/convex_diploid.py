import cvxpy as cvx
import numpy as np

def convex_diploid(lhs, ncs, pi_robust):
    n = len(lhs)
    w = cvx.Variable(n)
    print('building objective. . .')
    f = ((1-pi_robust) * np.asarray(lhs) + pi_robust) / np.asarray(ncs)
    print(f.shape)
    for i in range(f.shape[0]):
        print(np.sum(np.log(f[i,:])))
    obj = sum(cvx.log(f.T * w))
    # for i in range(len(lhs[0])):
    #     obj_sum += cvx.log(w.T * [(pi_robust + (1-pi_robust)*lh[i])/nc[i] for (lh, nc) in zip(lhs, ncs)])
    objective = cvx.Maximize(obj)

    constraints = [w >= 0,
                   sum(w) == 1]
    problem = cvx.Problem(objective, constraints)
    print('solving. . .')
    problem.solve()
    print(problem.status)
    return w.value, objective.value
