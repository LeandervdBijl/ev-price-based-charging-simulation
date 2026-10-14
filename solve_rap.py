import statistics
import numpy as np
from scipy.optimize import root_scalar
from decimal import Decimal, getcontext
import copy


#This file will be used to create classes for RAPs and solve it.

class RAP:
    def __init__(self, C, l_t, u_t, profile_others, d):
        #We assume g_t = x_t and f_t = a_t x_t^2 + b_t x_t with bounds l_t and u_t
        self.C = C
        self.n = len(profile_others)
        self.l_t = l_t
        self.u_t = u_t
        self.profile_others = profile_others
        self.d = d


    def evaluate_solution(self, x):
        objective = 0
        summ = 0
        for t in range(len(x)):
            objective += x[t]*(self.profile_others[t] + x[t])**self.d
            summ += x[t]

        return(objective, summ)

    

    def solve_rap_MB2(self):
        N = [i for i in range(0,self.n)]
        b_k = self.C

        #Initialize step
        x_opt = [0 for j in N]

        #Use water filling if d is infty
        if self.d == 'infty':
            #Step 1: if first iteration: find all minimum values of profile other, as well as the second lowest. Else, only find second lowest.
            #Make deep copy
            new_profile = copy.deepcopy(self.profile_others)

            stop = False

            while not stop:

                min_value = min(new_profile)
            
                # Find all indices where this minimum value occurs
                min_indices = [i for i, val in enumerate(new_profile) if val == min_value]

                # Find the second lowest value
                unique_values = sorted(set(new_profile))  # Sort unique values

                if len(unique_values) == 1:
                    stop = True

                    for index in range(len(x_opt)):
                        x_opt[index] += b_k/len(x_opt)

                else:
                    second_lowest = unique_values[1] if len(unique_values) > 1 else None  # Handle case where no second lowest exists

                    if len(min_indices)*(second_lowest - min_value) <= b_k:
                        for index in min_indices:
                            x_opt[index] += second_lowest - min_value
                            new_profile[index] += second_lowest - min_value
                        b_k -= len(min_indices)*(second_lowest - min_value)

                    else:
                        for index in min_indices:
                            x_opt[index] += b_k/len(min_indices)
                        stop = True

            #final step
            return x_opt


        #Calculate this according to (5) of Patricksson
        mu_l = [0 for j in N]
        mu_u = [0 for j in N]

        #avoid problems with zero
        for indx in range(len(self.profile_others)):
            self.profile_others[indx] += 1

        for j in range(self.n):
                mu_l[j] = (self.d)*self.l_t[j]*(self.profile_others[j] + self.l_t[j])**(self.d-1) + (self.l_t[j]+self.profile_others[j])**(self.d)
                mu_u[j] = (self.d)*self.u_t[j]*(self.profile_others[j] + self.u_t[j])**(self.d-1) + (self.u_t[j]+self.profile_others[j])**(self.d)

        mu_k = mu_l + mu_u
        mu_k = list(set(mu_k))

        iterations = 0
        while len(mu_k) > 0:
            iterations += 1

            mu_m = statistics.median(mu_k)

            #Calculate betas
            beta_l = 0
            beta_u = 0
            for j in N:
                if mu_m <= mu_l[j]:
                    beta_l += self.l_t[j]
                if mu_m >= mu_u[j]:
                    beta_u += self.u_t[j]


            #Calculate delta
            delta = 0
            for j in N:
                #Calculate x_mu_m from (5)
                if mu_u[j] >= mu_m:
                    if mu_m >= mu_l[j]:
                        x_mu_m = self.solve_x(mu_m, j)
                        delta += x_mu_m
            delta += beta_u + beta_l

            #Stopping condition
            if abs(b_k - delta)/(self.C) < 0.00001:
                mu_k = []

            if delta < b_k:
                #Update and fix lower bounds
                mu_k = list(filter(lambda x: x >= mu_m, mu_k))
                indexlist = []
                for j in N:
                    if mu_m >= mu_u[j]:
                        indexlist = indexlist + [j]
                        x_opt[j] = self.u_t[j]
                for j in indexlist:
                    N = list(filter(lambda x: x != j, N))
                b_k = b_k - beta_u
            else:
                if delta > b_k:
                    #Update and fix upper bounds
                    mu_k = list(filter(lambda x: x <= mu_m, mu_k))
                    indexlist = []
                    for j in N:
                        if mu_m <= mu_l[j]:
                            indexlist = indexlist + [j]
                            x_opt[j] = self.l_t[j]
                    for j in indexlist:
                        N = list(filter(lambda x: x != j, N))
                    b_k = b_k - beta_l
                else:
                    mu_k = []
            if len(mu_k) > 0:
                if mu_m not in mu_k:
                    mu_k = mu_k + [mu_m]

        for j in N:
            x_opt[j] = self.solve_x(mu_m, j)
            mu_opt = mu_m

        #remove the 1s added at the start
        for indx in range(len(self.profile_others)):
            self.profile_others[indx] -= 1

        return x_opt

    

    def solve_x(self, mu, j):
            def f(x):
                return self.d * x * (x + self.profile_others[j])**(self.d - 1) + (x + self.profile_others[j])**self.d - mu

            # return sol.root
            def binary_search(f, d, a, b, tol=1e-6):
                while b - a > tol:
                    mid = (a + b) / 2
                    if f(mid) < d:
                        a = mid  # Move right
                    else:
                        b = mid  # Move left
                return (a + b) / 2  # Approximate solution

            solution = binary_search(f, 0, self.l_t[j],self.u_t[j])
            return solution