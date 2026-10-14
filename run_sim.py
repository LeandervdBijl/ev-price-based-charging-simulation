import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from solve_rap import *
import math
import copy


import datetime as dt
from FOCS import FOCS, FlowNet, FlowOperations, FOCSinstance
from Price import Price, Greedy
#The FOCS and Price files are obtained from: https://github.com/lwinschermann/FlowbasedOfflineChargingScheduler


#Get the data

#This is for bus data
#data = pd.read_csv('Data/output/data.csv')

d_values = [1]
#d_values = [0.5]

run_results = True
#d_values = [5]
plot_result = False
show_iterations = False
savedata = True


#To look at a single instance
investigate_single_instance = False

#Save a single profile
save_single_profile = False
save_total_profile = False
save_d_value = 'infty'

#Time tax, meaning that future usage is more expensive.
run_and_save_time_tax = False
time_taxes = [1.01, 1.03, 1.05, 1.07]
time_tax_d = 1.4

#date_range = ['2024-1-15','2024-01-16','2024-01-17','2024-01-18','2024-01-19']
year = '2024'


def linear_price(bus_profile, total_profile_except_bus):
    sum = 0
    for i in range(len(total_profile_except_bus)):
        sum += (total_profile_except_bus[i] + bus_profile[i])*bus_profile[i]
    return sum

#Every bus will now be its own actor, they have full flexibility to charge within the interval as long as its max 30kW and they have to be full at the end. Will react to price stimulations. We now make the class bus
class EV:
    def __init__(self, bus_id, start_time, end_time, tocharge, max_power):
        self.bus_id = bus_id
        self.start_time = start_time
        self.end_time = end_time
        self.intervals = self.end_time - self.start_time

        self.tocharge = tocharge
        self.max_power = max_power
        self.bus_profile = [0 for i in range(self.end_time-self.start_time)]
        self.cost = 0
    
    def react_to_prices(self, a_t, b_t, l_t, u_t, d, profile_of_others):
        #Adjust the cost to take the other profiles into account
        profile_others_range = profile_of_others[self.start_time: self.end_time]

        individual_problem = RAP(self.tocharge, l_t, u_t, profile_others_range, d)

        self.bus_profile = individual_problem.solve_rap_MB2()

        #Adjust cost (split in case where d can be infty (as simply the max peak I guess))
        if d == 'infty':
            #find all indices in which this player participates
            timelist = []
            for index in range(0,len(self.bus_profile)):
                if self.bus_profile[index] > 0:
                    timelist += [index]

            #Find the maximum peak in the overall profile within these indices
            total_profile = [profile_others_range[i] + self.bus_profile[i] for i in timelist]
    
            #The cost is equal to this peak size
            cost = max(total_profile)
        else:
            cost = 0
            for j in range(len(profile_others_range)):
                cost += self.bus_profile[j]*(self.bus_profile[j] + profile_others_range[j])**d
        self.cost = cost

class EV_Charging_Case:
    #Class is used to create an instance of the problem and solve it with the different methods
    def __init__(self, data):
        #Create a list of all buses
        self.data = data

        self.number_EVs = 0
        self.EV_list = []

        start_times = []
        end_times = []
        for index, row in data.iterrows():
            #Get the data
            self.number_EVs += 1
            start_time = int(row['t0_900'])
            end_time = int(row['t1_900'])
            tocharge = row['total_energy_Wh']
            max_power = row['maxPower']

            #Bookkeeping for the time range
            start_times += [start_time]
            end_times += [end_time]

            #Create an EV object and add it to a list
            ev = EV(self.number_EVs, start_time, end_time, tocharge, max_power)
            self.EV_list.append(ev)

        self.start_time = min(start_times)
        self.end_time = max(end_times)
        self.time_intervals = self.end_time

    def reset_profiles(self):
        for EV in self.EV_list:
            EV.bus_profile = [0 for i in range(self.end_time-self.start_time)]
            EV.cost = 0

    def best_sol(self):
        #Find the best solution with complete control.

        instance = FOCSinstance(self.data, timeStep=900, periodicity= True)

        relevant_interval = [0,0]
        relevant_interval[0] = instance.intervals_start[0]

        #Why is the -2 here? Remove it for now
        relevant_interval[1] = instance.intervals_end[-1]

        powerprofile_df = pd.DataFrame(columns = [str(i) for i in range(relevant_interval[1])])

        flowNet = FlowNet()
        flowNet.focs_instance_to_network(instance)
        flowOp = FlowOperations(flowNet.G, instance)

        focs = FOCS(instance,flowNet,flowOp)

        #Here the problem is solved
        f = focs.solve_focs()


        #data output
        column_vec = ['wf', 'flat_cost', 'cotwocost', 'peak_size']
        cost_df = pd.DataFrame(columns=column_vec)
        cost_df.loc[len(cost_df)] = ['inf', 0, 0, max([f['i'+str(i)]['t']*(1/focs.instance.tau) for i in focs.instance.I_a])]

        timelist = [0 for i in range(relevant_interval[1])]
        indx2 = relevant_interval[0]
        for indx in range(relevant_interval[1]):
            if 'i'+str(indx) in f:
                timelist[indx2] += (f['i'+str(indx)]['t'])
                indx2 += 1
        powerprofile_df.loc[len(powerprofile_df)] = timelist


        return powerprofile_df
    
    def individual_opt(self, a_t, b_t, l_t, u_t, d):
        #Choose convergion parameter, percentage of improvement
        epsilon = 0.001
        total_profile = [0 for t in a_t]
        stop = False
        iterations = 0
        while not stop:
        # while iterations < 50:
            iterations += 1

            if show_iterations:
                print(iterations)
            stop = True
            indexx = 0
            for EV in self.EV_list:
                indexx += 1
                if indexx >= 1:
                    old_cost = EV.cost
                    # print(f'old_cost is : {old_cost}')

                    for j in range(EV.intervals):
                        total_profile[j + EV.start_time] -= EV.bus_profile[j]
                    
                    EV.react_to_prices(a_t, b_t, l_t, u_t, d, total_profile)

                    for j in range(EV.intervals):
                        total_profile[j + EV.start_time] += EV.bus_profile[j]

                    new_cost = EV.cost
                    # print(f'new_cost : {new_cost}')
                    if iterations == 1:
                        stop = False
                    else:
                        if abs(old_cost - new_cost)/old_cost > epsilon:
                            stop = False

        #return total profile
        return total_profile
    


    def single_round_opt(self, a_t, b_t, l_t, u_t, d):
        #Choose convergion parameter, percentage of improvement
        epsilon = 0.001
        total_profile = [0 for t in a_t]

        self.reset_profiles()

        for EV in self.EV_list:
            EV.react_to_prices(a_t, b_t, l_t, u_t, d, total_profile)

            for j in range(EV.intervals):
                total_profile[j + EV.start_time] += EV.bus_profile[j]

        return total_profile

    def single_round_opt_time_tax(self, a_t, b_t, l_t, u_t, d, tax_factor):
        #Choose convergion parameter, percentage of improvement
        epsilon = 0.001
        total_profile = [0 for t in a_t]



        self.reset_profiles()

        for EV in self.EV_list:
            a_adjusted = [t for t in a_t]
            multiplier = 0
            for t in range(EV.start_time, EV.end_time):
                a_adjusted[t] = a_adjusted[t]*((tax_factor)**multiplier)
                multiplier += 1


            EV.react_to_prices(a_adjusted, b_t, l_t, u_t, d, total_profile)

            for j in range(EV.intervals):
                total_profile[j + EV.start_time] += EV.bus_profile[j]

        return total_profile

    
    def greedy_opt(self):
        profile = [0 for t in range(self.time_intervals)]
        for EV in self.EV_list:
            tocharge = EV.tocharge
            time = EV.start_time
            endtime = EV.end_time
            while tocharge > EV.max_power and time < endtime:
                profile[time] += EV.max_power
                tocharge -= EV.max_power
                time += 1

            if tocharge > EV.max_power:
                profile[time] += EV.max_power
            else:
                profile[time] += tocharge

        return profile




average_resuts_df = pd.DataFrame(columns=['steering','max_peak','two_norm','d', 'relative difference peak list', 'relative difference two norm list'])

if run_results:
    for d in d_values:
        print(f'd value : {d}')

        if d == save_d_value:
            total_profile_greedy = []
            total_profile_optimal = []
            total_profile_online = []
            total_profile_offline = []


        peaklist_greedy = []
        two_norm_greedy = []

        peaklist_FOCS = []
        two_norm_FOCS = []

        peaklist_individ = []
        two_norm_individ = []

        peaklist_single = []
        two_norm_single = []


        #Create a set with all dates of the year
        #15-minute time intervals
        dis = 900
        data = pd.read_csv('Data/input/202410DatasetEVOfficeParking_v0.csv', delimiter=';', encoding='utf-8')
        

        date_set = set()
        for index, row in data.iterrows():
            if row['start_datetime'][:4] == year:
                date_set.add(row['start_datetime'][:10])

        print(f'There are {len(date_set)} days in {year}')


        #Here it iterates per day, I want this to remain but use the new days instead somehow.
        for date in date_set:

            #filter based on date
            data_day = data[data['start_datetime'].str[:10] == date].copy()


            data_day['start_datetime_seconds'] = data_day['start_datetime'].apply(lambda x: int(x[-2:]) + int(x[-5:-3])*60 + int(x[-8:-6])*3600)
            data_day['end_datetime_seconds'] = data_day['end_datetime'].apply(lambda x: int(x[-2:]) + int(x[-5:-3])*60 + int(x[-8:-6])*3600)


            df = pd.DataFrame()
            df['total_energy_Wh'] = data_day['total_energy']

            #df = df.rename(columns = {'total_energy' : 'total_energy_Wh'})


            df['t0_' + str(dis)] = (data_day['start_datetime_seconds']/dis).apply(lambda x: math.floor(x))
            df['t1_' + str(dis)] = (data_day['end_datetime_seconds']/dis).apply(lambda x: math.ceil(x))

            df['maxPower'] = 22000
            df['average_power_W'] = 22000

            data_day = df


            #FOCS
            EV_case = EV_Charging_Case(data_day)
            best_sol_profile = EV_case.best_sol()

            best_sol_profile = best_sol_profile.values.tolist()[0]



            maxpeak = max(best_sol_profile)

            peaklist_FOCS.append(maxpeak)

            two_norm = 0

            for val in best_sol_profile:
                two_norm += val**2

            two_norm_FOCS.append(two_norm)


            #Find the number of time intervals
            time_intervals = EV_case.time_intervals

            #Create the cost function
            a_t = [1 for t in range(time_intervals)]
            b_t = [0 for t in range(time_intervals)]
            l_t = [0 for t in range(time_intervals)]
            u_t = [22000 for t in range(time_intervals)]


            #Price is in the form of a*x^d + b where a and b can change per interval
            flatten_individual_profile = EV_case.individual_opt(a_t, b_t, l_t, u_t, d)


            # #Convert to kWh
            # for i in range(len(flatten_individual_profile)):
            #     flatten_individual_profile [i] = flatten_individual_profile[i] / 1000

            maxpeak = max(flatten_individual_profile)

            peaklist_individ.append(maxpeak)

            two_norm = 0

            for val in flatten_individual_profile:
                two_norm += val**2

            two_norm_individ.append(two_norm)


            #print(f'Largest peak of individual optimization: {maxpeak}')



            single_round_pofile = EV_case.single_round_opt(a_t, b_t, l_t, u_t, d)

            maxpeak = max(single_round_pofile)

            peaklist_single.append(maxpeak)

            two_norm = 0

            for val in single_round_pofile:
                two_norm += val**2

            two_norm_single.append(two_norm)

            #print(f'Largest peak single-round: {maxpeak}')



            greedy_profile = EV_case.greedy_opt()

            maxpeak = max(greedy_profile)

            peaklist_greedy.append(maxpeak)

            two_norm = 0

            for val in greedy_profile:
                two_norm += val**2

            two_norm_greedy.append(two_norm)
            #print(f'Largest peak of greedy: {maxpeak}')

            if save_total_profile:
                if d == save_d_value:
                    total_profile_greedy += greedy_profile
                    total_profile_offline += flatten_individual_profile
                    total_profile_online += single_round_pofile
                    total_profile_optimal += best_sol_profile


        relative_difference_peak_individ_list = [(peaklist_individ[index] - peaklist_FOCS[index])/peaklist_FOCS[index] for index in range(len(peaklist_individ))]
        relative_difference_peak_online_list = [(peaklist_single[index] - peaklist_FOCS[index])/peaklist_FOCS[index] for index in range(len(peaklist_single))]

        relative_difference_peak_greedy_list = [(peaklist_greedy[index] - peaklist_FOCS[index])/peaklist_FOCS[index] for index in range(len(peaklist_greedy))]

        print(f'For d={d}, we get a mean offline peak difference of {np.mean(relative_difference_peak_individ_list)}')
        print(f'For d={d}, we get a mean online peak difference of {np.mean(relative_difference_peak_online_list)}')

        relative_difference_two_norm_individ_list = [(two_norm_individ[index] - two_norm_FOCS[index])/two_norm_FOCS[index] for index in range(len(peaklist_individ))]
        relative_difference_two_norm_online_list = [(two_norm_single[index] - two_norm_FOCS[index])/two_norm_FOCS[index] for index in range(len(peaklist_single))]

        relative_difference_two_norm_greedy_list = [(two_norm_greedy[index] - two_norm_FOCS[index])/two_norm_FOCS[index] for index in range(len(peaklist_greedy))]

        print(f'For d={d}, we get a mean offline two norm difference of {np.mean(relative_difference_two_norm_individ_list)}')
        print(f'For d={d}, we get a mean online two norm difference of {np.mean(relative_difference_two_norm_online_list)}')



        average_resuts_df.loc[len(average_resuts_df)] = ['FOCS',np.mean(peaklist_FOCS),np.mean(two_norm_FOCS), d, [], []]
        average_resuts_df.loc[len(average_resuts_df)] = ['individual_equilibrium',np.mean(peaklist_individ),np.mean(two_norm_individ), d, relative_difference_peak_individ_list, relative_difference_two_norm_individ_list]
        average_resuts_df.loc[len(average_resuts_df)] = ['Online',np.mean(peaklist_single),np.mean(two_norm_single), d, relative_difference_peak_online_list, relative_difference_two_norm_online_list]
        average_resuts_df.loc[len(average_resuts_df)] = ['Greedy',np.mean(peaklist_greedy),np.mean(two_norm_greedy), d, relative_difference_peak_greedy_list, relative_difference_two_norm_greedy_list]

        #print(average_resuts_df)
        
        if save_single_profile:
            if d == save_d_value:
                optimal_profile_df = pd.DataFrame(best_sol_profile, columns=['power'])
                offline_profile_df = pd.DataFrame(flatten_individual_profile, columns=['power'])
                online_profile_df = pd.DataFrame(single_round_pofile, columns=['power'])
                greedy_profile_df = pd.DataFrame(greedy_profile, columns=['power'])

                optimal_profile_df.to_csv('Data/output/optimal_profile.csv', index=False)
                online_profile_df.to_csv('Data/output/online_profile.csv', index=False)
                offline_profile_df.to_csv('Data/output/offline_profile.csv', index=False)
                greedy_profile_df.to_csv('Data/output/greedy_profile.csv', index=False)

    if save_total_profile:
        total_profile_greedy_df = pd.DataFrame(total_profile_greedy, columns=['power'])
        total_profile_offline_df = pd.DataFrame(total_profile_offline, columns=['power'])
        total_profile_online_df = pd.DataFrame(total_profile_online, columns=['power'])
        total_profile_optimal_df = pd.DataFrame(total_profile_optimal, columns=['power'])

        total_profile_greedy_df.to_csv('Data/output/total_greedy_profile.csv', index=False)
        total_profile_offline_df.to_csv('Data/output/total_offline_profile.csv', index=False)
        total_profile_online_df.to_csv('Data/output/total_online_profile.csv', index=False)
        total_profile_optimal_df.to_csv('Data/output/total_optimal_profile.csv', index=False)







if savedata:
    #We want to save the df of all date, but also make separate dataframes for the plots that we want to make to make it easier to work with in Overleaf
    average_resuts_df.to_csv('Data/output/average_results.csv', index=False)  # Set index=False to avoid writing the index to the file

    #Create a new dataframe per steering method. We probably choose the two_norm as metric but keep both for now
    offline_dataframe = pd.DataFrame(columns=['d_value','max_peak','two_norm'])
    online_dataframe = pd.DataFrame(columns=['d_value','max_peak','two_norm'])

    #Fill the dataframes
    for index, row in average_resuts_df.iterrows():
        if row['steering'] == 'individual_equilibrium':
            offline_dataframe.loc[len(offline_dataframe)] = [row['d'],row['max_peak'],row['two_norm']]
        else:
            if row['steering'] == 'Online':
                online_dataframe.loc[len(online_dataframe)] = [row['d'],row['max_peak'],row['two_norm']]

    #Save the dataframes
    online_dataframe.to_csv('Data/output/online_results.csv', index=False)
    offline_dataframe.to_csv('Data/output/offline_results.csv', index=False)





if plot_result:
    x = [i for i in range(len(flatten_individual_profile))]
    # Plot the two lists
    plt.plot(x, best_sol_profile, label="FOCS", color='blue', marker='o')  # Line 1
    plt.plot(x, flatten_individual_profile, label="Price Steering", color='red', marker='s')  # Line 2
    plt.plot(x, greedy_profile, label="Greedy", color='green', marker='s')  # Line 2
    plt.plot(x, single_round_pofile, label="Single-round", color='yellow', marker='s')  # Line 2

    # Add labels, legend, and title
    plt.xlabel("Time inteval")
    plt.ylabel("kW")
    plt.title("Energy profiles for different strategies")
    plt.legend()

    # Show the graph
    plt.show()

if run_and_save_time_tax:
    online_taxed_results = pd.DataFrame(columns=['tax','max_peak','two_norm'])

    for tax in time_taxes:
        peaklist_single = []
        two_norm_single = []
        total_profile_online = []

        for date in date_range:

            #15-minute time intervals
            dis = 900
            data = pd.read_csv('Data/input/202410DatasetEVOfficeParking_v0.csv', delimiter=';', encoding='utf-8')

            #filter based on date
            data = data[data['start_datetime'].str[:10] == date]


            data['start_datetime_seconds'] = data['start_datetime'].apply(lambda x: int(x[-2:]) + int(x[-5:-3])*60 + int(x[-8:-6])*3600)
            data['end_datetime_seconds'] = data['end_datetime'].apply(lambda x: int(x[-2:]) + int(x[-5:-3])*60 + int(x[-8:-6])*3600)


            df = pd.DataFrame()
            df['total_energy_Wh'] = data['total_energy']

            #df = df.rename(columns = {'total_energy' : 'total_energy_Wh'})


            df['t0_' + str(dis)] = (data['start_datetime_seconds']/dis).apply(lambda x: math.floor(x))
            df['t1_' + str(dis)] = (data['end_datetime_seconds']/dis).apply(lambda x: math.ceil(x))

            df['maxPower'] = 22000
            df['average_power_W'] = 22000

            data = df

            EV_case = EV_Charging_Case(data)


            #Find the number of time intervals
            time_intervals = EV_case.time_intervals

            #Create the cost function
            a_t = [1 for t in range(time_intervals)]
            b_t = [0 for t in range(time_intervals)]
            l_t = [0 for t in range(time_intervals)]
            u_t = [22000 for t in range(time_intervals)]



            single_round_pofile = EV_case.single_round_opt_time_tax(a_t, b_t, l_t, u_t, time_tax_d, tax)

            maxpeak = max(single_round_pofile)

            peaklist_single.append(maxpeak)

            two_norm = 0

            for val in single_round_pofile:
                two_norm += val**2

            two_norm_single.append(two_norm)

            #print(f'Largest peak single-round: {maxpeak}')
            total_profile_online += single_round_pofile
        #Save the total profiles
        total_profile_online_df = pd.DataFrame(total_profile_online, columns=['power'])
        total_profile_online_df.to_csv(f'Data/output/total_online_profile_tax_{tax}.csv', index=False)

if investigate_single_instance:
    EV_case.reset_profiles()
    flatten_individual_profile = EV_case.individual_opt(a_t, b_t, l_t, u_t, d)


    #Isolate a single EV that has load in the peak and that can shift it to somewhere else (?)
    isolate_time_intervals_list = [36,37,38]
    peak1 = flatten_individual_profile[36:39]
    peak2 = flatten_individual_profile[55:61]
    dip = flatten_individual_profile[28:32]
    baseline = best_sol_profile[33]

    #The plan is to take some EV that can decrease in the peak and increase in the dip and see why this does not happen

    #Step 1: Find such an EV, meaning it is there from 28 until 38.
    #Notes: EV 1 and 2 can not help, EV 7 can switch stuff from busy (37) interval to not so busy (34). Why does it not do this? EV 9 can also do this. EV 13 puts everything is busy interval 38. Why?? EV 24 puts almost everythin in 38.
    indexx = 0
    for EV in EV_case.EV_list:
        indexx += 1
        #if EV.start_time < 38:
        if indexx == 13:
            print(indexx)
            print(f'Start time: {EV.start_time}, end time: {EV.end_time}, to charge: {EV.tocharge}, cost: {EV.cost}')

            profile_of_others = copy.deepcopy(flatten_individual_profile)
            for index in range(len(EV.bus_profile)):
                profile_of_others[index + EV.start_time] -= EV.bus_profile[index]

            print(profile_of_others[EV.start_time:EV.end_time])

            new_cost = 0
            for i in range(len(EV.bus_profile)):
                 new_cost += EV.bus_profile[i]*(profile_of_others[i + EV.start_time] + EV.bus_profile[i] )**d
            print(f'Old cost should be: {new_cost}')
                
            #EV.react_to_prices(a_t, b_t, l_t, u_t, d, profile_of_others)
            print(EV.bus_profile)

            EV.bus_profile = [0.96, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0]

            #Manually adjust profile
            print(EV.bus_profile)


            #print(f'New cost: {EV.cost}')

            #old_profile = [flatten_individual_profile[i] - profile_of_others[i] for i in range(len(profile_of_others))]

            # print(old_profile)

            new_cost = 0
            for i in range(len(EV.bus_profile)):
                 new_cost += EV.bus_profile[i]*(profile_of_others[i + EV.start_time] + EV.bus_profile[i] )**d

            print(f'New cost: {new_cost}')

            for index in range(len(EV.bus_profile)):
                profile_of_others[index + EV.start_time] += EV.bus_profile[index]

            # for i in range(len(flatten_individual_profile)):
            #     print(flatten_individual_profile[i] - profile_of_others[i])


            x = [i for i in range(len(flatten_individual_profile))]
            # Plot the two lists
            plt.plot(x, flatten_individual_profile, label="Old", color='blue', marker='o')  # Line 1
            plt.plot(x, profile_of_others, label="New", color='red', marker='s')  # Line 2
            plt.plot(x, best_sol_profile, label="FOCS", color='yellow', marker='s')

            # Add labels, legend, and title
            plt.xlabel("Time inteval")
            plt.ylabel("kW")
            plt.title("Energy profiles for different strategies")
            plt.legend()

            # Show the graph
            plt.show()

            # Something is changed, check whether the obj value is improved more than threshold. Forcing more iterations does not improve result
