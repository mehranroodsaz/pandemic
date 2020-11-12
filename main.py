import numpy as np

#  number of nodes
nodes_no= 3

# number of days you want to run
day_no = 2

# number of the timeslos in a day.
# e.g. if the time slot is 3 the day will be divided in 8 time slots
timeslot_no=4

# s : Susceptible (S)
# i: Infectious (I)
# r: Recovered (R)
# sp : susceptible after the infection reaction
# ip : infectious after the infection reaction
# rp : recovered after the infection reaction
s_index = 0
i_index = 1
r_index = 2
sp_index = 3
ip_index = 4
rp_index = 5


# read the initial values for Susceptible (S), Infectious (I), and Recovered (R).
initial_node = np.loadtxt("nodes2_input.txt")

# Load the values for F. F : expected volume of human mobility from one node to another during a time interval
f =  np.loadtxt("F2.txt")

# F forward factor is per node
f_forward_factor_per_node = [0,0.1,0.3,0.2,0.2,0.1,0.05,0.05]

# F return factor is per node.
f_return_factor_per_node = [0,0.3,0.05,0.1,0.2,0.4,0.2,0.05]


# F forward factor is per node
f_forward_factor_per_time_slot = [0,0.1,0.3,0.2]

# F return factor is per node.
f_return_factor_per_time_slot = [0,0,0.05,0.1]

# b seems to be per node but at this moment all nodes have the value 0.5
b=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]

# constant
v=0.1

# initial nodes matrix. This matrix contains all nodes with all information (s,i,r,sp,ip,rp)
# in the flow of this program values will be filled in
nodes = [[[0,0,0,0,0,0] for i in range(timeslot_no*day_no) ]for j in range(nodes_no) ]

# set the initial values for S, I and R for all nodes for the first time slot.
for node_index in range(nodes_no):
     s = initial_node[node_index][s_index]
     i = initial_node[node_index][i_index]
     r = initial_node[node_index][r_index]
     p = s+i+r
     nodes[node_index][0][s_index] = s
     nodes[node_index][0][i_index] = i
     nodes[node_index][0][r_index] = r
     nodes[node_index][0][sp_index] = s-b[node_index]*(s*i/p)
     nodes[node_index][0][ip_index] = i+b[node_index]*s*i/p-v*i
     nodes[node_index][0][rp_index] = r+v*i

# np.asarray(nodes).tofile('yourfile.txt',sep=" ",format="%s")

# p_var_index could be sp_index, ip_index or rp_index
def calc_value(current_time_index, time_slot, current_node_index, p_var_index):
     sigma_forward_enter = 0
     sigma_forward_exit = 0
     sigma_return_enter = 0
     sigma_return_exit = 0

     prime_n = nodes[current_node_index][current_time_index-1][p_var_index]

     for node_index in range(nodes_no):
          # calculat section sigma_forward
          forward_insolation_level_current_node = f_forward_factor_per_node[current_node_index]
          forward_insolation_level_node_index = f_forward_factor_per_node[node_index]

          return_insolation_level_current_node = f_return_factor_per_node[current_node_index]
          return_insolation_level_node_index = f_return_factor_per_node[node_index]

          f_value_column = f[node_index][current_node_index]
          f_value_row = f[current_node_index][node_index]

          prime_m = nodes[node_index][current_time_index-1][p_var_index]
          s_m = nodes[node_index][current_time_index-1][s_index]
          i_m = nodes[node_index][current_time_index-1][i_index]
          r_m = nodes[node_index][current_time_index-1][r_index]
          p_m = s_m + i_m + r_m

          s_n = nodes[current_node_index][current_time_index-1][s_index]
          i_n = nodes[current_node_index][current_time_index-1][i_index]
          r_n = nodes[current_node_index][current_time_index-1][r_index]
          p_n = s_n + i_n + r_n

          forward_value_column = f_value_column * f_forward_factor_per_time_slot[time_slot-1]
          forward_value_row = f_value_row * f_forward_factor_per_time_slot[time_slot-1]

          return_value_column = f_value_column * f_return_factor_per_time_slot[time_slot-1]
          return_value_row = f_value_row * f_return_factor_per_time_slot[time_slot-1]


          sigma_forward_enter += forward_value_column * (1-forward_insolation_level_node_index)
          sigma_forward_exit += forward_value_row * prime_m / p_m * (1-forward_insolation_level_current_node)

          sigma_return_enter += return_value_column * prime_m / p_m *(1- return_insolation_level_node_index)
          sigma_return_exit  += return_value_row * prime_n / p_n * (1- return_insolation_level_current_node)


     result = prime_n - sigma_forward_enter * prime_n / p_n + sigma_forward_exit + sigma_return_enter - sigma_return_exit
     return result

for day in range(day_no):
     for timeslot in range(timeslot_no):
          # index of the time slot in a specific day
          # e.g. assume we have 10 days and each has been devided into 8 time slots. the 2nd time slot of the 3rd day =  1 + 2 * 8 = 17
          timeslot_index = timeslot + day *  timeslot_no
          for node_index in range(nodes_no):
               # skeep the first row
               if not (day == 0 and timeslot_index == 0 ) :
                    s = calc_value(timeslot_index, timeslot, node_index,sp_index)
                    i = calc_value(timeslot_index, timeslot, node_index,ip_index)
                    r = calc_value(timeslot_index, timeslot, node_index,rp_index)
                    p = s + i + r
                    nodes[node_index][timeslot_index][s_index] = s
                    nodes[node_index][timeslot_index][i_index] = i
                    nodes[node_index][timeslot_index][r_index] = r

                    # prime values
                    nodes[node_index][timeslot_index][sp_index] = s-b[node_index]*(s*i/p)
                    nodes[node_index][timeslot_index][ip_index] = i+b[node_index]*(s*i/p)-v*i
                    nodes[node_index][timeslot_index][rp_index] = r+v*i

for node_index in range(nodes_no):
     print ("Node "+str(node_index+1)+"\n")
     for time_day_index in range(day_no * timeslot_no):
         print(*nodes[node_index][time_day_index], sep=' | ')
     print("\n\n")


