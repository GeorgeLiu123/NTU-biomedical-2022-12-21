import matplotlib.pyplot as plt

####################################################################################
# STEP1: Definite the functions
def f(step, subtrate, intermediate):   # The derivative function of S
    a1 = 100/60
    a2 = 600/60
    lamda = a1 * 1
    fun_s = -lamda*subtrate + (a1*subtrate + a2)*intermediate
    return fun_s

def g(step, subtrate, intermediate):   # The derivative function of C
    b1 = 100/60
    b2 = 600/60
    b3 = 150/60
    theta = b2 + b3
    lamda = b1 * 1
    fun_c = lamda*subtrate - (theta + b1*subtrate)*intermediate
    return fun_c

#######################################################################################
# STEP2: Input the initial value
E0 = 1         # the initial value
k3 = 150/60    # the initial value
h = 0.001      # step width
n = int(30/0.001)  # The # of steps in 30s
t = []         # The list of steps
temp = 0
c = [0]
s = [10]
for i in range(1, n+1):
    temp = h*i
    print('temp:', temp)
    t.append(temp)  # add initial value

#######################################################################################
# STEP3: Model simulation
for k in range(0, n-1):
    z0 = h * f(t[k], s[k], c[k])
    l0 = h * g(t[k], s[k], c[k])
    z1 = h * f(t[k]+h/2, s[k]+z0/2, c[k]+l0/2)
    l1 = h * g(t[k]+h/2, s[k]+z0/2, c[k]+l0/2)
    z2 = h * f(t[k]+h/2, s[k]+z1/2, c[k]+l1/2)
    l2 = h * g(t[k]+h/2, s[k]+z1/2, c[k]+l1/2)
    z3 = h * f(t[k]+h, s[k]+z2, c[k]+l2)
    l3 = h * g(t[k]+h, s[k]+z2, c[k]+l2)
    s_temp = s[k] + (z0+2*z1+2*z2+z3)/6
    c_temp = c[k] + (l0+2*l1+2*l2+l3)/6
    s.append(s_temp)   # add new s to list
    c.append(c_temp)   # add new c to list

########################################################################################
# STEP4: output the solutions
E = []
P = []
Vp = []

# The solution of E
for item in c:
    E_value = E0 - item
    E.append(E_value)

# The solution to P
the_sum = 0
for k in range(0, n):
    the_sum = c[k] + the_sum
    P_value = k3 * the_sum * h  # integrating discrete results
    P.append(P_value)

# The solution to the change rate of P
for item in c:
    Vp_value = k3*item
    Vp.append(Vp_value)


############################################################################################
# figure 1, plot the relationship between time (x-axis) and concentration (y-axis)
plt.figure(dpi=200)
plt.axis([0, 35, 0, 15])
plt.xlabel('Time(s)')
plt.ylabel('Concentration(uM)')
l_1, = plt.plot(t, E, color='green')
l_2, = plt.plot(t, s, color='orange')
l_3, = plt.plot(t, c, color='red')
l_4, = plt.plot(t, P, color='blue')
plt.legend(loc='best', handles=[l_1, l_2, l_3, l_4], fontsize=6,
           labels=['E(enzyme)', 'S(substrate)', 'ES(compound)', 'P(product)'])
plt.show()


# figure 2, plot the relationship between concentration of S (x-axis) and change rate of P (y-axis)
plt.figure(dpi=200)
plt.axis([0, 15, 0, 2])
plt.xlabel('Concentration of S(uM)')
plt.ylabel('Change rate of P(uM/s)')
l_5, = plt.plot(s, Vp, color='red')
plt.show()

