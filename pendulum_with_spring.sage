#
# Deriving the state space equations for when we put in a spring to our
# inverted pendulum model
#
# d denotes dot
# dd denotes double dot
var('g Mc bc Ks rd x xd xdd V')
var('R Jd bm Km rd thetam thetamd thetamdd')
var('mp1 l1 J1 theta1 theta1d theta1dd')
var('mp2 l2 J2 theta2 theta2d theta2dd')

soln = solve([
    (mp1*l1^2+J1)*theta1dd - mp1*g*l1*theta1 == mp1*l1*xdd,  # (14)
    (mp2*l2^2+J2)*theta2dd + mp2*g*l2*theta2 == -mp2*l2*xdd, # (15)
    (Mc+mp1+mp2)*xdd + bc*xd - Ks*x + Ks*rd*thetam - mp1*l1*theta1dd + mp2*l2*theta2dd == 0, # (16)
    #Jd*thetamdd + bm*thetamd - Ks*rd*thetam + Ks*x == 0, # (17)
    #V == Km*thetamd + R*Ks/Km*rd*thetam - R*Ks/Km*x # (18)
    Jd*thetamdd + (bm+Km^2/R)*thetamd - Km/R*V == 0 # (20)
], xdd, theta1dd, theta2dd, thetamdd)
print(soln)
