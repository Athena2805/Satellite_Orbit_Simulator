from vpython import*
#Note1:vpython library integrates Visual and python lang. to create 3D graphics. Basically 3D models which can be animated here we import sphere , vector, graph , mag, gcurve etc)

G = 6.67e-11
m1 = 4e30
m2 = 2e30
rdist = 2e10
M = m1 + m2

x1 = -(m2/M) * rdist
x2 = (m1/M) * rdist

star1 = sphere(pos=vector(x1,0,0), radius=2e9,color=color.yellow, make_trail=True)

star2 = sphere(pos=vector(x2,0,0), radius=8e8,color=color.cyan, make_trail=True)

Rcom=(star1.pos*m1 + star2.pos*m2)/M
r=star2.pos-star1.pos
v1circle=sqrt(G*m2*mag(star1.pos)/mag(r)**2)

star1.v=vector(0,0.5*v1circle,0)
# Note2 : when the y element of vector star1.v is 1*v1circle we get CONCENTRIC ORBITS but since her there is 0.5*v1circle we get two orbits which cross each other
star1.p=m1*star1.v
star2.p=-star1.p

mu=m1*m2/M
l=mag(cross(star1.pos,star1.p)+cross(star2.pos,star2.p))

t=0
dt=100
#Increasing the dt crashed the model rather gave it straight line at start or end or extreme points of the potetial between m1 and m2

tgraph=graph(xtitle="r [m]",ytitle="U [J]")
fU=gcurve(color=color.blue, dot=True)

while t<1e6:
# This line is to Keep the window alive, how many calculaion per second
# Doubt if this is where rendering was causing problem previosly
    rate(100)
  #rate of 1000 created straigh line motion of the bodies in the start or end of the simulation and varyingly at some extreme points of potential between the two bodies
    r=star2.pos-star1.pos
    F2=-G*m1*m2*norm(r)/mag(r)**2
    star2.p=star2.p+F2*dt
    star1.p=star1.p-F2*dt
    star1.pos=star1.pos+star1.p*dt/m1
    star2.pos=star2.pos+star2.p*dt/m2
    Ug=-G*m1*m2/mag(r)
    Uc=l**2/(2*mu*mag(r)**2)
    Uef=Ug+Uc3 #effective potential
    fU.plot(mag(r),Uef) #Plot of the effective potential (Uef(J)) between the two bodies and dist between the two bodies (r(m)) or (r=star2.pos-star1.pos).
    t=t+dt
