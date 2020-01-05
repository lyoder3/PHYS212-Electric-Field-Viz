from vpython import *
# some constants defined here
k = 9e9
charge_e = 1.60217662e-19
mass_e = 9.10938356e-31


# Charge class is what I have developed for storing point charges
# It features some useful calculations and makes it easy to repeat them
# It is also nice because we can easily create many point charge objects
class Charge:
    """init method is what is called when you declare a new object of a class.
    It takes the different object attributes as arguments along with 'self'
    self is important because it allows us to set the arguments passed to the attributes
    to the current object. 
    
    charge can be any number positive or negative and can have decimals
    position must be passed as a vpython vector() of the form vector(x,y,z) where (x,y,z)
    is an ordered pair representing the charge's position in space"""

    # Our charge object currently features two attributes: charge and position
    def __init__(self, charge, position):
        self.charge = charge
        self.position = position

    # This function just finds the relative position vector from the charge to another point 'v'
    # V must be passed as a vpython vector() object
    def delta_r(self, v):
        return v - self.position

    """The function below calculates the electric field due to the charge at "point"
     Once again, point must be passed as a vpython vector object"""

    def efield(self, point):
        return k * (self.charge) / (self.delta_r(point).mag ** 3) * self.delta_r(point)
def runge_kutta_2(charges, r_0, ds):
    net_efield = vector(0,0,0)
    for charge in charges:
        net_efield += charge.efield(r_0)
    k1 = ds * hat(net_efield)
    net_efield = vector(0,0,0)
    for charge in charges:
        net_efield+= charge.efield(r_0 + k1/2)
    k2 = ds * hat(net_efield)
    net_efield = vector(0,0,0)
    for charge in charges:
        net_efield+= charge.efield(r_0 + k2/2)
    k3 = ds * hat(net_efield)
    net_efield = vector(0,0,0)
    for charge in charges:
        net_efield+= charge.efield(r_0 + k3)
    k4 = ds * hat(net_efield)
    return r_0 + ((1/6) * (k1 + 2*k2 + 2*k3 + k4))
def euler(charges, r_0, ds):
    net_efield = vector(0,0,0)
    for charge in charges:
        net_efield+=charge.efield(r_0)
    return r_0 + hat(net_efield)*ds


